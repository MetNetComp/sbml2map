"""Node avoidance prefers a direct U, with checked endpoint-step fallbacks."""
import contextlib
import copy
import io
import math
import unittest
from parallel_routing import avoid_nodes, _legs, validate_route_offset


def path(vertical=False, reverse=False, **kwargs):
    points=[(0,0),(3,0)]
    if vertical:points=[(y,x) for x,y in points]
    if reverse:points.reverse()
    return dict(rid='r',priority=2,backbone=True,points=points,**kwargs)


class NodeRoutingTests(unittest.TestCase):
    def test_rotated_u_shape_keeps_endpoints_and_moves_backbone_edges(self):
        for vertical in (False,True):
            for reverse in (False,True):
                with self.subTest(vertical=vertical,reverse=reverse):
                    p=path(vertical,reverse);before=copy.deepcopy(p)
                    node=(0,1) if vertical else (1,0)
                    result=avoid_nodes([p],obstacles=[node])
                    self.assertEqual(result['shifted_segments'],1)
                    self.assertEqual(p['points'][0],before['points'][0])
                    self.assertEqual(p['points'][-1],before['points'][-1])
                    self.assertEqual(len(p['points']),4)
                    legs=list(_legs(p['points']))
                    self.assertEqual(len(legs),3)
                    self.assertAlmostEqual(legs[0]['hi']-legs[0]['lo'],.3)
                    self.assertAlmostEqual(legs[-1]['hi']-legs[-1]['lo'],.3)
                    self.assertEqual(legs[0]['axis'],legs[-1]['axis'])
                    self.assertNotEqual(legs[0]['axis'],legs[1]['axis'])

    def test_offset_value_and_zero(self):
        for offset in (.2,.3,.5,0):
            with self.subTest(offset=offset):
                p=path();original=copy.deepcopy(p)
                report=avoid_nodes([p],obstacles=[(1,0)],route_offset=offset)
                if offset:
                    self.assertEqual(report['shifted_segments'],1)
                    self.assertAlmostEqual(abs(p['points'][1][1]),offset)
                else:self.assertEqual(p,original)

    def test_endpoints_gray_lines_and_noncolliding_lines_stay_put(self):
        for obstacles,priority in [([(0,0),(3,0)],2), ([(1,.5)],2), ([(1,0)],0)]:
            p=path();p['priority']=priority;before=copy.deepcopy(p)
            self.assertEqual(avoid_nodes([p],obstacles=obstacles)['shifted_segments'],0)
            self.assertEqual(p,before)

    def test_both_connector_legs_are_checked(self):
        p=path();before=copy.deepcopy(p)
        report=avoid_nodes([p],obstacles=[(1,0),(0,.15),(0,-.15),(.15,.15),(.15,-.15)])
        self.assertEqual(p,before)
        self.assertEqual(report['shifted_segments'],0)
        self.assertEqual(len(report['unresolved']),1)

    def test_other_line_selects_opposite_side(self):
        p=path();other=dict(rid='gray',priority=0,points=[(0,.3),(3,.3)])
        report=avoid_nodes([p,other],obstacles=[(1,0)])
        self.assertEqual(report['shifted_segments'],1)
        self.assertAlmostEqual(p['points'][1][1],-.3)

    def test_shared_trunks_move_together_for_multiple_obstacles(self):
        a=path();b=path();b['points'][-1]=(4,0)
        report=avoid_nodes([a,b],obstacles=[(1,0),(2,0)])
        self.assertEqual(report['shifted_segments'],2)
        self.assertEqual(a['points'][1][1],b['points'][1][1])
        self.assertEqual(len(a['points']),4)
        self.assertEqual(len(b['points']),4)

    def test_scale_invariance_determinism_and_idempotency(self):
        a=path();b=path();b['points']=[(x*180,y*180) for x,y in b['points']]
        avoid_nodes([a],obstacles=[(1,0)])
        avoid_nodes([b],obstacles=[(180,0)],grid=180)
        for p,q in zip(a['points'],b['points']):
            for x,y in zip(p,q):self.assertAlmostEqual(x,y/180)
        before=copy.deepcopy(a);avoid_nodes([a],obstacles=[(1,0)]);self.assertEqual(a,before)
        duplicate=path();avoid_nodes([duplicate],obstacles=[(1,0)]);self.assertEqual(a,duplicate)

    def test_invalid_offsets(self):
        for value in (-1,math.nan,math.inf,-math.inf):
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'route-offset'):
                validate_route_offset(value)

    def test_cli_default_custom_zero_and_invalid_offsets(self):
        from sbml2map import parser
        p=parser()
        self.assertEqual(p.parse_args([]).route_offset,.3)
        self.assertEqual(p.parse_args(['--route-offset','.5']).route_offset,.5)
        self.assertEqual(p.parse_args(['--route-offset','0']).route_offset,0)
        for value in ('-1','nan','inf','-inf'):
            with self.subTest(value=value),contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):
                p.parse_args(['--route-offset='+value])


class EscherNodeRoutingTests(unittest.TestCase):
    def test_custom_offset_is_persisted_and_zero_or_fba_disables_node_detours(self):
        from export_escher import create_map, validate_map_geometry
        from escher_routing import separate_routes
        from reaction_network import build_network
        model=dict(id='node_test',metabolites=[dict(id=x) for x in ('a','b')],reactions=[
            dict(id='red',metabolites={'a':1},lower_bound=0,upper_bound=10),
            dict(id='other',metabolites={'b':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'rxn:red':[0,0],'met:a':[4,0],'rxn:other':[2,0],'met:b':[2,1]}
        original=create_map(model,nodes,edges,xy,scale=180,padding=40)
        orange_pairs={frozenset(('rxn:other','met:b'))}
        for offset in (.3,.5,0):
            with self.subTest(offset=offset):
                data=copy.deepcopy(original)
                separate_routes(data,{'red'},orange_pairs,route_offset=offset)
                meta=data[0]['metadata']
                report=meta['overlap_routing']['node_detours']
                self.assertEqual(meta['route_offset_grid'],offset)
                self.assertEqual(report['route_offset_grid'],offset)
                self.assertEqual(report['shifted_segments'],1 if offset else 0)
                self.assertTrue(validate_map_geometry(data,model,nodes,xy)['orthogonal_segments_verified'])
                first=copy.deepcopy(data);separate_routes(data,{'red'},orange_pairs);self.assertEqual(data,first)
                separate_routes(data,enabled=False)
                self.assertIsNone(data[0]['metadata']['overlap_routing']['node_detours'])
                self.assertTrue(validate_map_geometry(data,model,nodes,xy)['orthogonal_segments_verified'])


class EndpointStepTests(unittest.TestCase):
    def test_formate_configuration_in_both_directions_and_rotations(self):
        for vertical in (False,True):
            for reverse in (False,True):
                with self.subTest(vertical=vertical,reverse=reverse):
                    transform=lambda p:(p[1],p[0]) if vertical else p
                    target=dict(rid='EX_for_e',priority=1,points=[(6,2),(6,0)])
                    if reverse:target['points'].reverse()
                    others=[dict(rid='FORt',priority=1,points=[(5,2),(6,2)]),
                            dict(rid='SUCCt2_2',priority=1,points=[(8,2),(6,2)])]
                    for p in [target,*others]:p['points']=[transform(x) for x in p['points']]
                    before=copy.deepcopy(target)
                    report=avoid_nodes([target,*others],obstacles=[transform(p) for p in [(6,1),(6,2),(6,0),(5,2),(8,2)]])
                    self.assertFalse(report['unresolved'])
                    self.assertEqual(report['endpoint_step_fallback_segments'],1)
                    self.assertEqual(len(target['points']),5)
                    self.assertEqual(target['points'][0],before['points'][0])
                    self.assertEqual(target['points'][-1],before['points'][-1])
                    self.assertEqual(len(list(_legs(target['points']))),4)
                    self.assertTrue(any(abs(x[0 if vertical else 1]-(2-.15))<1e-6 for x in target['points']))

    def test_both_endpoints_can_step_when_both_connectors_are_blocked(self):
        p=path()
        others=[dict(rid='left',priority=1,points=[(0,-1),(0,1)]),
                dict(rid='right',priority=1,points=[(3,-1),(3,1)])]
        report=avoid_nodes([p,*others],obstacles=[(1,0)])
        self.assertEqual(len(p['points']),6)
        self.assertEqual(report['endpoint_step_fallback_segments'],1)
        detail=report['detours'][0]
        self.assertTrue(detail['inset_start'] and detail['inset_end'])
        self.assertAlmostEqual(p['points'][1][0],.15)
        self.assertAlmostEqual(p['points'][-2][0],2.85)

    def test_fallback_still_checks_new_connector_nodes_and_offset(self):
        p=path()
        others=[dict(rid='near',priority=1,points=[(0,-1),(0,1)])]
        report=avoid_nodes([p,*others],obstacles=[(1,0),(.15,.3)],route_offset=.3)
        self.assertEqual(report['endpoint_step_fallback_segments'],1)
        self.assertAlmostEqual(p['points'][2][1],-.3)
        self.assertAlmostEqual(abs(report['detours'][0]['offset_grid']),.3)


class FormateArtifactRegressionTests(unittest.TestCase):
    def test_saved_escher_and_svg_formate_branch_clear_succinate(self):
        import csv
        import json
        import re
        from pathlib import Path
        import xml.etree.ElementTree as ET
        from parallel_routing import _distance
        folder=Path(__file__).resolve().parent/'results/e_coli_core'
        data=json.loads((folder/'escher_map.json').read_text())
        nodes=data[1]['nodes'];scale=data[0]['metadata']['coordinate_transform']['scale']
        succ=next(n for n in nodes.values() if n.get('bigg_id')=='succ_c' and n.get('node_is_primary'))
        reaction=next(r for r in data[1]['reactions'].values() if r['bigg_id']=='EX_for_e')
        for segment in reaction['segments'].values():
            points=[tuple(nodes[segment[k]][axis] for axis in ('x','y')) for k in ('from_node_id','to_node_id')]
            for leg in _legs(points):self.assertGreaterEqual(_distance((succ['x'],succ['y']),leg),.29*scale)
        with (folder/'coordinates.csv').open() as f:coords=list(csv.DictReader(f))
        with (folder/'edges.csv').open() as f:edges=list(csv.DictReader(f))
        svg=ET.parse(folder/'layout.svg')
        def position(mid):
            index=next(i for i,c in enumerate(coords) if c['id']=='met:'+mid)
            node=svg.find(f'.//*[@id="node-{index}"]')
            use=node.find('.//{*}use')
            return float(use.get('x')),float(use.get('y'))
        succ_pos=position('succ_c');for_pos=position('for_e')
        step=abs(succ_pos[1]-for_pos[1])
        index=next(i for i,e in enumerate(edges) if e['reaction_id']=='EX_for_e')
        d=svg.find(f'.//*[@id="edge-{index}"]').find('.//{*}path').get('d')
        numbers=list(map(float,re.findall(r'[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?',d)))
        points=list(zip(numbers[::2],numbers[1::2]))
        self.assertGreaterEqual(len(points),5)
        for leg in _legs(points):self.assertGreaterEqual(_distance(succ_pos,leg),.29*step)
