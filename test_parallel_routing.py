"""Rectangular detours preserve endpoints, backbones and same-reaction trunks."""
import copy
import unittest
from parallel_routing import separate_parallel, _legs


def path(rid, y=0, length=3, priority=1, vertical=False, backbone=False):
    points=[(0,y),(length,y)]
    if vertical:points=[(b,a) for a,b in points]
    return dict(rid=rid,priority=priority,points=points,backbone=backbone)


class ParallelRoutingTests(unittest.TestCase):
    def test_horizontal_and_vertical_u_detours(self):
        for vertical in (False,True):
            with self.subTest(vertical=vertical):
                paths=[path('red',priority=2,vertical=vertical,backbone=True),path('orange',y=.1,vertical=vertical)]
                before=copy.deepcopy(paths)
                report=separate_parallel(paths)
                self.assertEqual(paths[0],before[0])
                self.assertEqual(report['shifted_segments'],1)
                self.assertFalse(report['unresolved'])
                self.assertEqual(paths[1]['points'][0],before[1]['points'][0])
                self.assertEqual(paths[1]['points'][-1],before[1]['points'][-1])
                legs=list(_legs(paths[1]['points']))
                shaft=max(legs,key=lambda x:x['hi']-x['lo'])
                self.assertAlmostEqual(abs(shaft['level']),.3)
                self.assertEqual(len(legs),5)

    def test_one_grid_inclusive_and_short_or_distant_paths_unchanged(self):
        for length,y,shifted in [(1,.1,1),(.99,.1,0),(3,.3,0),(3,.31,0)]:
            with self.subTest(length=length,y=y):
                paths=[path('a',priority=2),path('b',length=length,y=y)]
                self.assertEqual(separate_parallel(paths)['shifted_segments'],shifted)

    def test_only_highlighted_different_reactions_separate(self):
        for rid,priority in [('red',1),('gray',0)]:
            paths=[path('red',priority=2),path(rid,priority=priority)]
            self.assertEqual(separate_parallel(paths)['shifted_segments'],0)

    def test_shared_trunks_take_same_lane(self):
        paths=[path('red',priority=2),path('orange'),path('orange',length=4)]
        separate_parallel(paths)
        levels=[max(_legs(p['points']),key=lambda x:x['hi']-x['lo'])['level'] for p in paths[1:]]
        self.assertEqual(levels[0],levels[1])

    def test_three_reactions_get_distinct_lanes(self):
        paths=[path('red',priority=2),path('a'),path('b')]
        report=separate_parallel(paths)
        self.assertFalse(report['unresolved'])
        levels=[max(_legs(p['points']),key=lambda x:x['hi']-x['lo'])['level'] for p in paths]
        for i,a in enumerate(levels):
            for b in levels[i+1:]:self.assertGreaterEqual(abs(a-b),.3-1e-7)

    def test_obstacle_selects_other_side_and_both_blocked_are_reported(self):
        for obstacles,shifted in [([(1,.3)],1), ([(1,.3),(1,-.3)],0)]:
            paths=[path('red',priority=2),path('orange')]
            report=separate_parallel(paths,obstacles=obstacles)
            self.assertEqual(report['shifted_segments'],shifted)
            if shifted:self.assertAlmostEqual(paths[1]['points'][2][1],-.3)
            else:self.assertEqual(report['unresolved'][0]['reaction_id'],'orange')

    def test_fixed_backbones_are_never_moved(self):
        paths=[path('a',priority=2,backbone=True),path('b',priority=2,backbone=True)]
        before=copy.deepcopy(paths);report=separate_parallel(paths)
        self.assertEqual(paths,before)
        self.assertEqual(report['unresolved'][0]['reason'],'fixed backbone')

    def test_deterministic_scale_invariant_and_orthogonal(self):
        paths=[path('red',priority=2),path('orange',y=.1)]
        other=copy.deepcopy(paths)
        for p in other:p['points']=[(x*180,y*180) for x,y in p['points']]
        separate_parallel(paths);separate_parallel(other,grid=180)
        for a,b in zip(paths,other):
            self.assertEqual(len(a['points']),len(b['points']))
            for p,q in zip(a['points'],b['points']):
                for x,y in zip(p,q):self.assertAlmostEqual(x,y/180)
        duplicate=[path('red',priority=2),path('orange',y=.1)]
        separate_parallel(duplicate);self.assertEqual(paths,duplicate)


class EscherParallelIntegrationTests(unittest.TestCase):
    def test_routing_is_connected_deterministic_and_can_restore_neutral_fba(self):
        from export_escher import create_map, validate_map_geometry
        from escher_routing import separate_routes
        from reaction_network import build_network
        model=dict(id='lanes',metabolites=[dict(id=x) for x in ('a','b')],reactions=[
            dict(id='red',metabolites={'a':1},lower_bound=0,upper_bound=10),
            dict(id='orange',metabolites={'b':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={key:([0,0] if key.startswith('rxn:') else [4,0]) for key in nodes}
        data=create_map(model,nodes,edges,xy,scale=180,padding=40)
        original=copy.deepcopy(data)
        pairs={frozenset(('rxn:orange','met:b'))}
        separate_routes(data,{'red'},pairs)
        report=validate_map_geometry(data,model,nodes,xy)
        self.assertTrue(report['orthogonal_segments_verified'])
        self.assertEqual(data[0]['metadata']['overlap_routing']['parallel_detours']['shifted_segments'],1)
        for nid in original[0]['metadata']['source_node_ids']:
            self.assertEqual(data[1]['nodes'][nid],original[1]['nodes'][nid])
        routed=copy.deepcopy(data)
        separate_routes(data,{'red'},pairs)
        self.assertEqual(data,routed)
        separate_routes(data,enabled=False)
        self.assertTrue(validate_map_geometry(data,model,nodes,xy)['orthogonal_segments_verified'])
        self.assertIsNone(data[0]['metadata']['overlap_routing']['parallel_detours'])
        self.assertEqual(data[0]['map_id'],original[0]['map_id'])
