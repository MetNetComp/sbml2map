"""Important routes prefer few bends; gray geometry yields instead."""
import copy
import csv
import json
from pathlib import Path
import re
import unittest
import xml.etree.ElementTree as ET
from parallel_routing import avoid_nodes,_legs,_distance
from routing_simplification import simplify_important,yield_gray,gray_conflicts


def fixture(rotated=False,reverse=False):
    important=dict(rid='EX_ac_e',priority=1,points=[(0,0),(0,.3),(1.85,.3),(1.85,0),(2,0)])
    gray=dict(rid='PGL',priority=0,points=[(2,-1),(2,1)])
    obstacles=[(0,0),(1,0),(2,0)]
    if rotated:
        for p in [important,gray]:p['points']=[(y,x) for x,y in p['points']]
        obstacles=[(y,x) for x,y in obstacles]
    if reverse:important['points'].reverse()
    return important,gray,obstacles


class ImportantSimplificationTests(unittest.TestCase):
    def test_acetate_uses_two_bends_and_gray_yields_in_both_orientations(self):
        for rotated in (False,True):
            for reverse in (False,True):
                with self.subTest(rotated=rotated,reverse=reverse):
                    important,gray,obstacles=fixture(rotated,reverse)
                    endpoints=(important['points'][0],important['points'][-1])
                    report=simplify_important([important,gray],obstacles)
                    self.assertEqual(len(important['points']),4)
                    self.assertEqual(report['removed_bends'],1)
                    self.assertEqual((important['points'][0],important['points'][-1]),endpoints)
                    before=copy.deepcopy(important)
                    report=yield_gray([important,gray],obstacles)
                    self.assertEqual(important,before)
                    self.assertEqual(report['changed_paths'],1)
                    self.assertEqual(gray_conflicts(gray,gray['points'],[important]),0)
                    for p in obstacles:
                        if p not in (gray['points'][0],gray['points'][-1]):
                            for leg in _legs(gray['points']):self.assertGreaterEqual(_distance(p,leg),.12-1e-7)

    def test_gray_lines_do_not_force_endpoint_step_fallbacks(self):
        important,gray,obstacles=fixture()
        important['points']=[(0,0),(2,0)]
        report=avoid_nodes([important,gray],obstacles=obstacles)
        self.assertEqual(len(important['points']),4)
        self.assertEqual(report['endpoint_step_fallback_segments'],0)

    def test_shortest_route_wins_among_equal_bend_counts(self):
        p=dict(rid='orange',priority=1,points=[(0,0),(0,2),(3,2),(3,0)])
        simplify_important([p],obstacles=[(1.5,0)])
        self.assertEqual(len(p['points']),4)
        self.assertAlmostEqual(sum(abs(b[0]-a[0])+abs(b[1]-a[1]) for a,b in zip(p['points'],p['points'][1:])),3.6)

    def test_preallocated_red_routes_are_unchanged(self):
        p=dict(rid='red',priority=2,routing_locked=True,points=[(0,0),(-1,0),(-1,3),(0,3)])
        gray=dict(rid='gray',priority=0,points=[(-1,0),(-1,3)])
        before=copy.deepcopy(p)
        simplify_important([p,gray]);yield_gray([p,gray])
        self.assertEqual(p,before)

    def test_positive_length_shared_trunks_remain_bundled(self):
        a,_,obstacles=fixture();b=copy.deepcopy(a);b['points'][-1]=(3,0)
        before=copy.deepcopy([a,b]);simplify_important([a,b],obstacles)
        self.assertEqual([a,b],before)

    def test_gray_yield_failure_does_not_add_important_bends(self):
        p,gray,obstacles=fixture();simplify_important([p,gray],obstacles)
        gray2=copy.deepcopy(gray);gray2['points'][-1]=(2,2)
        before=copy.deepcopy(p)
        report=yield_gray([p,gray,gray2],obstacles)
        self.assertEqual(p,before)
        self.assertTrue(report['unresolved'])

    def test_zero_offset_preserves_existing_paths(self):
        p,gray,obstacles=fixture();before=copy.deepcopy([p,gray])
        simplify_important([p,gray],obstacles,route_offset=0);yield_gray([p,gray],obstacles,route_offset=0)
        self.assertEqual([p,gray],before)


class AcetateArtifactRegressionTests(unittest.TestCase):
    def test_saved_svg_and_escher_acetate_have_exactly_two_bends(self):
        folder=Path(__file__).resolve().parent/'results/e_coli_core'
        data=json.loads((folder/'escher_map.json').read_text());reaction=next(r for r in data[1]['reactions'].values() if r['bigg_id']=='EX_ac_e')
        self.assertEqual(len(reaction['segments']),3)
        with (folder/'edges.csv').open() as f:edges=list(csv.DictReader(f))
        index=next(i for i,e in enumerate(edges) if e['reaction_id']=='EX_ac_e')
        svg=ET.parse(folder/'layout.svg')
        d=svg.find(f'.//*[@id="edge-{index}"]').find('.//{*}path').get('d')
        numbers=list(map(float,re.findall(r'[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?',d)))
        points=list(zip(numbers[::2],numbers[1::2]))
        self.assertEqual(len(points),4)
        legs=list(_legs(points));self.assertEqual(len(legs),3)
        self.assertEqual(legs[0]['axis'],legs[2]['axis'])
        self.assertNotEqual(legs[0]['axis'],legs[1]['axis'])
        self.assertAlmostEqual(points[0][1],points[-1][1])
        meta=data[0]['metadata'];pri=meta['segment_priorities']
        important=[]
        for r in data[1]['reactions'].values():
            if any(pri[sid]>0 for sid in r['segments']):
                for sid,s in r['segments'].items():
                    if pri[sid]>0:important.append(dict(rid=r['bigg_id'],priority=pri[sid],points=[tuple(data[1]['nodes'][s[k]][c] for c in ('x','y')) for k in ('from_node_id','to_node_id')]))
        pgl=next(r for r in data[1]['reactions'].values() if r['bigg_id']=='PGL')
        for s in pgl['segments'].values():
            gray=dict(rid='PGL',priority=0,points=[tuple(data[1]['nodes'][s[k]][c] for c in ('x','y')) for k in ('from_node_id','to_node_id')])
            # At least the acetate overlap must be gone; unrelated gray crossings are allowed.
            self.assertEqual(gray_conflicts(gray,gray['points'],[r for r in important if r['rid']=='EX_ac_e']),0)
