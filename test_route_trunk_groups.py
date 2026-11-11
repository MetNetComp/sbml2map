"""Shared endpoints alone must not bundle unrelated straight intervals."""
import copy
import csv
import json
from pathlib import Path
import re
import unittest
import xml.etree.ElementTree as ET
from parallel_routing import separate_parallel,avoid_nodes


class TrunkGroupTests(unittest.TestCase):
    def test_node_avoidance_leaves_touching_safe_branch_straight(self):
        for rotated in (False,True):
            for reverse in (False,True):
                with self.subTest(rotated=rotated,reverse=reverse):
                    transform=lambda p:(p[1],p[0]) if rotated else p
                    safe=dict(rid='same',priority=2,points=[(0,0),(1,0)])
                    blocked=dict(rid='same',priority=2,points=[(1,0),(4,0)])
                    for p in [safe,blocked]:
                        p['points']=[transform(x) for x in p['points']]
                        if reverse:p['points'].reverse()
                    before=copy.deepcopy(safe)
                    result=avoid_nodes([safe,blocked],obstacles=[transform((2,0))])
                    self.assertEqual(safe,before)
                    self.assertEqual(result['shifted_segments'],1)
                    self.assertEqual(len(blocked['points']),4)

    def test_parallel_separation_leaves_touching_safe_branch_straight(self):
        for rotated in (False,True):
            with self.subTest(rotated=rotated):
                transform=lambda p:(p[1],p[0]) if rotated else p
                red=dict(rid='red',priority=2,backbone=True,points=[(1.5,.1),(3.5,.1)])
                safe=dict(rid='same',priority=1,points=[(0,0),(1,0)])
                blocked=dict(rid='same',priority=1,points=[(1,0),(4,0)])
                for p in [red,safe,blocked]:p['points']=[transform(x) for x in p['points']]
                before=copy.deepcopy(safe)
                result=separate_parallel([red,safe,blocked])
                self.assertEqual(safe,before)
                self.assertEqual(result['shifted_segments'],1)
                self.assertGreater(len(blocked['points']),2)

    def test_actual_overlapping_trunks_stay_bundled(self):
        a=dict(rid='same',priority=2,points=[(0,0),(3,0)])
        b=dict(rid='same',priority=2,points=[(2,0),(4,0)])
        report=avoid_nodes([a,b],obstacles=[(1,0)])
        self.assertEqual(report['shifted_segments'],2)
        self.assertEqual(a['points'][1][1],b['points'][1][1])
        self.assertEqual(report['trunk_grouping'],'positive_length_overlap')


class GlucoseBranchRegressionTests(unittest.TestCase):
    def test_saved_glucose_branch_is_straight_in_svg_and_escher(self):
        folder=Path(__file__).resolve().parent/'results/e_coli_core'
        with (folder/'edges.csv').open() as f:edges=list(csv.DictReader(f))
        index=next(i for i,e in enumerate(edges) if e['reaction_id']=='GLCpts' and e['metabolite_id']=='glc__D_e')
        svg=ET.parse(folder/'layout.svg')
        d=svg.find(f'.//*[@id="edge-{index}"]').find('.//{*}path').get('d')
        numbers=list(map(float,re.findall(r'[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?',d)))
        points=list(zip(numbers[::2],numbers[1::2]))
        self.assertEqual(len(points),2)
        self.assertAlmostEqual(points[0][0],points[1][0])
        self.assertNotEqual(points[0][1],points[1][1])
        data=json.loads((folder/'escher_map.json').read_text());meta=data[0]['metadata'];nodes=data[1]['nodes']
        reaction=next(r for r in data[1]['reactions'].values() if r['bigg_id']=='GLCpts')
        segments=[s for sid,s in reaction['segments'].items() if 'met:glc__D_e' in meta['segment_source_edges'][sid]]
        self.assertEqual(len(segments),1)
        a,b=[nodes[segments[0][key]] for key in ('from_node_id','to_node_id')]
        self.assertAlmostEqual(a['x'],b['x'])
        self.assertNotEqual(a['y'],b['y'])
