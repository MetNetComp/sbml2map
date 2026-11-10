"""The oxidative TCA order must be clockwise in both output coordinate systems."""
import csv
import json
from pathlib import Path
import unittest
import xml.etree.ElementTree as ET

from major_backbones2 import backbone_plan, pathway_definitions, validate_backbones, clockwise_square
from model_input import load_model
from reaction_network import build_network
from sbml2map import BacterialStrategy

BASE = Path(__file__).resolve().parent
MODEL = BASE/'models/e_coli_core.xml' if (BASE/'models').is_dir() else BASE/'e_coli_core.xml'


def area2(points):
    return sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(points,points[1:]+points[:1]))


class TCAOrientationTests(unittest.TestCase):
    def test_square_starts_at_top_left_and_moves_right(self):
        ring = clockwise_square(4, -20, 5)
        self.assertEqual(ring[:3], [(4,-20),(5,-20),(6,-20)])
        self.assertEqual(len(set(ring)), 20)
        self.assertEqual(area2(ring), -50)

    def test_discovered_and_reference_routes_are_clockwise(self):
        model=load_model(MODEL);nodes,_=build_network(model,[])
        strategy=BacterialStrategy(model)
        routes,_=pathway_definitions(model)
        for fixed,plan in (strategy.plan(strategy.routes,nodes),backbone_plan(routes,nodes)):
            with self.subTest(plan=plan['geometry_scope']):
                route=next(r for r in plan['routes'] if r['id']=='tca')
                self.assertLess(area2(list(route['node_offsets'].values())),0)
                self.assertGreater(fixed['met:cit_c'][0],fixed['met:oaa_c'][0])
                xy={k:[x,y+40] for k,(x,y) in fixed.items()}
                self.assertTrue(validate_backbones(plan,xy,41)['all_displayed_backbone_positions_verified'])

    def test_validation_rejects_reflected_counterclockwise_cycle(self):
        model=load_model(MODEL);nodes,_=build_network(model,[]);strategy=BacterialStrategy(model)
        _,plan=strategy.plan(strategy.routes,nodes)
        route=next(r for r in plan['routes'] if r['id']=='tca')
        square=plan['square'];axis=2*square['left']+square['side_length']
        for key,(x,y) in route['node_offsets'].items():
            route['node_offsets'][key]=[axis-x,y]
            plan['fixed_node_offsets'][key]=[axis-x,y]
        xy={k:[x,y+40] for k,(x,y) in plan['fixed_node_offsets'].items()}
        with self.assertRaisesRegex(RuntimeError,'clockwise'):
            validate_backbones(plan,xy,41)

    def test_saved_coordinates_and_escher_are_clockwise(self):
        folder=BASE/'results/e_coli_core'
        if not folder.is_dir():folder=BASE/'public/results/e_coli_core'
        plan=json.loads((folder/'major_backbones.json').read_text())
        route=next(r for r in plan['routes'] if r['id']=='tca')
        keys=route['displayed_node_ids']
        with (folder/'coordinates.csv').open() as stream:
            rows=list(csv.DictReader(stream))
        xy={r['id']:[float(r['x']),float(r['y'])] for r in rows}
        self.assertLess(area2([xy[k] for k in keys]),0)
        data=json.loads((folder/'escher_map.json').read_text())
        exported={source:data[1]['nodes'][node] for node,source in data[0]['metadata']['source_node_ids'].items()}
        self.assertGreater(area2([[exported[k]['x'],exported[k]['y']] for k in keys]),0)
        self.assertEqual(plan['square']['orientation'],'clockwise')
        indices={r['id']:i for i,r in enumerate(rows)}
        for filename in ('layout.svg','index.html'):
            with self.subTest(output=filename):
                text=(folder/filename).read_text()
                tree=ET.fromstring(text[text.index('<svg'):text.index('</svg>')+6])
                groups={e.get('id'):e for e in tree.iter() if e.get('id')}
                points=[]
                for key in keys:
                    group=groups['node-'+str(indices[key])]
                    mark=next(e for e in group.iter() if e.tag.endswith('}use') and 'x' in e.attrib)
                    points.append([float(mark.get('x')),float(mark.get('y'))])
                self.assertGreater(area2(points),0)


if __name__=='__main__':
    unittest.main()
