"""Canonical routes survive every export without renderer-specific rerouting."""
import copy
import json
import re
import tempfile
from pathlib import Path
import unittest
from export_escher import create_map, validate_map, write_pathway_preview
from export_sbml_layout import _incidences
from escher_routing import separate_routes
from reaction_network import build_network
from routing_geometry import transformed_paths, orient_for_edges
from routing_metrics import measure
from sbml2map import parser


class SharedGeometryTests(unittest.TestCase):
    def fixture(self):
        model=dict(id='shared',metabolites=[dict(id=m) for m in ('a','b')],reactions=[
            dict(id='r',metabolites={'a':-1,'b':1},lower_bound=-10,upper_bound=0)])
        nodes,edges=build_network(model,set())
        xy={'met:a':[0,0],'met:b':[3,3],'rxn:r':[1,1]}
        paths=[dict(rid='r',mid='a',source='met:a',target='rxn:r',priority=0,
                    points=[(0,0),(-1,0),(-1,1),(1,1)]),
               dict(rid='r',mid='b',source='rxn:r',target='met:b',priority=0,
                    points=[(1,1),(1,4),(3,4),(3,3)])]
        return model,nodes,edges,xy,paths

    def test_escher_and_sbml_incidence_reader_preserve_canonical_geometry(self):
        model,nodes,edges,xy,paths=self.fixture()
        data=create_map(model,nodes,edges,xy)
        separate_routes(data,route_geometry=paths,routing_report={'optimization': {'after':measure(paths)}})
        self.assertTrue(validate_map(data,model,nodes,xy)['orthogonal_segments_verified'])
        expected={(p['rid'],p['mid']):p for p in transformed_paths(paths,data[0]['metadata']['coordinate_transform'])}
        for rid,target,keys in _incidences(data):
            mid=data[1]['nodes'][target]['bigg_id'];p=expected[(rid,mid)]
            points=list(p['points'])
            if p['source']!='rxn:'+rid:points.reverse()
            actual=[tuple(data[1]['nodes'][k][c] for c in ('x','y')) for k in keys]
            self.assertEqual(actual,points)
        ordered=orient_for_edges(list(reversed(paths)),edges,xy)
        self.assertEqual(len(ordered),len(edges))

    def test_pathway_preview_does_not_replace_optimized_routes(self):
        model,nodes,edges,xy,paths=self.fixture();data=create_map(model,nodes,edges,xy)
        separate_routes(data,route_geometry=paths,routing_report={'routing_preference':'crossings-first'})
        before=copy.deepcopy(data[1])
        with tempfile.TemporaryDirectory() as tmp:
            file=Path(tmp)/'preview.html';write_pathway_preview(file,data)
            embedded=json.loads(re.search(r'const mapData=(.*);\n',file.read_text()).group(1))
        self.assertEqual(embedded[1],before)

    def test_objective_options_and_defaults(self):
        p=parser();args=p.parse_args([])
        self.assertEqual(args.routing_objective,'crossings-first')
        args=p.parse_args(['--routing-objective','weighted','--crossing-weight','100','--bend-weight','2','--routing-time-limit','5','--no-routing-node-moves'])
        self.assertEqual((args.crossing_weight,args.bend_weight,args.routing_time_limit),(100,2,5))
        self.assertTrue(args.no_routing_node_moves)
