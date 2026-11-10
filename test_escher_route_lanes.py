"""Lane routing keeps chemistry and endpoints while separating reaction paths."""
import copy
import json
from pathlib import Path
import unittest

from escher_routing import separate_routes
from export_escher import create_map,validate_map
from reaction_network import build_network
from run_fba_escher import read_fba_layout


class RouteLaneTests(unittest.TestCase):
    def fixture(self):
        model=dict(id='lanes',metabolites=[dict(id=x) for x in ('a','b','c')],reactions=[
            dict(id='red',metabolites={'a':1},lower_bound=0,upper_bound=10),
            dict(id='orange',metabolites={'b':1},lower_bound=0,upper_bound=10),
            dict(id='gray',metabolites={'c':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={key:([0,0] if key.startswith('rxn:') else [10,0]) for key in nodes}
        return model,nodes,edges,xy,create_map(model,nodes,edges,xy,scale=18,padding=40)

    def test_priority_offsets_and_original_endpoints(self):
        model,nodes,edges,xy,data=self.fixture();original=copy.deepcopy(data)
        separate_routes(data,{'red'},{frozenset(('rxn:orange','met:b'))})
        check=validate_map(data,model,nodes,xy)
        self.assertTrue(check['orthogonal_segments_verified'])
        self.assertGreater(check['routing_nodes'],0)
        for node in data[0]['metadata']['source_node_ids']:
            self.assertEqual(data[1]['nodes'][node],original[1]['nodes'][node])
        offsets=data[0]['metadata']['overlap_routing']['max_offset_by_priority']
        self.assertEqual(offsets['red'],0)
        self.assertGreater(offsets['orange'],0)
        self.assertGreater(offsets['gray'],offsets['orange'])

    def test_deterministic_and_neutral_reset(self):
        model,nodes,edges,xy,data=self.fixture();other=copy.deepcopy(data)
        separate_routes(data,{'red'},{frozenset(('rxn:orange','met:b'))})
        separate_routes(other,{'red'},{frozenset(('rxn:orange','met:b'))})
        self.assertEqual(data,other)
        first=copy.deepcopy(data)
        separate_routes(data,{'red'},{frozenset(('rxn:orange','met:b'))})
        self.assertEqual(data,first)
        separate_routes(data,enabled=False)
        self.assertTrue(validate_map(data,model,nodes,xy)['orthogonal_segments_verified'])
        self.assertEqual(data[0]['metadata']['overlap_routing']['shifted_legs'],0)
        self.assertFalse(any(data[0]['metadata']['segment_priorities'].values()))

    def test_e_coli_backbones_stay_fixed_and_fba_has_no_priority(self):
        root=Path(__file__).resolve().parent;p=root/'results/e_coli_core'
        if not p.exists():p=root/'public/results/e_coli_core'
        data=json.loads((p/'escher_map.json').read_text());source=json.loads((p/'source_model.json').read_text())
        model,nodes,edges,xy=read_fba_layout(p,source)
        self.assertTrue(validate_map(data,model,nodes,xy)['orthogonal_segments_verified'])
        self.assertEqual(data[0]['metadata']['overlap_routing']['primary_backbone_max_offset'],0)
        report=json.loads((p/'fba_report.json').read_text())
        self.assertFalse(report['pathway_highlighting'])
        self.assertAlmostEqual(report['objective_value'],0.8739215069684304,places=8)


if __name__=='__main__':unittest.main()
