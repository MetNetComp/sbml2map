"""Only structurally red/orange nodes constrain display routing."""
import copy
import json
from pathlib import Path
import unittest
from parallel_routing import structural_node_priorities
from export_escher import create_map, validate_map_geometry
from escher_routing import separate_routes
from reaction_network import build_network


class StructuralNodePolicyTests(unittest.TestCase):
    def test_reaction_metabolite_and_cofactor_colors_match_display_policy(self):
        nodes={
            'rxn:red':dict(kind='reaction',original_id='red'),
            'rxn:orange':dict(kind='reaction',original_id='orange'),
            'rxn:gray':dict(kind='reaction',original_id='gray'),
            'met:backbone':dict(kind='metabolite',original_id='backbone'),
            'met:connector':dict(kind='metabolite',original_id='connector'),
            'met:side':dict(kind='metabolite',original_id='side'),
            'copy:red:atp':dict(kind='currency',original_id='atp'),
            'copy:orange:nad':dict(kind='currency',original_id='nad')}
        pairs={frozenset(('rxn:orange','met:connector')),
               frozenset(('rxn:orange','copy:orange:nad')),
               frozenset(('rxn:red','met:backbone'))}
        expected={'rxn:red':2,'rxn:orange':1,'rxn:gray':0,'met:backbone':2,
                  'met:connector':1,'met:side':0,'copy:red:atp':0,'copy:orange:nad':1}
        actual=structural_node_priorities(nodes,{'red'},{'backbone','atp'},pairs)
        self.assertEqual(actual,expected)
        # UI selections and visible/hidden status do not affect structure.
        for node in nodes.values():node.update(selected=True,color='purple',hidden=True)
        self.assertEqual(structural_node_priorities(nodes,{'red'},{'backbone','atp'},pairs),actual)

    def fixture(self,obstacle_kind='reaction'):
        model=dict(id='obstacles',metabolites=[dict(id=x) for x in ('a','b')],reactions=[
            dict(id='red',metabolites={'a':1},lower_bound=0,upper_bound=10),
            dict(id='other',metabolites={'b':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'rxn:red':[0,0],'met:a':[4,0],'rxn:other':[2,0],'met:b':[2,1]}
        if obstacle_kind=='metabolite':xy['rxn:other'],xy['met:b']=xy['met:b'],xy['rxn:other']
        return model,nodes,xy,create_map(model,nodes,edges,xy,scale=180,padding=40)

    def test_gray_nodes_can_be_crossed_but_orange_nodes_are_avoided(self):
        for kind in ('reaction','metabolite'):
            for orange in (False,True):
                with self.subTest(kind=kind,orange=orange):
                    model,nodes,xy,data=self.fixture(kind)
                    pairs={frozenset(('rxn:other','met:b'))} if orange else set()
                    separate_routes(data,{'red'},pairs)
                    meta=data[0]['metadata'];report=meta['overlap_routing']
                    self.assertEqual(report['node_obstacle_policy'],'red_orange_structural_nodes')
                    self.assertEqual(report['node_detours']['shifted_segments'],1 if orange else 0)
                    self.assertEqual('rxn:other' in report['obstacle_node_ids'],orange)
                    self.assertEqual('met:b' in report['obstacle_node_ids'],orange)
                    self.assertTrue(validate_map_geometry(data,model,nodes,xy)['orthogonal_segments_verified'])

    def test_red_metabolites_are_avoided_even_on_otherwise_gray_reactions(self):
        model,nodes,xy,data=self.fixture('metabolite')
        data[0]['metadata']['sbml2map']=dict(major_backbones=dict(routes=[dict(metabolite_ids=['a','b'],steps=[])]))
        separate_routes(data,{'red'})
        report=data[0]['metadata']['overlap_routing']
        self.assertIn('met:b',report['obstacle_node_ids'])
        self.assertNotIn('rxn:other',report['obstacle_node_ids'])
        self.assertEqual(report['node_detours']['shifted_segments'],1)
        self.assertTrue(validate_map_geometry(data,model,nodes,xy)['orthogonal_segments_verified'])

    def test_saved_svg_and_escher_use_identical_protected_source_nodes(self):
        folder=Path(__file__).resolve().parent/'results/e_coli_core'
        svg=json.loads((folder/'parallel_routing.json').read_text())
        escher=json.loads((folder/'escher_map.json').read_text())[0]['metadata']
        self.assertEqual(svg['obstacle_node_ids'],escher['overlap_routing']['obstacle_node_ids'])
        self.assertEqual(svg['structural_node_priorities'],escher['structural_node_priorities'])
        self.assertIn('met:succ_c',svg['obstacle_node_ids'])
        self.assertTrue(any(rank==0 for rank in svg['structural_node_priorities'].values()))
