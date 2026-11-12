"""Bounded compaction keeps red fixed and rejects an orange-length increase."""
import copy
import unittest
from unittest.mock import patch
from reaction_network import build_network
from exchange_corridors import place_corridors
from exchange_alignment import validate_alignment
from boundary_compaction import compact_boundary
from escher_grid import attach_coenzymes
from sbml2map import parser


class BoundaryCompactionTests(unittest.TestCase):
    def fixture(self):
        model=dict(id='small',metabolites=[dict(id=k,name=k,compartment='e' if k.endswith('_e') else 'c') for k in ['a_e','a_c','b_c']],
                   reactions=[dict(id='EX_a',name='EX_a',metabolites={'a_e':-1},lower_bound=-10,upper_bound=10),
                              dict(id='T',name='T',metabolites={'a_e':-1,'a_c':1},lower_bound=-10,upper_bound=10),
                              dict(id='R',name='R',metabolites={'b_c':-1,'a_c':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'rxn:EX_a':[10.,4.],'met:a_e':[9.25,4.],'rxn:T':[8.5,4.],
            'rxn:R':[3.,4.],'met:b_c':[3.,5.],'met:a_c':[5.,4.]}
        red=[dict(rid='R',mid=mid,source='rxn:R',target='met:'+mid,priority=2,backbone=True,points=[xy['rxn:R'],xy['met:'+mid]]) for mid in ['a_c','b_c']]
        b=dict(fixed_node_offsets={k:xy[k] for k in ['rxn:R','met:b_c','met:a_c']},
               routes=[dict(id='red',reaction_ids=['R'],metabolite_ids=['a_c','b_c'],steps=[])],
               red_first=dict(routes=red),tca_exclusion=dict(region=None),routing_boundary=[0,10,0,10],
               exchange_connectors=dict(edge_indices=[i for i,e in enumerate(edges) if e['reaction_id']!='R']))
        groups=[dict(reaction_id='EX_a',exchange='rxn:EX_a',external='met:a_e',internal='met:a_c')]
        b['exchange_corridors']=place_corridors(model,nodes,xy,b,{'rxn:EX_a'},groups,11,11)
        return model,nodes,edges,xy,b,groups

    def run_compaction(self,fixture,**kwargs):
        m,n,e,xy,b,g=fixture
        return compact_boundary(m,n,e,xy,b,{'rxn:EX_a'},g,11,11,{'R'},attach_coenzymes,attempts=0,progress=lambda *a,**k:None,**kwargs)

    def test_lower_and_right_compaction_keeps_red_and_reduces_length(self):
        fixture=self.fixture();xy,b=fixture[3],fixture[4];before=copy.deepcopy(xy);red=copy.deepcopy(b['red_first'])
        report=self.run_compaction(fixture)
        self.assertTrue(report['accepted']);self.assertEqual(report['after_boundary'],[0,6,3,10])
        self.assertLess(report['after_orange_length_grid'],report['before_orange_length_grid'])
        self.assertEqual(b['red_first'],red)
        for k in ['rxn:R','met:a_c','met:b_c']:self.assertEqual(xy[k],before[k])
        self.assertTrue(validate_alignment(fixture[5],xy,11,11,b['routing_boundary'])['verified'])
        self.assertTrue(b['crossing_routing']['report']['boundary_validation']['verified'])

    def test_length_increase_is_rejected_without_moving_nodes(self):
        fixture=self.fixture();before=copy.deepcopy(fixture[3])
        with patch('boundary_compaction.orange_length',side_effect=[5.]+[100.]*3):
            report=self.run_compaction(fixture,max_candidates=3)
        self.assertFalse(report['accepted']);self.assertEqual(fixture[3],before)
        self.assertTrue(all('increase' in c['reason'] for c in report['rejected_candidates']))

    def test_red_movement_is_rejected(self):
        fixture=self.fixture();before=copy.deepcopy(fixture[3])
        def moving_attach(nodes,positions,**kwargs):positions['met:a_c'][0]+=.25
        m,n,e,xy,b,g=fixture
        report=compact_boundary(m,n,e,xy,b,{'rxn:EX_a'},g,11,11,{'R'},moving_attach,attempts=0,max_candidates=1,progress=lambda *a,**k:None)
        self.assertFalse(report['accepted']);self.assertEqual(xy,before)
        self.assertIn('Frozen red node moved',report['rejected_candidates'][0]['reason'])

    def test_cofactor_candidates_respect_shifted_boundary(self):
        nodes={'rxn:R':dict(kind='reaction'),'copy:a':dict(kind='currency',parent_id='rxn:R'),'copy:b':dict(kind='currency',parent_id='rxn:R')}
        xy={'rxn:R':[5.9,3.2]}
        attach_coenzymes(nodes,xy,grid_shape=(11,11),boundary=[0,6,3,10])
        for k in ['copy:a','copy:b']:
            x,y=xy[k];self.assertTrue(0<=x<=6 and 3<=y<=10)
            self.assertAlmostEqual(abs(x-5.9)+abs(y-3.2),.45)

    def test_compaction_can_be_disabled(self):
        self.assertTrue(parser().parse_args([]).compact_boundary)
        self.assertFalse(parser().parse_args(['--no-compact-boundary']).compact_boundary)


    def test_candidate_constraint_failure_is_rejected(self):
        fixture=self.fixture();before=copy.deepcopy(fixture[3])
        with patch('boundary_compaction.validate_exclusion',side_effect=RuntimeError('TCA exclusion failed')):
            report=self.run_compaction(fixture,max_candidates=1)
        self.assertFalse(report['accepted']);self.assertEqual(fixture[3],before)
        self.assertIn('TCA exclusion failed',report['rejected_candidates'][0]['reason'])
