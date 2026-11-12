"""Overlap relocation preserves biology and constraints; point crossings stay free."""
import copy
import unittest
from unittest.mock import patch
from reaction_network import build_network
from crossing_routing import plan_routes
from exchange_corridors import validate_corridors
from orange_reposition import orange_overlaps, relocation_candidates, reposition_orange, place_gray_inside_tca
from sbml2map import parser


class OrangeRepositionTests(unittest.TestCase):
    def path(self, rid, points, priority=1, mid=None):
        return dict(rid=rid, mid=mid or rid, source='rxn:'+rid, target='met:'+(mid or rid), priority=priority, points=points)

    def fixture(self):
        model=dict(id='small',metabolites=[dict(id=k,name=k,compartment='e' if k.endswith('_e') else 'c') for k in ['x_e','sugar_c','next_c']],
                   reactions=[dict(id='EX_x',name='EX_x',metabolites={'x_e':-1},lower_bound=-10,upper_bound=10),
                              dict(id='IMPORT',name='IMPORT',metabolites={'x_e':-1,'sugar_c':1},lower_bound=0,upper_bound=10),
                              dict(id='RED',name='RED',metabolites={'sugar_c':-1,'next_c':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'rxn:EX_x':[2.,10.],'met:x_e':[2.,4.],'rxn:IMPORT':[3.,5.],
            'met:sugar_c':[1.,5.],'rxn:RED':[1.,6.],'met:next_c':[1.,7.]}
        red=[dict(rid='RED',mid=mid,source='rxn:RED',target='met:'+mid,priority=2,points=[xy['rxn:RED'],xy['met:'+mid]]) for mid in ['sugar_c','next_c']]
        group=dict(kind='pair',reaction_id='EX_x',exchange='rxn:EX_x',external='met:x_e',primary_transport=None,
                   node_ids=['rxn:EX_x','met:x_e'],corridor=[xy['rxn:EX_x'],xy['met:x_e']])
        b=dict(fixed_node_offsets={p['target']:xy[p['target']] for p in red},routes=[],red_first=dict(routes=red),
               routing_boundary=[0,10,0,10],exchange_corridors=dict(groups=[group]),
               exchange_connectors=dict(edge_indices=[i for i,e in enumerate(edges) if e['reaction_id']!='RED']))
        b['crossing_routing']=plan_routes(nodes,edges,xy,b,{'RED'},attempts=0,progress=lambda *a,**k:None)
        for p in b['crossing_routing']['routes']:
            if p['rid']=='IMPORT' and p['mid']=='x_e':
                p['points']=[[3.,5.],[2.,5.],[2.,4.]]
        return nodes,edges,xy,b

    def test_union_overlap_includes_same_reaction_but_not_point_crossings(self):
        paths=[self.path('A',[[0,0],[4,0]]),self.path('A',[[1,0],[3,0]],mid='other'),
               self.path('B',[[2,0],[5,0]]),self.path('C',[[3,-1],[3,1]]),
               self.path('RED',[[0,0],[5,0]],2),self.path('GRAY',[[0,0],[5,0]],0)]
        overlap=orange_overlaps(paths)
        self.assertEqual(overlap['length_grid'],4.)
        self.assertEqual(len(overlap['segments']),2)
        self.assertTrue(any(s['same_reaction'] for s in overlap['segments']))

    def test_generic_external_relocation_removes_overlap_and_keeps_red_and_pair_straight(self):
        nodes,edges,xy,b=self.fixture();before=copy.deepcopy(xy);red=copy.deepcopy(b['red_first'])
        self.assertGreater(orange_overlaps(b['crossing_routing']['routes'])['length_grid'],0)
        result=reposition_orange(nodes,edges,xy,b,{'RED'},[],11,11,attempts=0,max_candidates=6,progress=lambda *a,**k:None)
        self.assertTrue(result['accepted_moves'])
        self.assertEqual(result['after_overlap_grid'],0.)
        self.assertEqual(xy['met:x_e'][0],xy['rxn:EX_x'][0])
        self.assertEqual(xy['rxn:EX_x'],before['rxn:EX_x'])
        self.assertEqual(b['red_first'],red)
        self.assertTrue(b['crossing_routing']['report']['red_routes_preserved'])
        self.assertTrue(validate_corridors(b['crossing_routing']['routes'],b['exchange_corridors'])['verified'])
        self.assertGreater(result['attempted_candidates'],0)
        self.assertEqual(len(edges),5)

    def test_no_overlap_leaves_layout_unchanged(self):
        nodes,edges,xy,b=self.fixture();before=copy.deepcopy(xy)
        with patch('orange_reposition.orange_overlaps',return_value=dict(length_grid=0.,segments=[])):
            result=reposition_orange(nodes,edges,xy,b,{'RED'},[],11,11,attempts=0,progress=lambda *a,**k:None)
        self.assertEqual(xy,before);self.assertEqual(result['attempted_candidates'],0)

    def test_infeasible_candidates_are_discarded_without_mutating_layout(self):
        nodes,edges,xy,b=self.fixture();before=copy.deepcopy(xy);spec=copy.deepcopy(b['exchange_corridors'])
        with patch('orange_reposition.plan_routes',side_effect=ValueError('red crossing')):
            result=reposition_orange(nodes,edges,xy,b,{'RED'},[],11,11,attempts=0,max_candidates=3,progress=lambda *a,**k:None)
        self.assertEqual(xy,before);self.assertEqual(b['exchange_corridors'],spec)
        self.assertGreaterEqual(result['attempted_candidates'],3)
        self.assertEqual(len(result['rejected_candidates']),result['attempted_candidates'])

    def test_triplet_and_frozen_nodes_are_never_relocation_candidates(self):
        nodes,edges,xy,b=self.fixture()
        b['exchange_corridors']['groups'][0]['kind']='triplet'
        candidates=relocation_candidates(nodes,xy,b,b['crossing_routing']['routes'])
        fixed={k for p in b['red_first']['routes'] for k in (p['source'],p['target'])}
        self.assertTrue(all(k not in fixed|{'met:x_e','rxn:EX_x'} for k,p in candidates))

    def test_runs_when_boundary_compaction_is_disabled(self):
        nodes,edges,xy,b=self.fixture();b.pop('routing_boundary')
        result=reposition_orange(nodes,edges,xy,b,{'RED'},[],11,11,attempts=0,max_candidates=6,progress=lambda *a,**k:None)
        self.assertEqual(b['routing_boundary'],[0,10,0,10])
        self.assertEqual(result['after_overlap_grid'],0.)

    def test_same_reaction_junction_moves_left_to_remove_shared_trunk(self):
        nodes,edges,xy,b=self.fixture()
        xy['met:x_e']=[2.,6.]
        b['exchange_corridors']['groups'][0]['corridor']=[xy['rxn:EX_x'],xy['met:x_e']]
        b['crossing_routing']=plan_routes(nodes,edges,xy,b,{'RED'},attempts=1,progress=lambda *a,**k:None)
        for p in b['crossing_routing']['routes']:
            if p['rid']=='IMPORT':
                p['points']=[[3.,5.],[2.,5.],[2.,6.]] if p['mid']=='x_e' else [[3.,5.],[1.,5.]]
        self.assertTrue(any(s['same_reaction'] for s in orange_overlaps(b['crossing_routing']['routes'])['segments']))
        result=reposition_orange(nodes,edges,xy,b,{'RED'},[],11,11,attempts=1,max_candidates=4,max_passes=1,progress=lambda *a,**k:None)
        self.assertEqual(xy['rxn:IMPORT'],[2.,5.])
        self.assertEqual(result['after_overlap_grid'],0.)

    def test_search_budget_is_per_overlap_and_expands_when_needed(self):
        nodes,edges,xy,b=self.fixture()
        with patch('orange_reposition.plan_routes',side_effect=ValueError('red crossing')):
            result=reposition_orange(nodes,edges,xy,b,{'RED'},[],11,11,attempts=0,max_candidates=1,progress=lambda *a,**k:None)
        self.assertGreater(result['attempted_candidates'],1)
        self.assertEqual(result['overlap_searches'][0]['radii_tried'],[1,3,6])
        self.assertIn('bounded search',result['overlap_searches'][0]['reason'])

    def test_gray_node_can_move_inside_without_moving_red_anchors(self):
        model=dict(metabolites=[dict(id=k,name=k,compartment='c') for k in ['a','b','c']],
                   reactions=[dict(id='RED',metabolites={'a':-1,'b':1,'c':1},lower_bound=0,upper_bound=10),
                              dict(id='GRAY',metabolites={'a':-1,'b':1,'c':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'rxn:RED':[3.,1.],'met:a':[5.,3.],'met:b':[3.,5.],'met:c':[1.,3.],'rxn:GRAY':[7.,3.]}
        before=copy.deepcopy(xy)
        region=dict(shape='circle',center=[3.,3.],radius=2.,node_ids=['rxn:RED','met:a','met:b','met:c'],interior_policy='allow')
        b=dict(routes=[dict(id='red',reaction_ids=['RED'],steps=[])],routing_boundary=[0,8,0,8],tca_exclusion=dict(region=region))
        result=place_gray_inside_tca(nodes,edges,xy,b)
        self.assertEqual(xy['rxn:GRAY'],[3.,3.])
        self.assertTrue(result['moves'])
        for k in region['node_ids']:self.assertEqual(xy[k],before[k])

    def test_can_disable_automatic_relocation(self):
        self.assertTrue(parser().parse_args([]).reposition_orange)
        self.assertFalse(parser().parse_args(['--no-reposition-orange']).reposition_orange)
