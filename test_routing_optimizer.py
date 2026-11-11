"""Bounded search improves actual geometry while preserving constraints."""
import copy
import unittest
from routing_optimizer import optimize, simple
from routing_metrics import measure, objective
from test_routing_metrics import path


class OptimizerTests(unittest.TestCase):
    def test_crossing_priority_accepts_extra_bends_to_remove_crossing(self):
        paths=[path('horizontal',[(-1,0),(1,0)]),path('vertical',[(0,-1),(0,1)],routing_locked=True)]
        result,xy,report=optimize(paths,{},time_limit=2,route_offset=.3,node_moves=False)
        self.assertEqual(report['before']['crossings'],1)
        self.assertEqual(report['after']['crossings'],0)
        self.assertGreater(report['after']['bends'],report['before']['bends'])
        self.assertTrue(report['nonworsening_verified'])
        self.assertEqual(report['after'],measure(result))
        self.assertFalse(report['global_optimality_proven'])

    def test_weighted_mode_respects_crossing_weight(self):
        paths=[path('horizontal',[(-1,0),(1,0)]),path('vertical',[(0,-1),(0,1)],routing_locked=True)]
        _,_,high=optimize(paths,{},mode='weighted',crossing_weight=10,bend_weight=1,time_limit=1,node_moves=False)
        _,_,low=optimize(paths,{},mode='weighted',crossing_weight=1,bend_weight=10,time_limit=1,node_moves=False)
        self.assertEqual(high['after']['crossings'],0)
        self.assertEqual(low['after']['crossings'],1)

    def test_empty_graph(self):
        paths,xy,report=optimize([],{},time_limit=0)
        self.assertEqual(paths,[])
        self.assertEqual(report['after']['crossings_plus_bends'],0)

    def test_sum_preserves_single_crossing_if_detour_adds_two_bends(self):
        paths=[path('horizontal',[(-1,0),(1,0)]),path('vertical',[(0,-1),(0,1)],routing_locked=True)]
        result,_,report=optimize(paths,{},mode='sum',time_limit=1,node_moves=False)
        self.assertLessEqual(report['after']['crossings_plus_bends'],1)

    def test_zero_budget_and_zero_offset_keep_geometry(self):
        paths=[path('a',[(0,0),(2,0),(2,2)])]
        for budget,offset in ((0,.3),(1,0)):
            result,_,report=optimize(paths,{},time_limit=budget,route_offset=offset)
            self.assertEqual(result,paths)
            self.assertEqual(report['before'],report['after'])

    def test_locked_routes_and_backbone_are_unchanged(self):
        paths=[path('a',[(0,0),(2,0),(2,2)],backbone=True),path('b',[(1,-1),(1,1)],routing_locked=True)]
        result,_,report=optimize(paths,{},time_limit=1)
        self.assertEqual(result,paths)

    def test_no_new_node_hits_or_overlaps_are_introduced(self):
        paths=[path('a',[(0,0),(0,2),(2,2)]),path('b',[(1,-1),(1,1)],routing_locked=True)]
        xy={'block':[.5,0]}
        _,_,report=optimize(paths,xy,time_limit=1,node_moves=False)
        self.assertEqual(report['remaining_node_hits'],[])
        self.assertEqual(report['after']['coincident_intervals'],0)

    def test_no_backtracking_or_diagonal_candidates(self):
        self.assertFalse(simple([(0,0),(2,0),(1,0)]))
        self.assertFalse(simple([(0,0),(1,1)]))
        self.assertTrue(simple([(0,0),(1,0),(1,1)]))

    def test_node_move_improves_without_moving_fixed_nodes(self):
        # A legal node move removes a bend after a route-only detour.
        paths=[path('a',[(0,0),(1,0)],source='a',target='b'),
               path('locked',[(.5,-.5),(.5,.5)],routing_locked=True)]
        xy={'a':[0,0],'b':[1,0]}
        nodes={key:dict(kind='metabolite') for key in xy}
        result,positions,report=optimize(paths,xy,nodes,movable={'a'},cols=3,rows=3,
            mode='crossings-first',time_limit=2,valid_position=lambda key,pos:pos[1]>0)
        self.assertEqual(positions['b'],xy['b'])
        self.assertGreater(report['accepted_node_moves'],0)
        self.assertLessEqual(objective(report['after'],'sum'),objective(report['before'],'sum'))
        self.assertEqual(report['after'],measure(result))

    def test_shared_trunks_stay_shared(self):
        paths=[path('a',[(0,0),(2,0),(2,2)]),path('a',[(0,0),(2,0),(2,3)],target='other')]
        result,_,report=optimize(paths,{},time_limit=1)
        from routing_optimizer import shares
        self.assertTrue(shares(*result))
        self.assertLessEqual(objective(report['after']),objective(report['before']))
