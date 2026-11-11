"""Layout-stage red channels satisfy the same spacing in both solver modes."""
import contextlib
import copy
import io
import math
from pathlib import Path
import json
import unittest
from unittest.mock import patch
from red_route_layout import plan_red_routes,spacing_report,validate_red_spacing
from parallel_routing import _legs,_distance
from grid_layout2 import optimize_grid
from sbml2map import parser


def edge(rid,mid,a,b,**kwargs):
    return dict(rid=rid,mid=mid,source=a,target=b,**kwargs)


class RedRouteLayoutTests(unittest.TestCase):
    def fixture(self):
        xy=dict(a=[1,0],b=[1,4],c=[1,1],d=[1,3])
        return xy,[edge('one','b','a','b',backbone=True),edge('two','d','c','d')]

    def test_default_custom_disabled_and_invalid_cli(self):
        p=parser();self.assertEqual(p.parse_args([]).red_parallel_spacing,1)
        for value in (0,.5,1,2):self.assertEqual(p.parse_args(['--red-parallel-spacing',str(value)]).red_parallel_spacing,value)
        for value in ('-1','nan','inf','-inf'):
            with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):p.parse_args(['--red-parallel-spacing='+value])
        for value in (-1,math.nan,math.inf):
            with self.assertRaises(ValueError):validate_red_spacing(value)

    def test_vertical_horizontal_and_reversed_edges_have_required_spacing(self):
        for rotated in (False,True):
            for reverse in (False,True):
                for spacing in (1,1.5,2):
                    with self.subTest(rotated=rotated,reverse=reverse,spacing=spacing):
                        xy,edges=self.fixture()
                        if rotated:xy={k:list(reversed(p)) for k,p in xy.items()}
                        for e in edges:e['reverse']=reverse
                        before=copy.deepcopy(xy)
                        result=plan_red_routes(xy,edges,spacing)
                        self.assertTrue(result['verified'],result['unresolved'])
                        self.assertEqual(xy,before)
                        self.assertEqual(result['stage'],'grid layout routing-channel allocation')
                        self.assertTrue(spacing_report(result['routes'],spacing)['verified'])
                        for r in result['routes']:
                            self.assertEqual(tuple(r['points'][0]),tuple(xy[r['source']]))
                            self.assertEqual(tuple(r['points'][-1]),tuple(xy[r['target']]))

    def test_short_parallel_and_endpoint_contact_do_not_require_a_lane(self):
        for length in (.99,1):
            xy=dict(a=[0,0],b=[0,3],c=[0,1],d=[0,1+length])
            r=plan_red_routes(xy,[edge('one','b','a','b'),edge('two','d','c','d')])
            self.assertTrue(r['verified'])
            self.assertEqual(any(a['offset'] for a in r['allocations']),length>=1)
        r=plan_red_routes(dict(a=[0,0],b=[0,2],c=[0,4]),[edge('one','b','a','b'),edge('two','c','b','c')])
        self.assertFalse(any(a['offset'] for a in r['allocations']))

    def test_same_reaction_trunks_are_exempt_and_disabled_option_is_fast(self):
        xy,edges=self.fixture();edges[1]['rid']='one'
        r=plan_red_routes(xy,edges)
        self.assertTrue(r['verified']);self.assertFalse(any(a['offset'] for a in r['allocations']))
        r=plan_red_routes(xy,edges,spacing=0)
        self.assertFalse(r['enabled']);self.assertEqual(r['routes'],[])

    def test_protected_nodes_are_avoided_when_channels_are_allocated(self):
        xy=dict(a=[0,0],b=[0,4],blocked=[0,2])
        r=plan_red_routes(xy,[edge('red','b','a','b')],obstacle_node_ids=['blocked'])
        self.assertTrue(r['verified'])
        for leg in _legs(r['routes'][0]['points']):self.assertGreaterEqual(_distance((0,2),leg),.12)

    def test_safe_touching_glucose_branch_stays_straight(self):
        xy=dict(glucose=[1,5],marker=[1,4],product=[1,0],a=[1,1],b=[1,3])
        r=plan_red_routes(xy,[edge('pts','glucose','glucose','marker',backbone=True),
            edge('pts','product','marker','product'),edge('other','b','a','b')])
        self.assertTrue(r['verified'])
        glucose=next(p for p in r['routes'] if p['mid']=='glucose')
        self.assertEqual(glucose['points'],[(1,5),(1,4)])

    def test_time_budget_reports_unresolved_without_claiming_success(self):
        xy,edges=self.fixture()
        with patch('red_route_layout.time.monotonic',side_effect=[0,1,1,1,1,1]):
            r=plan_red_routes(xy,edges,time_limit=.05)
        self.assertFalse(r['verified']);self.assertTrue(r['unresolved'])

    def test_layout_wrapper_applies_to_milp_and_heuristic(self):
        ids=['a','b','c','d'];pairs=[('a','b'),('c','d')]
        fixed={'a':[1,-4],'b':[1,-1],'c':[1,-3],'d':[1,-2]}
        red=[edge('one','b','a','b'),edge('two','d','c','d')]
        for budget in (None,1):
            with self.subTest(max_milp_nodes=budget):
                xy,report=optimize_grid(ids,pairs,4,5,time_limit=2,fixed_offsets=fixed,
                    max_milp_nodes=budget,red_edges=red,red_obstacle_node_ids=ids,
                    progress=lambda *a,**k:None)
                plan=report['red_parallel_layout']
                self.assertTrue(plan['verified'],plan['unresolved'])
                self.assertTrue(spacing_report(plan['routes'],1)['verified'])
                self.assertTrue(report['position_constraints']['fixed_backbone_positions_verified'])

    def test_saved_layout_and_both_viewers_verify_one_grid_spacing(self):
        folder=Path(__file__).resolve().parent/'results/e_coli_core'
        svg=json.loads((folder/'parallel_routing.json').read_text())
        meta=json.loads((folder/'escher_map.json').read_text())[0]['metadata']
        self.assertTrue(svg['red_parallel_spacing']['verified'])
        self.assertGreaterEqual(svg['red_parallel_spacing']['minimum_parallel_gap_grid'],1-1e-7)
        check=meta['overlap_routing']['red_parallel_spacing']
        self.assertTrue(check['verified']);self.assertGreaterEqual(check['minimum_parallel_gap_grid'],1-1e-7)
