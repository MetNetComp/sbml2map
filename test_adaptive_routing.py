"""New routing rules are checked against actual geometry, not stored counts."""
import copy
import math
import unittest
from adaptive_routing import Grid,route_paths,validate_paths,segment_clearance,integer,incident_exempt,parallel_conflicts,validate_node_grid,attach_grid_cofactors
from adaptive_geometry import compact_tca
from tca_geometry import configure_plan
from sbml2map import parser,BacterialStrategy
from model_input import load_model
from pathlib import Path


def path(rid,source,target,a,b,priority):
    return dict(rid=rid,mid=target,source=source,target=target,points=[a,b],priority=priority,backbone=False)


class AdaptiveRoutingTests(unittest.TestCase):
    def test_long_uncrossed_route_allows_audited_red_crossing(self):
        paths=[path('RED','r1','r2',(1,2),(3,2),2),path('ORANGE','o1','o2',(2,1),(2,3),1)]
        xy={k:p for route in paths for k,p in [(route['source'],route['points'][0]),(route['target'],route['points'][-1])]}
        before=copy.deepcopy(xy);exceptions=route_paths(paths,xy,[-2,6,-2,6],None)
        self.assertEqual(len(exceptions),1)
        self.assertEqual(exceptions[0]['crossed_length_grid'],2.)
        self.assertGreaterEqual(exceptions[0]['uncrossed_length_grid'],6.)
        report=validate_paths(paths,xy,[-2,6,-2,6],None,exceptions)
        self.assertEqual(report['metrics']['red_orange'],1)
        self.assertEqual(xy,before)
        with self.assertRaisesRegex(ValueError,'Unaudited'):validate_paths(paths,xy,[-2,6,-2,6],None,[])

    def test_orange_branches_do_not_share_a_positive_length_trunk(self):
        paths=[path('R','r','a',(2,2),(4,2),1),path('R','r','b',(2,2),(4,3),1)]
        xy={'r':(2,2),'a':(4,2),'b':(4,3)}
        exceptions=route_paths(paths,xy,[0,6,0,6],None)
        self.assertEqual(validate_paths(paths,xy,[0,6,0,6],None,exceptions)['orange_overlap_grid'],0)

    def test_circle_clearance_checks_segment_interior(self):
        region=dict(shape='circle',center=[4,4],radius=2,node_ids=[])
        self.assertEqual(segment_clearance([0,4],[8,4],region),0.)
        self.assertAlmostEqual(segment_clearance([0,1],[8,1],region),1.)
        p=path('O','a','b',(1,2),(7,2),1)
        with self.assertRaisesRegex(ValueError,'clearance'):validate_paths([p],{'a':(1,2),'b':(7,2)},[0,8,0,8],region,[])

    def test_non_grid_lane_is_rejected(self):
        p=path('O','a','b',(1,.15),(4,.15),1)
        with self.assertRaisesRegex(ValueError,'Non-grid'):validate_paths([p],{'a':(1,.15),'b':(4,.15)},[0,6,0,6],None,[])

    def test_auxiliary_and_cofactor_flags_cannot_hide_diagonal_segments(self):
        for flag in ('auxiliary_lane','local_cofactor'):
            p=path('R','a','b',(.25,1),(2,3),0);p[flag]=True
            with self.assertRaisesRegex(ValueError,'Non-orthogonal'):
                validate_paths([p],{'a':(.25,1),'b':(2,3)},[0,6,0,6],None,[])

    def test_short_uncrossed_route_keeps_red_crossing_forbidden(self):
        paths=[path('R','a','b',(2,2),(3,2),2),path('O','c','d',(1,1),(4,3),1)]
        xy={k:p for route in paths for k,p in [(route['source'],route['points'][0]),(route['target'],route['points'][-1])]}
        exceptions=route_paths(paths,xy,[0,6,0,6],None)
        self.assertFalse(exceptions)
        self.assertEqual(validate_paths(paths,xy,[0,6,0,6],None,exceptions)['metrics']['red_orange'],0)

    def test_auxiliary_grid_can_audit_crossing_when_strict_route_is_absent(self):
        red=path('RED','a','b',(3,0),(3,6),2);red['red_fixed']=True
        orange=path('ORANGE','c','d',(1.25,1.25),(5,1.25),1)
        xy={'a':(3,0),'b':(3,6),'c':(1.25,1.25),'d':(5,1.25)}
        exceptions=route_paths([red,orange],xy,[0,6,0,6],None)
        self.assertEqual(len(exceptions),1)
        self.assertTrue(exceptions[0]['no_uncrossed_path'])
        self.assertEqual(validate_paths([red,orange],xy,[0,6,0,6],None,exceptions)['metrics']['red_orange'],1)

    def test_square_fanout_validation_survives_straight_segment_merging(self):
        region=dict(shape='square',left=0,top=5,side_length=5,node_ids=['ring'])
        route=path('R','other','ring',(9,0),(5,3),1);xy={'ring':(5,3),'other':(9,0)}
        self.assertTrue(incident_exempt((8,2),(5.25,2),route,xy,region))
        self.assertFalse(incident_exempt((8,0),(5.25,0),route,xy,region))

    def test_circle_fanout_is_local_and_cannot_cross_the_cycle(self):
        region=dict(shape='circle',center=[0,0],radius=4,node_ids=['ring'])
        route=path('R','other','ring',(8,1),(4,0),1);xy={'ring':(4,0),'other':(8,1)}
        self.assertTrue(incident_exempt((4.25,.5),(7,.5),route,xy,region))
        self.assertFalse(incident_exempt((4.25,2),(7,2),route,xy,region))
        self.assertFalse(incident_exempt((3.75,.5),(7,.5),route,xy,region))

    def test_locked_corridor_detour_is_not_a_feasible_uncrossed_candidate(self):
        red=path('RED','a','b',(1.75,2),(2.25,2),2);red['red_fixed']=True
        orange=path('EX','c','d',(2,0),(2,4),1);orange['straight_locked']=True
        xy={'a':(1.75,2),'b':(2.25,2),'c':(2,0),'d':(2,4)}
        exceptions=route_paths([red,orange],xy,[0,5,0,5],None)
        self.assertEqual(len(exceptions),1)
        self.assertTrue(exceptions[0]['no_uncrossed_path'])
        self.assertEqual(orange['points'],[(2.,0.),(2.,4.)])

    def test_gray_locked_corridor_is_reserved_before_free_orange_routes(self):
        gray=path('T','t','c',(5,2),(1,2),0);gray['straight_locked']=True
        orange=path('O','r','c',(3,3),(1,2),1)
        xy={'t':(5,2),'c':(1,2),'r':(3,3)}
        route_paths([orange,gray],xy,[0,7,0,6],None)
        self.assertEqual(gray['points'],[(5.,2.),(1.,2.)])
        self.assertFalse(parallel_conflicts([gray,orange]))

    def test_shared_node_ports_do_not_consume_thin_corridor_capacity(self):
        endpoints={'a':(4,2),'b':(4,3),'c':(4,1),'d':(4,4)};xy=dict(endpoints,r=(2,2))
        points=[[(2,2),(4,2)],[(2,2),(2.25,2),(2.25,3),(4,3)],
                [(2,2),(2.25,2),(2.25,1),(4,1)],
                [(2,2),(2.25,2),(2.25,2.25),(2.75,2.25),(2.75,4),(4,4)]]
        paths=[]
        for key,vertices in zip(endpoints,points):
            p=path(key,'r',key,vertices[0],vertices[-1],1);p['points']=vertices
            if key!='a':p['auxiliary_lane']=True
            paths.append(p)
        report=validate_paths(paths,xy,[0,6,0,6],None,[])
        self.assertLessEqual(report['maximum_thin_corridor_load'],2)
        self.assertEqual(report['orange_overlap_grid'],0)

    def test_short_thin_corridor_allows_multiple_lines_below_two_grids(self):
        paths=[path(str(i),'a'+str(i),'b'+str(i),(x,1),(x,2.5),1) for i,x in enumerate([2.,2.25,2.75])]
        for p in paths[1:]:p['auxiliary_lane']=True
        xy={k:p for route in paths for k,p in [(route['source'],route['points'][0]),(route['target'],route['points'][-1])]}
        self.assertEqual(validate_paths(paths[:2],xy,[0,6,0,6],None,[])['maximum_thin_corridor_load'],2)
        report=validate_paths(paths,xy,[0,6,0,6],None,[])
        self.assertEqual(report['maximum_thin_corridor_load'],3)
        self.assertIsNone(report['thin_lane_limit'])
        for p in paths:p['points'][-1]=(p['points'][-1][0],3.)
        with self.assertRaisesRegex(ValueError,'Long close parallel'):validate_paths(paths,xy,[0,6,0,6],None,[])

    def test_parallel_limit_includes_two_grids_and_all_colors(self):
        for priority in (0,1,2):
            first=path('A','a','b',(1,1),(5,1),priority)
            second=path('B','c','d',(2,1.75),(4,1.75),0)
            self.assertEqual(len(parallel_conflicts([first,second])),1)
            second['points']=[(2,2),(4,2)]
            self.assertFalse(parallel_conflicts([first,second]))
            second['points']=[(2,1.75),(3.99,1.75)]
            self.assertFalse(parallel_conflicts([first,second]))
        first['points']=[(1,1),(2,1),(3,1),(5,1)]
        second['points']=[(2,1),(4,1)]
        self.assertTrue(parallel_conflicts([first,second]))

    def test_search_reroutes_long_close_parallel_gray_line(self):
        paths=[path('A','a','b',(1,2),(5,2),1),path('B','c','d',(2,2.25),(5,4),0)]
        xy={'a':(1,2),'b':(5,2),'c':(2,2.25),'d':(5,4)}
        route_paths(paths,xy,[0,7,0,7],None)
        self.assertFalse(parallel_conflicts(paths))

    def test_grid_exempts_cofactors_and_tca_only(self):
        nodes={'ring':dict(kind='metabolite'),'cofactor':dict(kind='currency'),'ordinary':dict(kind='reaction')}
        xy={'ring':(1.2,2),'cofactor':(3.2,4.3),'ordinary':(3,4)}
        self.assertTrue(validate_node_grid(nodes,xy,dict(node_ids=['ring']))['verified'])
        xy['ordinary']=(3.5,4)
        with self.assertRaisesRegex(ValueError,'grid intersections'):validate_node_grid(nodes,xy,dict(node_ids=['ring']))

    def test_cofactors_can_leave_both_grid_lines(self):
        nodes={'r':dict(kind='reaction')}
        nodes.update({k:dict(kind='currency',parent_id='r') for k in ('a','b','c')})
        xy={'r':[5.,5.],'a':[5,6],'b':[5,7],'c':[5,8]}
        attach_grid_cofactors(nodes,xy,None,[0,10,0,10])
        self.assertTrue(all(not any(integer(v) for v in xy[k]) for k in ('a','b','c')))
        self.assertTrue(all(math.dist(xy[k],xy['r'])<1 for k in ('a','b','c')))

    def test_parallel_run_cannot_reset_by_retracing_a_short_spur(self):
        xy={'r':(3,44),'pep':(3,27),'pyr':(3,25),'helper':(3.45,43.25)}
        saved=path('GLCpts','r','pep',xy['r'],xy['pep'],2)
        saved.update(points=[(3,44),(3.25,44),(3.25,27),(3,27)],auxiliary_lane=True)
        following=path('GLCpts','r','pyr',xy['r'],xy['pyr'],2)
        grid=Grid(xy,[0,8,20,46],None,auxiliary=True);grid.reserve(saved)
        points,length=grid.search(following)
        following['points']=points
        self.assertFalse(parallel_conflicts([saved,following]))

    def test_shared_connection_is_confined_to_node_glyph(self):
        first=path('A','r','a',(2,2),(4,3),1);second=path('B','r','b',(2,2),(4,1),1)
        first.update(points=[(2,2),(2.25,2),(2.25,3),(4,3)],auxiliary_lane=True)
        second.update(points=[(2,2),(2.25,2),(2.25,1),(4,1)],auxiliary_lane=True)
        xy={'r':(2,2),'a':(4,3),'b':(4,1)}
        report=validate_paths([first,second],xy,[0,6,0,6],None,[])
        self.assertEqual(report['orange_overlap_grid'],0.)
        self.assertEqual(report['raw_orange_overlap_grid'],.25)
        self.assertEqual(report['node_port_nodes'],['r'])
        first['points']=[(2,2),(3,2),(3,3),(4,3)]
        second['points']=[(2,2),(3,2),(3,1),(4,1)]
        with self.assertRaisesRegex(ValueError,'Visible orange'):validate_paths([first,second],xy,[0,6,0,6],None,[])

    def test_tca_attachment_cannot_reenter_opposite_clearance_band(self):
        region=dict(shape='circle',center=[0,0],radius=2,node_ids=['ring'])
        route=path('R','ring','other',(2,0),(0,0),1);xy={'ring':(2,0),'other':(0,0)}
        self.assertTrue(incident_exempt((2,0),(0,0),route,xy,region))
        self.assertFalse(incident_exempt((2,0),(-2,0),route,xy,region))

    def test_compact_circle_anchors_have_one_integer_coordinate(self):
        base=Path(__file__).parent;model=load_model(base/'e_coli_core.xml' if (base/'e_coli_core.xml').exists() else base/'models/e_coli_core.xml');strategy=BacterialStrategy(model)
        keys={k for route in strategy.routes for mid in route['metabolite_ids'] for k in ['met:'+mid]}
        keys.update('rxn:'+r for route in strategy.routes for r in route['reaction_ids'])
        fixed,b=strategy.plan(strategy.routes,keys);fixed,b=configure_plan(fixed,b,'circle','allow');old=b['tca_geometry']['radius']
        fixed,b=compact_tca(fixed,b);region=b['tca_geometry'];self.assertAlmostEqual(region['radius']/old,2/3,delta=.005)
        for k in region['node_ids']:
            self.assertTrue(any(integer(v) for v in fixed[k]))
            self.assertAlmostEqual(math.dist(fixed[k],region['center_offset']),region['radius'])
        self.assertEqual(len(set(tuple(fixed[k]) for k in region['node_ids'])),18)
        self.assertEqual(parser().parse_args([]).crossing_policy,'adaptive-grid')
