"""Crossing constraints include touching, overlaps and circular Beziers."""
import copy
import unittest
from crossing_routing import measure,route,orange_length,validate_boundary
from tca_curves import arc_controls


def path(rid,a,b,start,end,priority=2,points=None):
    return dict(rid=rid,mid=b,source=a,target=b,priority=priority,backbone=False,points=points or [start,end])


class CrossingRoutingTests(unittest.TestCase):
    def test_unrelated_crossings_contacts_and_overlaps_are_hard(self):
        a=path('A','a','b',(0,1),(2,1))
        for p in [path('B','c','d',(1,0),(1,2)),path('B','c','d',(1,1),(1,2)),path('B','c','d',(1,1),(3,1))]:
            self.assertFalse(measure([a,p])['hard_verified'])
        p['priority']=1
        self.assertGreater(measure([a,p])['red_orange'],0)

    def test_common_biological_endpoint_and_connected_reaction_trunk(self):
        a=path('R','r','a',(0,0),(2,0))
        b=path('S','a','b',(2,0),(2,1))
        self.assertTrue(measure([a,b])['hard_verified'])
        c=path('R','r','c',(0,0),(1,1),points=[(0,0),(1,0),(1,1)])
        self.assertTrue(measure([a,c])['hard_verified'])
        # Same reaction alone cannot excuse a disconnected crossing.
        c=path('R','r','c',(0,0),(2,2),points=[(0,0),(0,-1),(3,-1),(3,1),(1,1),(1,2),(2,2)])
        a=path('R','r','a',(0,0),(2,0),points=[(0,0),(2,0)])
        self.assertTrue(measure([a,c])['hard_verified'])
        c['points']=[(0,0),(0,-1),(1,-1),(1,2),(2,2)]
        self.assertFalse(measure([a,c])['hard_verified'])

    def test_curve_interior_and_tangency_are_detected(self):
        curve=path('TCA','a','b',(1,0),(0,1));curve['curve_controls']=arc_controls((1,0),(0,1),(0,0))
        line=path('B','c','d',(-1,.5),(2,.5),priority=1)
        self.assertGreater(measure([curve,line])['red_orange'],0)
        line['points']=[(-1,1),(2,1)]
        self.assertGreater(measure([curve,line])['red_orange'],0)

    def test_red_and_orange_reroute_without_crossings(self):
        for priority in (1,2):
            paths=[path('A','a','b',(0,1),(2,1)),path('B','c','d',(1,0),(1,2),priority=priority)]
            xy={k:p for q in paths for k,p in [(q['source'],q['points'][0]),(q['target'],q['points'][-1])]}
            before=copy.deepcopy(xy)
            result=route(paths,xy,attempts=2,progress=lambda *a,**k:None)
            self.assertTrue(result['verified']);self.assertTrue(measure(paths)['hard_verified'])
            self.assertEqual(xy,before)
            if priority==1:self.assertEqual(result['metrics']['orange_orange'],0)

    def test_closed_curve_is_a_real_barrier(self):
        points=[(1,0),(0,1),(-1,0),(0,-1)];paths=[];xy={}
        for i,(a,b) in enumerate(zip(points,points[1:]+points[:1])):
            p=path('R'+str(i),'n'+str(i),'n'+str((i+1)%4),a,b);p['curve_controls']=arc_controls(a,b,(0,0));paths.append(p)
            xy[p['source']]=a;xy[p['target']]=b
        p=path('X','inside','outside',(0,0),(2,0),priority=1);paths.append(p);xy.update(inside=(0,0),outside=(2,0))
        result=route(paths,xy,attempts=0,progress=lambda *a,**k:None)
        self.assertFalse(result['verified'])

    def test_orange_crossing_is_reported_without_hard_failure(self):
        paths=[path('A','a','b',(0,1),(2,1),priority=1),path('B','c','d',(1,0),(1,2),priority=1)]
        r=measure(paths);self.assertTrue(r['hard_verified']);self.assertEqual(r['orange_orange'],1)


    def test_orange_crossing_does_not_trigger_a_detour(self):
        paths=[path('A','a','b',(0,1),(2,1),priority=1),path('B','c','d',(1,0),(1,2),priority=1)]
        xy={k:p for q in paths for k,p in [(q['source'],q['points'][0]),(q['target'],q['points'][-1])]}
        result=route(paths,xy,attempts=3,progress=lambda *a,**k:None)
        self.assertTrue(result['verified'])
        self.assertFalse(result['orange_crossings_optimized'])
        self.assertEqual(result['metrics']['orange_orange'],1)
        self.assertAlmostEqual(orange_length(paths),4)
        self.assertNotIn('orange crossings',result['objective_order'])

    def test_shared_orange_reaction_trunk_counts_once(self):
        paths=[path('R','r','a',(0,0),(2,0),priority=1),
               path('R','r','b',(0,0),(1,2),priority=1,points=[(0,0),(1,0),(1,2)])]
        self.assertAlmostEqual(orange_length(paths),4)
        paths.append(path('OTHER','c','d',(0,0),(2,0),priority=1))
        self.assertAlmostEqual(orange_length(paths),6)


    def test_boundary_covers_gray_lines_and_curve_extrema(self):
        gray=path('G','a','b',(0,.5),(1,.5),priority=0,points=[(0,.5),(-.1,.5),(-.1,.8),(1,.8),(1,.5)])
        with self.assertRaisesRegex(ValueError,'exchange boundary'):validate_boundary([gray],[0,1,0,1])
        curve=path('TCA','a','b',(.2,.5),(.8,.5))
        curve['curve_controls']=[(.2,2),(.8,2)]
        with self.assertRaisesRegex(ValueError,'exchange boundary'):validate_boundary([curve],[0,1,0,1])
        # Controls can lie outside while the actual curve stays inside.
        curve['curve_controls']=[(.2,1.1),(.8,1.1)]
        self.assertTrue(validate_boundary([curve],[0,1,0,1])['verified'])

    def test_bounded_detour_keeps_red_crossings_zero(self):
        paths=[path('A','a','b',(.5,1),(1.5,1)),path('B','c','d',(1,.5),(1,1.5),priority=1)]
        xy={k:p for q in paths for k,p in [(q['source'],q['points'][0]),(q['target'],q['points'][-1])]}
        result=route(paths,xy,boundary=[0,2,0,2],attempts=1,progress=lambda *a,**k:None)
        self.assertTrue(result['verified']);self.assertTrue(result['boundary_validation']['verified'])
        self.assertEqual(result['boundary_validation']['outside_path_count'],0)

    def test_bounded_infeasibility_never_escapes_outside(self):
        paths=[path('A','a','b',(0,1),(2,1)),path('B','c','d',(1,0),(1,2),priority=1)]
        xy={k:p for q in paths for k,p in [(q['source'],q['points'][0]),(q['target'],q['points'][-1])]}
        result=route(paths,xy,boundary=[0,2,0,2],attempts=0,progress=lambda *a,**k:None)
        self.assertFalse(result['verified'])
        validate_boundary(paths,[0,2,0,2])

    def test_exchange_boundary_entry_points_inward(self):
        paths=[path('EX','exchange','external',(1,2),(1,1))]
        xy=dict(exchange=(1,2),external=(1,1))
        result=route(paths,xy,boundary=[0,2,0,2],attempts=0,progress=lambda *a,**k:None)
        self.assertTrue(result['verified']);self.assertEqual(paths[0]['points'][1][0],1)
        self.assertLess(paths[0]['points'][1][1],2)

    def test_exchange_envelope_matches_routing_rectangle(self):
        from exchange_alignment import anchor_exchange_boundary,validate_alignment
        xy={'e0':[0,2],'a0':[1,2],'i0':[2,2], 'e1':[0,4],'a1':[1,4],'i1':[2,4],
            'e2':[6,3],'a2':[5,3],'i2':[4,3], 'e3':[3,6],'a3':[3,5],'i3':[3,4]}
        before=copy.deepcopy(xy);groups=[dict(reaction_id=str(i),exchange='e'+str(i),external='a'+str(i),internal='i'+str(i)) for i in range(4)]
        changes=anchor_exchange_boundary(groups,{'e0','e1','e2','e3'},xy,7,7)
        self.assertEqual(len(changes),1);self.assertEqual(changes[0]['face'],'bottom')
        self.assertTrue(validate_alignment(groups,xy,7,7)['verified'])
        for i in range(4):self.assertEqual(xy['i'+str(i)],before['i'+str(i)])
