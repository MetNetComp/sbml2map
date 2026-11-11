"""Crossing/bend counts describe displayed geometry, including corner contacts."""
import copy
import unittest
from routing_metrics import measure, MetricIndex, objective, validate_options


def path(rid, points, source=None, target=None, **kwargs):
    return dict(rid=rid, points=points, source=source or rid+':a', target=target or rid+':b', priority=0, **kwargs)


class MetricsTests(unittest.TestCase):
    def test_crossing_at_bend_and_unrelated_endpoint_counts_once(self):
        for points in ([(-1,0),(1,0)], [(-1,0),(0,0),(1,0)], [(-1,0),(0,0),(0,1)]):
            m=measure([path('a',points),path('b',[(0,-1),(0,0)])])
            self.assertEqual(m['crossings'],1)

    def test_common_biological_endpoint_is_not_crossing(self):
        paths=[path('a',[(-1,0),(0,0)],target='met:x'),path('b',[(0,0),(0,1)],source='met:x')]
        self.assertEqual(measure(paths)['crossings'],0)
        paths[1]['source']='met:y'
        self.assertEqual(measure(paths)['crossings'],1)

    def test_shared_corner_and_trunk_count_once_and_ignore_segmentation(self):
        paths=[path('a',[(0,0),(2,0),(2,2)]),path('a',[(0,0),(1,0),(2,0),(2,3)])]
        m=measure(paths)
        self.assertEqual(m['bends'],1)
        self.assertEqual(m['total_length'],5)
        self.assertEqual(m['crossings'],0)

    def test_coincident_unrelated_lines_are_explicit_violations(self):
        m=measure([path('a',[(0,0),(3,0)]),path('b',[(1,0),(2,0)])])
        self.assertEqual(m['coincident_intervals'],1)

    def test_reversal_rotation_and_scaling_preserve_counts(self):
        paths=[path('a',[(-1,0),(1,0),(1,2)]),path('b',[(0,-1),(0,1)])]
        expected=measure(paths)
        for reverse in (False,True):
            new=copy.deepcopy(paths)
            for p in new:
                p['points']=[(10*y,10*x) for x,y in p['points']]
                if reverse:p['points'].reverse();p['source'],p['target']=p['target'],p['source']
            actual=measure(new)
            self.assertEqual(actual['crossings'],expected['crossings'])
            self.assertEqual(actual['bends'],expected['bends'])
            self.assertAlmostEqual(actual['total_length'],10*expected['total_length'])

    def test_incremental_metrics_equal_full_recomputation_for_two_reactions(self):
        paths=[path('a',[(-1,0),(1,0)]),path('b',[(0,-1),(0,1)]),path('c',[(2,0),(3,0)])]
        index=MetricIndex(paths)
        changed={'a':[path('a',[(-1,0),(-1,2),(1,2),(1,0)])],'b':[path('b',[(0,-1),(2,-1),(2,1),(0,1)])]}
        actual,events,geometries=index.evaluate(changed)
        expected=measure(changed['a']+changed['b']+[paths[2]])
        self.assertEqual(actual,expected)
        index.commit(changed,events,geometries)
        self.assertEqual(index.metrics(),expected)

    def test_crossings_first_never_trades_a_crossing_for_fewer_bends(self):
        a=dict(crossings=0,bends=100,crossings_plus_bends=100,total_length=100)
        b=dict(crossings=1,bends=0,crossings_plus_bends=1,total_length=1)
        self.assertLess(objective(a),objective(b))
        self.assertGreater(objective(a,'sum'),objective(b,'sum'))
        self.assertGreater(objective(a,'weighted',10,1),objective(b,'weighted',10,1))
        self.assertLess(objective(a,'weighted',1000,1),objective(b,'weighted',1000,1))

    def test_invalid_objective_weights_and_budget_are_rejected(self):
        for mode,c,b,t in [('bad',1,1,1),('sum',-1,1,1),('sum',1,float('nan'),1),('sum',0,0,1),('sum',1,1,-1)]:
            with self.assertRaises(ValueError):validate_options(mode,c,b,t)
