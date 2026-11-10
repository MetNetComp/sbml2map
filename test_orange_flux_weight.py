"""Flux-weighted connectors and placement share the same positive objective."""
import itertools
import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from exchange_routes import orange_flux_costs, connector_plan, connector_length, connector_cost
from grid_layout2 import optimize_grid, edge_length
from sbml2map import parser
from escher_grid import validate_sampling_options


class OrangeFluxTests(unittest.TestCase):
    def test_cli_default_and_explicit_weights(self):
        self.assertEqual(parser().parse_args([]).orange_flux_weight, 1)
        for weight in (0, .5, 1, 2):
            args=parser().parse_args(['--orange-flux-weight',str(weight)])
            self.assertEqual(args.orange_flux_weight, weight)

    def fixture(self):
        nodes={key:dict(kind=kind) for key,kind in [('root','metabolite'),('low','reaction'),
              ('high','reaction'),('b','metabolite'),('EX','reaction'),('EX2','reaction')]}
        edges=[dict(source=a,target=b,reaction_id=r) for a,b,r in [
            ('root','low','low'),('low','b','low'),('root','high','high'),('high','b','high'),
            ('b','EX','EX'),('b','EX2','EX2')]]
        xy=dict(root=[0,0],low=[2,0],high=[2,1],b=[4,0],EX=[5,0],EX2=[4,-1])
        records=[dict(reaction_id=rid,median=median) for rid,median in [('low',.1),('high',-10),('EX',1),('EX2',1)]]
        return nodes,edges,xy,records

    def test_longer_high_flux_route_can_win_and_shared_edges_count_once(self):
        nodes,edges,xy,records=self.fixture()
        weights=orange_flux_costs(nodes,edges,set(),records,1)['reaction_weights']
        plain=connector_plan(nodes,edges,{'root'},{'EX','EX2'},set(),xy)
        weighted=connector_plan(nodes,edges,{'root'},{'EX','EX2'},set(),xy,reaction_weights=weights)
        self.assertEqual(plain['edge_indices'],[0,1,4,5])
        self.assertEqual(weighted['edge_indices'],[2,3,4,5])
        self.assertEqual(connector_length(plain,edges,xy),6)
        self.assertEqual(connector_length(weighted,edges,xy),8)
        self.assertLess(connector_cost(weighted,edges,xy),connector_cost(dict(plain,reaction_weights=weights),edges,xy))
        self.assertTrue(weighted['all_exchanges_connected'])
        self.assertEqual(len(weighted['edge_indices']),len(set(weighted['edge_indices'])))

    def test_midrank_ties_sign_and_zero_weight(self):
        nodes,edges,xy,records=self.fixture()
        result=orange_flux_costs(nodes,edges,set(),records,1)
        self.assertEqual(result['normalized_flux']['low'],0)
        self.assertEqual(result['normalized_flux']['high'],1)
        self.assertEqual(result['normalized_flux']['EX'],result['normalized_flux']['EX2'])
        self.assertEqual(result['reaction_weights']['high'],1)
        self.assertEqual(result['reaction_weights']['low'],2)
        weights=orange_flux_costs(nodes,edges,set(),[],0)['reaction_weights']
        plain=connector_plan(nodes,edges,{'root'},{'EX'},set(),xy)
        zero=connector_plan(nodes,edges,{'root'},{'EX'},set(),xy,reaction_weights=weights)
        self.assertEqual(plain['edge_indices'],zero['edge_indices'])
        self.assertEqual(connector_length(zero,edges,xy),connector_cost(zero,edges,xy))

    def test_invalid_strength_and_missing_medians(self):
        nodes,edges,_,records=self.fixture()
        for invalid in (-1,float('nan'),float('inf')):
            with self.assertRaisesRegex(ValueError,'finite and nonnegative'):
                orange_flux_costs(nodes,edges,set(),records,invalid)
            args=parser().parse_args(['--orange-flux-weight',str(invalid)])
            with self.assertRaisesRegex(ValueError,'finite and nonnegative'):
                validate_sampling_options(args)
        with self.assertRaisesRegex(ValueError,'median missing'):
            orange_flux_costs(nodes,edges,set(),[],1)
        records[0]['median']=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            orange_flux_costs(nodes,edges,set(),records,1)

    def test_major_paths_are_free_and_do_not_enter_flux_normalization(self):
        nodes,edges,xy,records=self.fixture()
        result=orange_flux_costs(nodes,edges,{'low'},records,1)
        self.assertNotIn('low',result['reaction_weights'])
        plan=connector_plan(nodes,edges,{'root'},{'EX'}, {'low'},xy,reaction_weights=result['reaction_weights'])
        self.assertEqual(plan['edge_indices'],[4])


class WeightedGridTests(unittest.TestCase):
    def test_weighted_milp_matches_exhaustive_grid_search(self):
        ids=['a','b','c','d'];indices=np.array(list(itertools.combinations(range(4),2)))
        pairs=[(ids[a],ids[b]) for a,b in indices];weights=np.array([1.,1.3,2.,1.8,1.1,1.9])
        cells=[(x,y) for x in range(3) for y in range(2)]
        exact=min(edge_length(np.array(points),indices,weights) for points in itertools.permutations(cells,4))
        xy,report=optimize_grid(ids,pairs,3,2,time_limit=5,edge_weights=weights,progress=lambda *a,**k:None)
        actual=edge_length(np.array([xy[k] for k in ids]),indices,weights)
        self.assertAlmostEqual(report['objective'],exact)
        self.assertAlmostEqual(actual,exact)
        self.assertTrue(report['optimality_proven'])
        self.assertLessEqual(report['global_lower_bound'],actual+1e-7)

    def test_heuristic_and_solver_failure_report_weighted_cost(self):
        ids=['a','b','c'];pairs=[('a','b'),('b','c')];weights=[1.4,2.]
        for kwargs in ({'max_milp_nodes':1},{}):
            with patch('grid_layout2.milp',return_value=SimpleNamespace(x=None,status=1,message='limit',mip_dual_bound=None)):
                xy,report=optimize_grid(ids,pairs,3,2,edge_weights=weights,progress=lambda *a,**k:None,**kwargs)
            physical=sum(abs(xy[a][axis]-xy[b][axis])*w for (a,b),w in zip(pairs,weights) for axis in (0,1))
            self.assertAlmostEqual(report['objective'],physical)
            self.assertAlmostEqual(report['global_lower_bound'],sum(weights))
        for invalid in ([1],[-1,1],[0,1],[1,float('inf')]):
            with self.assertRaisesRegex(ValueError,'finite, positive'):
                optimize_grid(ids,pairs,3,2,edge_weights=invalid)


if __name__=='__main__':unittest.main()
