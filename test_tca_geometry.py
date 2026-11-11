"""Exact circular anchors and hard empty-interior constraints."""
import copy
import math
from pathlib import Path
import unittest
import numpy as np
from sbml2map import parser,BacterialStrategy
from model_input import load_model
from major_backbones2 import validate_backbones
from tca_geometry import configure_plan,absolute_region,allowed_points,validate_exclusion,on_perimeter
from escher_grid import attach_coenzymes
from grid_layout2 import optimize_grid

BASE=Path(__file__).resolve().parent
MODEL=BASE/'models/e_coli_core.xml' if (BASE/'models').is_dir() else BASE/'e_coli_core.xml'


class TCAGeometryTests(unittest.TestCase):
    def test_defaults_and_explicit_choices(self):
        args=parser().parse_args([])
        self.assertEqual(args.tca_shape,'circle')
        self.assertEqual(args.exchange_alignment,'collinear')
        args=parser().parse_args(['--tca-shape','square','--exchange-alignment','off'])
        self.assertEqual((args.tca_shape,args.exchange_alignment),('square','off'))
        self.assertFalse(hasattr(args,'routing_objective'))

    def test_real_cycle_is_clockwise_on_requested_perimeter(self):
        model=load_model(MODEL);strategy=BacterialStrategy(model)
        keys={k for route in strategy.routes for mid in route['metabolite_ids'] for k in ['met:'+mid]}
        keys.update('rxn:'+r for route in strategy.routes for r in route['reaction_ids'])
        for shape in ('circle','square'):
            fixed,plan=strategy.plan(strategy.routes,keys)
            fixed,plan=configure_plan(fixed,plan,shape)
            rows=50;xy={k:[x,rows-1+y] for k,(x,y) in fixed.items()}
            region=absolute_region(plan['tca_geometry'],rows)
            self.assertTrue(validate_backbones(plan,xy,rows)['all_displayed_backbone_positions_verified'])
            for key in region['node_ids']:self.assertTrue(on_perimeter(xy[key],region))
            if shape=='circle':self.assertTrue(any(not float(v).is_integer() for key in region['node_ids'] for v in xy[key]))
            self.assertEqual(next(r for r in plan['routes'] if r['id']=='tca')['shape'],shape)

    def test_interiors_and_clearance_are_forbidden(self):
        circle={'shape':'circle','center':[3,3],'radius':2,'node_ids':['ring']}
        square={'shape':'square','left':1,'top':5,'side_length':4,'node_ids':['ring']}
        for region in (circle,square):
            self.assertFalse(allowed_points([[3,3]],region)[0])
            self.assertFalse(allowed_points([[5.1,3]],region)[0])
            self.assertTrue(allowed_points([[5.2,3]],region)[0])
            with self.assertRaisesRegex(RuntimeError,'TCA'):validate_exclusion({}, {'other':[3,3]},region)

    def test_tca_local_cofactors_are_outside_the_ring(self):
        region={'shape':'circle','center':[3,3],'radius':2,'node_ids':['rxn:tca']}
        nodes={'rxn:tca':dict(kind='reaction')}
        nodes.update({f'copy:{i}':dict(kind='currency',parent_id='rxn:tca') for i in range(8)})
        xy={'rxn:tca':[5,3]}
        attach_coenzymes(nodes,xy,grid_shape=(8,8),forbidden_region=region)
        self.assertTrue(validate_exclusion(nodes,xy,region)['verified'])
        for key in nodes:
            if key!='rxn:tca':self.assertAlmostEqual(sum(abs(a-b) for a,b in zip(xy[key],xy['rxn:tca'])),.45)

    def test_both_solver_modes_enforce_allowed_cells_and_fractional_alignment(self):
        radius=1.;angle=.7;anchor=[3+math.cos(angle),3+math.sin(angle)]
        region={'shape':'circle','center_offset':[3,-3],'radius':radius,'node_ids':['a_c']}
        group=dict(reaction_id='EX_a',exchange='EX_a',external='a_e',internal='a_c')
        for budget in (1,None):
            xy,report=optimize_grid(['EX_a','a_e','a_c','other'],[('EX_a','a_e'),('a_e','a_c')],7,7,
                fixed_offsets={'a_c':[anchor[0],anchor[1]-6]},boundary_nodes=['EX_a'],exclusive_boundary=True,
                reserved_region=region,alignment_groups=[group],max_milp_nodes=budget,time_limit=1,
                progress=lambda *a,**k:None)
            self.assertTrue(report['position_constraints']['exchange_alignment']['verified'])
            self.assertTrue(validate_exclusion({},xy,absolute_region(region,7))['verified'])
            self.assertTrue(np.allclose(xy['a_c'],anchor))
            if budget is None:self.assertIn('MILP',report['method'])

    def test_explicit_impossible_grid_fails_without_relaxing_constraints(self):
        region={'shape':'square','left':0,'top_offset':0,'side_length':4,'node_ids':[]}
        with self.assertRaisesRegex(ValueError,'Explicit grid'):
            optimize_grid(['n'],[],3,3,reserved_region=region,exclusive_boundary=True,time_limit=0,progress=lambda *a,**k:None)

    def test_full_fba_added_cofactor_stays_outside_and_preserves_saved_nodes(self):
        from full_fba_layout import complete_layout
        from reaction_network import build_network
        model=dict(id='fba',metabolites=[dict(id=mid,compartment='c',formula='C3H6O3') for mid in ('a_c','b_c','atp_c')],
            reactions=[dict(id='TCA',metabolites={'a_c':-1,'b_c':1,'atp_c':-1},lower_bound=0,upper_bound=10)])
        all_nodes,_=build_network(model,{'atp_c'})
        saved={k:v for k,v in all_nodes.items() if v['kind']!='currency'}
        xy={'rxn:TCA':[4,3],'met:a_c':[3,4],'met:b_c':[6,3]}
        region={'shape':'circle','center':[3,3],'center_offset':[3,-4],'radius':1,'node_ids':['rxn:TCA','met:a_c']}
        nodes,edges,positions,report=complete_layout(model,saved,xy,region,grid_shape=(8,8))
        self.assertFalse(report['existing_coordinates_changed'])
        self.assertEqual(report['tca_interior']['interior_node_count'],0)
        self.assertTrue(any(n['kind']=='currency' for n in nodes.values()))
        self.assertTrue(validate_exclusion(nodes,positions,region)['verified'])

    def test_no_tca_does_not_create_a_region(self):
        fixed,plan=configure_plan({},dict(routes=[],square=None),'circle')
        self.assertIsNone(plan['tca_geometry'])
