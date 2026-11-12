"""Lattice capacity, circular tangents and restricted curved map incidences."""
import copy
import math
import unittest
from export_escher import create_map, validate_map
from escher_routing import separate_routes
from reaction_network import build_network
from tca_geometry import lattice_circle
from tca_curves import arc_controls, sample_curve


class TCACurveTests(unittest.TestCase):
    def test_lattice_capacity_and_clockwise_order(self):
        for count in (8,12,16,20,24,28):
            points,radius=lattice_circle(count)
            self.assertGreaterEqual(len(points),count)
            self.assertEqual(len(points),len(set(points)))
            origin=points[0]
            for point in points:
                self.assertAlmostEqual(math.hypot(*point),radius)
                for axis in (0,1):self.assertAlmostEqual(point[axis]-origin[axis],round(point[axis]-origin[axis]))
            self.assertLess(sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(points,points[1:]+points[:1])),0)

    def test_curve_endpoints_tangents_and_reverse(self):
        start,end,center=[3,4],[4,3],[0,0]
        controls=arc_controls(start,end,center)
        points=sample_curve(start,controls,end)
        self.assertEqual(points[0],start);self.assertEqual(points[-1],end)
        self.assertLess(max(abs(math.dist(p,center)-5) for p in points),.001)
        for point,control in zip((start,end),controls):
            self.assertAlmostEqual(sum((point[i]-center[i])*(control[i]-point[i]) for i in (0,1)),0)
        reverse=arc_controls(end,start,center)
        for actual,expected in zip(reverse,reversed(controls)):
            for a,b in zip(actual,expected):self.assertAlmostEqual(a,b)

    def make_map(self,shape='circle'):
        model=dict(id='test',metabolites=[dict(id=k,name=k,compartment='c') for k in ('a_c','b_c','other_c')],
                   reactions=[dict(id='R',name='R',metabolites={'a_c':-1,'b_c':1,'other_c':1},lower_bound=-10,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'met:a_c':[0,5],'rxn:R':[3,4],'met:b_c':[4,3],'met:other_c':[6,6]}
        data=create_map(model,nodes,edges,xy,scale=180,padding=80)
        backbone=dict(tca_geometry=dict(shape=shape),tca_exclusion=dict(region=dict(center=[0,0],radius=5)),
                      routes=[dict(id='tca',metabolite_ids=['a_c','b_c'],steps=[dict(reaction_id='R',source='a_c',target='b_c')])])
        data[0]['metadata']['sbml2map']=dict(major_backbones=backbone)
        return model,nodes,xy,separate_routes(data,['R'])

    def test_curves_only_follow_supported_tca_steps_and_keep_chemistry(self):
        model,nodes,xy,data=self.make_map()
        report=validate_map(data,model,nodes,xy)
        self.assertEqual(report['tca_curved_segments'],2)
        self.assertTrue(report['non_tca_orthogonal_segments_verified'])
        for reaction in data[1]['reactions'].values():
            for sid,seg in reaction['segments'].items():
                if seg['b1'] is not None:
                    self.assertNotIn('met:other_c',data[0]['metadata']['segment_source_edges'][sid])
        broken=copy.deepcopy(data)
        seg=next(s for r in broken[1]['reactions'].values() for s in r['segments'].values() if s['b1'] is not None)
        seg['b1']['x']+=1
        with self.assertRaisesRegex(ValueError,'valid TCA'):validate_map(broken,model,nodes,xy)
        # Idempotent regeneration, including reversal of substrate incidences.
        again=separate_routes(copy.deepcopy(data),['R'])
        self.assertEqual(again,data)

    def test_square_keeps_all_connections_orthogonal(self):
        model,nodes,xy,data=self.make_map('square')
        report=validate_map(data,model,nodes,xy)
        self.assertTrue(report['orthogonal_segments_verified'])
        self.assertEqual(report['tca_curved_segments'],0)
