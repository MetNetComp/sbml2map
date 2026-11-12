import copy
import math
import unittest

from red_objective import red_length,route_length,translation_seed,translated_core,preservation_report,reference_signature,relaxation_level,optimize_red_core
from test_adaptive_routing import path
from tca_curves import arc_controls,sample_curve


class RedObjectiveTests(unittest.TestCase):
    def test_actual_detour_length_and_color_selection(self):
        red=path('R','a','b',(0,0),(2,0),2)
        red['points']=[(0,0),(0,3),(2,3),(2,0)]
        orange=path('O','c','d',(0,0),(100,0),1)
        self.assertEqual(red_length([red,orange]),8)

    def test_cubic_length_is_not_chord_or_endpoint_l1(self):
        red=path('R','a','b',(1,0),(0,1),2)
        red['curve_controls']=arc_controls((1,0),(0,1),(0,0))
        actual=route_length(red)
        sampled=sample_curve((1,0),red['curve_controls'],(0,1),4096)
        expected=sum(math.dist(a,b) for a,b in zip(sampled,sampled[1:]))
        self.assertAlmostEqual(actual,expected,places=7)
        self.assertGreater(actual,math.sqrt(2));self.assertLess(actual,2)

    def test_milp_moves_tca_toward_red_attachment(self):
        red=path('R','fixed','ring',(0,0),(10,0),2)
        seed,report=translation_seed([red],{'fixed':(0,0),'ring':(10,0)},{'ring'},[(-8,0),(0,0)],1)
        self.assertEqual(seed,(-8,0));self.assertEqual(report['objective'],2)
        self.assertTrue(report['optimality_proven'])

    def test_rigid_translation_updates_all_coordinate_representations(self):
        region=dict(shape='circle',center=[5,4],center_offset=[5,-5],radius=1,node_ids=['ring'])
        b=dict(tca_exclusion=dict(region=region),tca_geometry=copy.deepcopy(region),
               circle=dict(center_offset=[5,-5]),fixed_node_offsets={'ring':[5,-4]},
               routes=[dict(node_offsets={'ring':[5,-4]})])
        xy={'ring':[5,5],'fixed':[1,5]};before=copy.deepcopy(b)
        positions,trial=translated_core(xy,b,{'ring'},(-1,2),10)
        self.assertEqual(positions['ring'],[4,7]);self.assertEqual(positions['fixed'],[1,5])
        self.assertEqual(trial['tca_exclusion']['region']['center'],[4,6])
        self.assertEqual(trial['circle']['center_offset'],[4,-3])
        self.assertEqual(trial['fixed_node_offsets']['ring'],[4,-2])
        self.assertEqual(b,before)

    def test_global_translation_keeps_lengths_and_the_top_boundary(self):
        b=dict(fixed_node_offsets={'r':[1,0],'a':[1,-2]},routes=[],tca_exclusion=dict(region=None))
        xy={'r':[1,9],'a':[1,7]}
        positions,trial=translated_core(xy,b,set(),(0,0),10,global_dx=4,red_keys=set(xy))
        self.assertEqual(positions,{'r':[5,9],'a':[5,7]})
        self.assertEqual(math.dist(*positions.values()),math.dist(*xy.values()))
        self.assertEqual(trial['fixed_node_offsets']['r'],[5,0])

    def test_red_stage_excludes_unplaced_nonred_obstacles(self):
        nodes={'rxn:R':dict(kind='reaction'),'met:a':dict(kind='metabolite'),
               'met:b':dict(kind='metabolite'),'other':dict(kind='reaction')}
        xy={'rxn:R':[1,1],'met:a':[3,1],'met:b':[1,3],'other':[2,1]}
        edges=[dict(reaction_id='R',metabolite_id=mid,source='rxn:R',target='met:'+mid) for mid in ('a','b')]
        b=dict(fixed_node_offsets={},routes=[],tca_exclusion=dict(region=None))
        cores,report=optimize_red_core(nodes,edges,xy,b,{'R'},10,10,2,2,lambda *a,**k:None)
        self.assertEqual(cores[0]['length'],4)
        self.assertNotIn('other',cores[0]['positions'])
        self.assertFalse(report['global_optimality_proven'])

    def test_ring_cofactors_use_distinct_long_ports(self):
        from adaptive_routing import attach_grid_cofactors,make_paths,parallel_conflicts,validate_paths
        parent='rxn:R';keys=['copy:R:q8','copy:R:q8h2']
        nodes={parent:dict(kind='reaction')}
        nodes.update({k:dict(kind='currency',parent_id=parent) for k in keys})
        xy={parent:[3.,3.],keys[0]:[9.,9.],keys[1]:[9.,8.]}
        region=dict(shape='circle',center=[0,0],radius=math.sqrt(18),node_ids=[parent],interior_policy='allow',clearance_grid=1.)
        attach_grid_cofactors(nodes,xy,region,[-10,10,-10,10])
        edges=[dict(reaction_id='R',metabolite_id=k,source=parent,target=k) for k in keys]
        b=dict(routes=[],tca_geometry=dict(shape='circle'),tca_exclusion=dict(region=region))
        paths=make_paths(nodes,edges,xy,b,{'R'})
        self.assertFalse(parallel_conflicts(paths))
        self.assertTrue(validate_paths(paths,xy,[-10,10,-10,10],region,[])['verified'])

    def test_snapshot_detects_position_route_and_length_changes(self):
        red=path('R','a','b',(0,0),(2,0),2);xy={'a':(0,0),'b':(2,0)}
        reference=dict(positions=copy.deepcopy(xy),routes=[copy.deepcopy(red)],length=2,
                       signature=reference_signature(xy,[red]))
        self.assertTrue(preservation_report(reference,xy,[red])['exactly_preserved'])
        red['points']=[(0,0),(0,1),(2,1),(2,0)]
        result=preservation_report(reference,xy,[red])
        self.assertEqual(result['changed_route_count'],1);self.assertEqual(result['length_increase_grid'],2)
        result=preservation_report(reference,dict(xy,b=(3,0)),[red])
        self.assertEqual(result['moved_node_ids'],['b'])
        self.assertEqual(reference['routes'][0]['points'],[(0,0),(2,0)])

    def test_relaxation_requires_equal_length_before_relocation_or_growth(self):
        red=path('R','a','b',(0,0),(2,2),2);red['points']=[(0,0),(0,2),(2,2)]
        xy={'a':[0,0],'b':[2,2]}
        reference=dict(positions=xy,routes=[red],length=4,signature='reference')
        self.assertEqual(relaxation_level(reference,reference),0)
        reroute=copy.deepcopy(reference);reroute['routes'][0]['points']=[(0,0),(2,0),(2,2)]
        self.assertEqual(relaxation_level(reference,reroute),1)
        relocate=copy.deepcopy(reference);relocate['positions']['a']=[1,0]
        relocate['routes'][0]['points']=[(1,0),(1,3),(2,3),(2,2)]
        self.assertEqual(relaxation_level(reference,relocate),3)
        relocate['routes'][0]['points']=[(1,0),(0,0),(0,2),(2,2)]
        relocate['length']=5
        self.assertEqual(relaxation_level(reference,relocate),3)
        relocate['positions']={'a':[1,0],'b':[3,2]}
        relocate['routes'][0]['points']=[(1,0),(1,2),(3,2)]
        self.assertEqual(relaxation_level(reference,relocate),2)


if __name__=='__main__':unittest.main()
