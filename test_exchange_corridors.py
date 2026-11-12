"""Straight exchange geometry is arranged after red routes are frozen."""
import copy
import unittest
from reaction_network import build_network
from exchange_corridors import discover_corridors,place_corridors,lock_corridor_paths,validate_corridors
from crossing_routing import initial_paths,route


class StraightExchangeTests(unittest.TestCase):
    def fixture(self,triplet=True,barrier=False):
        metabolites=[dict(id='a_e',compartment='e',name='a_e')]
        reactions=[dict(id='EX_a',name='EX_a',metabolites={'a_e':-1},lower_bound=-10,upper_bound=10)]
        if triplet:
            metabolites.append(dict(id='a_c',compartment='c',name='a_c'))
            reactions.append(dict(id='T',name='T',metabolites={'a_e':-1,'a_c':1},lower_bound=-10,upper_bound=10))
        model=dict(id='test',metabolites=metabolites,reactions=reactions)
        nodes,edges=build_network(model,[])
        xy={'rxn:EX_a':[0,3],'met:a_e':[1,3]}
        if triplet:xy.update({'met:a_c':[4,3],'rxn:T':[3,4]})
        groups=[dict(reaction_id='EX_a',exchange='rxn:EX_a',external='met:a_e',internal='met:a_c')] if triplet else []
        b=dict(fixed_node_offsets={'met:a_c':[4,3]} if triplet else {},routes=[],square=None)
        if barrier:
            xy.update(red_a=[2,0],red_b=[2,8]);nodes.update(red_a=dict(kind='reaction'),red_b=dict(kind='reaction'))
            b['red_first']=dict(routes=[dict(rid='RED',mid='red_b',source='red_a',target='red_b',priority=2,points=[(2,0),(2,8)])])
        return model,nodes,edges,xy,b,groups

    def test_red_first_selects_clear_outer_face_and_preserves_fixed_nodes(self):
        model,nodes,edges,xy,b,groups=self.fixture(barrier=True);before=copy.deepcopy(xy)
        spec=place_corridors(model,nodes,xy,b,{'rxn:EX_a'},groups,9,9)
        self.assertNotEqual(spec['groups'][0]['face'],'left')
        self.assertEqual(xy['met:a_c'],before['met:a_c'])
        for k in ['red_a','red_b']:self.assertEqual(xy[k],before[k])
        points=[xy[k] for k in ['rxn:EX_a','met:a_e','rxn:T','met:a_c']]
        self.assertTrue(len({p[1] for p in points})==1 or len({p[0] for p in points})==1)
        self.assertTrue(all(abs(a[0]-c[0])+abs(a[1]-c[1])>.15 for a,c in zip(points,points[1:])))

    def test_only_real_transport_incidences_are_locked(self):
        model,nodes,edges,xy,b,groups=self.fixture();spec=place_corridors(model,nodes,xy,b,{'rxn:EX_a'},groups,9,9)
        paths=initial_paths(nodes,edges,xy,b,set());before={(p['rid'],p['mid']) for p in paths}
        locked=lock_corridor_paths(paths,spec)
        self.assertEqual(len(locked),3)
        self.assertEqual(before,{(p['rid'],p['mid']) for p in paths})
        self.assertEqual(validate_corridors(paths,spec)['triplets'],1)
        locked[0]['points'].insert(1,[xy['rxn:EX_a'][0]-.1,xy['rxn:EX_a'][1]+.2])
        with self.assertRaisesRegex(ValueError,'not straight'):validate_corridors(paths,spec)

    def test_adaptive_exchange_nodes_use_integer_grid_intersections(self):
        model,nodes,edges,xy,b,groups=self.fixture()
        b['adaptive_grid']=True
        place_corridors(model,nodes,xy,b,{'rxn:EX_a'},groups,9,9)
        for point in xy.values():
            self.assertTrue(all(abs(v-round(v))<1e-7 for v in point))

    def test_periplasm_transport_chain_keeps_real_species_and_incidences(self):
        model,nodes,edges,xy,b,groups=self.fixture()
        model['metabolites'].append(dict(id='a_p',name='a_p',compartment='p'))
        model['reactions'][1]['metabolites']={'a_p':-1,'a_c':1}
        model['reactions'].append(dict(id='OUT',name='OUT',metabolites={'a_e':-1,'a_p':1},lower_bound=-10,upper_bound=10))
        before=copy.deepcopy(model);nodes,edges=build_network(model,[])
        xy={'rxn:EX_a':[0,8],'met:a_e':[1,8],'rxn:OUT':[2,8],'met:a_p':[3,8],'rxn:T':[4,8],'met:a_c':[6,8]}
        b.update(adaptive_grid=True,fixed_node_offsets={'met:a_c':[6,8]})
        specification=place_corridors(model,nodes,xy,b,{'rxn:EX_a'},groups,13,13)
        group=specification['groups'][0]
        self.assertEqual(group['transport_path'],['met:a_e','rxn:OUT','met:a_p','rxn:T','met:a_c'])
        paths=initial_paths(nodes,edges,xy,b,set());locked=lock_corridor_paths(paths,specification)
        self.assertEqual(len(locked),5)
        report=validate_corridors(paths,specification)
        self.assertEqual(report['multistep_groups'],1)
        self.assertEqual(report['straight_incidences'],5)
        self.assertEqual(model,before)
        self.assertTrue(all(abs(v-round(v))<1e-7 for point in xy.values() for v in point))

    def test_transport_chain_cannot_invent_a_missing_intermediate_reaction(self):
        model,nodes,edges,xy,b,groups=self.fixture()
        model['metabolites'].append(dict(id='a_p',name='a_p',compartment='p'))
        model['reactions'][1]['metabolites']={'a_p':-1,'a_c':1}
        nodes,edges=build_network(model,[])
        with self.assertRaisesRegex(ValueError,'No displayed source transport path'):
            discover_corridors(model,nodes,{'rxn:EX_a'},groups)

    def test_missing_cytosolic_species_uses_two_points_without_invention(self):
        model,nodes,edges,xy,b,groups=self.fixture(triplet=False);before=copy.deepcopy(model)
        spec=place_corridors(model,nodes,xy,b,{'rxn:EX_a'},groups,9,9)
        self.assertEqual(spec['groups'][0]['kind'],'pair');self.assertNotIn('met:a_c',nodes)
        self.assertEqual(model,before)
        paths=initial_paths(nodes,edges,xy,b,set());lock_corridor_paths(paths,spec)
        self.assertEqual(validate_corridors(paths,spec)['pairs'],1)

    def test_present_but_undisplayed_internal_is_not_silently_omitted(self):
        model,nodes,edges,xy,b,groups=self.fixture();nodes.pop('met:a_c')
        with self.assertRaisesRegex(ValueError,'absent from display'):discover_corridors(model,nodes,{'rxn:EX_a'},[])

    def test_red_fixed_polyline_is_not_recomputed(self):
        paths=[dict(rid='R',mid='b',source='a',target='b',priority=2,backbone=True,red_fixed=True,points=[(1,1),(1,3),(3,3),(3,1)]),
               dict(rid='O',mid='d',source='c',target='d',priority=1,backbone=False,points=[(2,0),(2,.5)])]
        xy=dict(a=(1,1),b=(3,1),c=(2,0),d=(2,.5))
        # Add intermediate lane coordinates as biological helper endpoints.
        xy.update(helper1=(1,3),helper2=(3,3))
        before=copy.deepcopy(paths[0]['points'])
        result=route(paths,xy,attempts=0,progress=lambda *a,**k:None)
        self.assertTrue(result['verified']);self.assertEqual(paths[0]['points'],before)
