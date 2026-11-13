import copy
import unittest

from sequential_orange import organic_graph, hop_distances, candidate_paths, place_exchange_branches, branch_milp


class SequentialOrangeTests(unittest.TestCase):
    def fixture(self):
        nodes={k:dict(kind='reaction' if k.startswith('rxn:') else 'metabolite')
               for k in ('rxn:R','met:a','rxn:Q','met:m','rxn:EX_A','rxn:EX_B')}
        edges=[]
        for a,b in [('rxn:R','met:a'),('rxn:Q','met:a'),('rxn:Q','met:m'),('rxn:EX_A','met:m'),('rxn:EX_B','met:m')]:
            edges.append(dict(source=a,target=b,reaction_id=a[4:],metabolite_id=b[4:]))
        red=[dict(source='rxn:R',target='met:a',rid='R',mid='a',priority=2,points=[[4.,4.],[6.,4.]],red_fixed=True)]
        return nodes,edges,red

    def test_shared_cofactor_is_excluded_even_if_displayed_as_ordinary(self):
        nodes,edges,red=self.fixture()
        nodes['met:m']['shared_cofactor']=True
        adjacency,_=organic_graph(nodes,edges)
        self.assertNotIn('met:m',adjacency)
        self.assertNotIn('rxn:EX_A',hop_distances(adjacency,{'met:a'}))

    def test_order_uses_initial_red_distance_not_current_orange_distance(self):
        nodes,edges,_=self.fixture();adjacency,_=organic_graph(nodes,edges)
        initial=hop_distances(adjacency,{'met:a','rxn:R'})
        self.assertEqual(initial['rxn:EX_A'],3)
        self.assertEqual(initial['rxn:EX_B'],3)
        chains=candidate_paths(adjacency,'rxn:EX_B',{'met:a','rxn:R','met:m'}, {'rxn:EX_A','rxn:EX_B'})
        self.assertEqual(chains,[['rxn:EX_B','met:m']])

    def test_milp_reuses_orange_without_charging_existing_length_and_keeps_red(self):
        nodes,edges,red=self.fixture();before=copy.deepcopy(red)
        xy={'rxn:R':[4.,4.],'met:a':[6.,4.]}
        result=place_exchange_branches(nodes,edges,xy,red,{'rxn:EX_A','rxn:EX_B'},[0,10,0,10],None,lambda *a,**k:None,5.)
        first,second=result['records']
        self.assertEqual([r['red_hops'] for r in result['records']],[3,3])
        self.assertEqual(first['attachment_color'],'red')
        self.assertEqual(second['attachment_color'],'orange')
        self.assertEqual(second['attachment_node'],'met:m')
        self.assertTrue(first['solver']['optimality_proven'])
        self.assertTrue(second['solver']['optimality_proven'])
        self.assertEqual(result['routes'][:1],before)
        self.assertFalse(result['disconnected_exchange_node_ids'])
        self.assertLess(second['new_length_grid'],first['new_length_grid']+second['new_length_grid'])

    def test_fifth_branch_shares_only_the_short_common_node_port(self):
        from adaptive_routing import validate_paths
        nodes,edges,_=self.fixture()
        xy={'met:a':[5.,5.]};frozen=[]
        for i,point in enumerate(([3.,5.],[5.,3.],[5.,7.],[7.,5.])):
            key='rxn:fixed'+str(i);nodes[key]=dict(kind='reaction');xy[key]=point
            frozen.append(dict(rid=key[4:],mid='a',source=key,target='met:a',priority=2 if i<3 else 1,points=[point,[5.,5.]]))
        _,lookup=organic_graph(nodes,edges)
        branch,solver=branch_milp(['rxn:EX_A','met:m','rxn:Q','met:a'],nodes,edges,lookup,xy,frozen,[0,10,0,10],None,10.)
        self.assertIsNotNone(branch,solver)
        self.assertTrue(any(p.get('auxiliary_lane') for p in branch['routes']))
        self.assertTrue(validate_paths(frozen+branch['routes'],dict(xy,**branch['positions']),[0,10,0,10],None,[])['verified'])


if __name__=='__main__':unittest.main()
