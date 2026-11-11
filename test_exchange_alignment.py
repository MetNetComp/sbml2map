"""EX triplets require chemical/compartment evidence and ordered straight positions."""
import unittest
from exchange_alignment import discover_alignment,validate_alignment


class ExchangeAlignmentTests(unittest.TestCase):
    def fixture(self):
        model=dict(metabolites=[dict(id='a_e',compartment='e'),dict(id='a_c',compartment='c')],
                   reactions=[dict(id='EX_a',metabolites={'a_e':-1})])
        nodes={'rxn:EX_a':dict(original_id='EX_a'),'met:a_e':{},'met:a_c':{}}
        return model,nodes

    def test_matching_source_ids_and_source_chemistry_are_preserved(self):
        import copy
        model,nodes=self.fixture();before=copy.deepcopy(model)
        groups,skipped=discover_alignment(model,nodes,{'rxn:EX_a'})
        self.assertEqual(len(groups),1);self.assertEqual(skipped,[])
        self.assertEqual(groups[0]['internal'],'met:a_c')
        self.assertEqual(model,before)

    def test_exact_annotations_match_different_ids(self):
        model,nodes=self.fixture();model['metabolites'][1]['id']='different_c';nodes['met:different_c']=nodes.pop('met:a_c')
        for m in model['metabolites']:m['annotations']=[dict(qualifier='{http://biomodels.net/biology-qualifiers/}is',uri='https://identifiers.org/chebi/CHEBI:123')]
        groups,_=discover_alignment(model,nodes,{'rxn:EX_a'})
        self.assertEqual(groups[0]['internal'],'met:different_c')

    def test_ambiguous_or_absent_matches_are_reported(self):
        model,nodes=self.fixture();nodes.pop('met:a_c')
        groups,skipped=discover_alignment(model,nodes,{'rxn:EX_a'})
        self.assertEqual(groups,[]);self.assertEqual(len(skipped),1)
        model,nodes=self.fixture()
        for m in model['metabolites']:m['annotations']=[dict(qualifier='is',uri='https://identifiers.org/chebi/CHEBI:123')]
        model['metabolites'].append(dict(model['metabolites'][1],id='other_c'));nodes['met:other_c']={}
        groups,skipped=discover_alignment(model,nodes,{'rxn:EX_a'})
        self.assertEqual(groups,[]);self.assertEqual(len(skipped),1)

    def test_all_four_faces_have_ordered_collinear_triplets(self):
        group=dict(reaction_id='EX_a',exchange='ex',external='e',internal='c')
        cases={'left':[(0,2),(1,2),(3,2)],'right':[(4,2),(3,2),(1,2)],'bottom':[(2,0),(2,1),(2,3)],'top':[(2,4),(2,3),(2,1)]}
        for face,positions in cases.items():
            report=validate_alignment([group],dict(zip(('ex','e','c'),positions)),5,5)
            self.assertEqual(report['aligned_triplets'][0]['face'],face)
        with self.assertRaises(RuntimeError):validate_alignment([group],dict(ex=[0,2],e=[2,2],c=[1,2]),5,5)
