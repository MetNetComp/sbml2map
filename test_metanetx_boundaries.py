"""Chemical identity and exchange discovery across MNXref representations."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from auto_pathways import annotation_ids
from bacterial_pathways import REFERENCE, discover, exchange_reactions, map_roles
from chemical_identity import MetaNetXReference, metanetx_id
from annotation_context import reaction_labels
from escher_grid import classify_metabolites, display_graph
from model_input import load_model
from sbml2map import discovery_context, BacterialStrategy


class MetaNetXIdentityTests(unittest.TestCase):
    def setUp(self):self.reference=json.loads(REFERENCE.read_text())

    def test_uri_namespace_and_chebi_variants(self):
        variants=('https://identifiers.org/CHEBI:57579',
                  'https://identifiers.org/chebi/CHEBI:57579',
                  'urn:miriam:chebi:57579')
        for uri in variants:
            entity=dict(annotations=[dict(qualifier='{q}is',uri=uri)])
            self.assertEqual(annotation_ids(entity,'chebi'),['CHEBI:57579'])
        self.assertEqual(annotation_ids(dict(annotation={'CHEBI':['57579']}),'chebi'),['CHEBI:57579'])

    def test_mnx_id_with_compartment_and_source_notes(self):
        mid=self.reference['roles']['pyr']['aliases']['metanetx.chemical'][0]
        model=dict(metabolites=[dict(id=mid+'@MNXC3',name='')])
        mapped,records=map_roles(model,self.reference)
        self.assertEqual(mapped[mid+'@MNXC3'],'pyr')
        self.assertIsNone(metanetx_id('MNXM23_unrecognized_suffix'))

    def tables(self,base,xrefs='bigg.metabolite:f6p\tMNXM900001\tF6P\n',depr='MNXM900000\tMNXM900001\t4.5\n'):
        (base/'chem_xref.tsv').write_text('# MNXref fixture release 4.5\n'+xrefs)
        (base/'chem_depr.tsv').write_text('# release 4.5\n'+depr)
        return MetaNetXReference(self.reference,base)

    def test_local_xrefs_old_ids_and_release_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory);resolver=self.tables(base)
            model=dict(id='opaque',metabolites=[dict(id='MNXM900000@cell',name='')],reactions=[])
            self.assertEqual(map_roles(model,self.reference,mnxref=resolver)[0]['MNXM900000@cell'],'f6p')
            _,report,_,_=discover(model,mnxref_dir=base)
            files=report['metanetx_reference']['files']
            self.assertEqual(len(files['chem_xref.tsv']['sha256']),64)
            self.assertIn('4.5',files['chem_xref.tsv']['headers'][0])

    def test_split_ids_stay_ambiguous_and_manual_override_wins(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory)
            resolver=self.tables(base,'bigg.metabolite:f6p\tMNXM900001\tF6P\nbigg.metabolite:pyr\tMNXM900002\tPyruvate\n',
                                 'MNXM900000\tMNXM900001\t4.5\nMNXM900000\tMNXM900002\t4.5\n')
            model=dict(metabolites=[dict(id='MNXM900000@cell',name='')])
            mapped,records=map_roles(model,self.reference,mnxref=resolver)
            self.assertFalse(mapped);self.assertEqual(records[0]['status'],'ambiguous')
            self.assertEqual(set(records[0]['candidate_roles']),{'f6p','pyr'})
            self.assertEqual(map_roles(model,self.reference,{'MNXM900000@cell':'f6p'},resolver)[0]['MNXM900000@cell'],'f6p')

    def test_bad_local_tables_are_explicit_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory)
            with self.assertRaisesRegex(ValueError,'Missing MetaNetX'):MetaNetXReference(self.reference,base)
            with self.assertRaisesRegex(ValueError,'Cycle'):
                self.tables(base,depr='MNXM900000\tMNXM900001\t4.5\nMNXM900001\tMNXM900000\t4.5\n')


def boundary_xml(sbo='SBO:0000627',other='MNXM900001',note='SOURCE: EX_glc__D_e'):
    return f'''<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"><model id="boundary">
    <listOfCompartments><compartment id="MNXC2" name="extracellular region"/>
    <compartment id="BOUNDARY" name="model boundary"/></listOfCompartments>
    <listOfSpecies><species id="MNXM900001__64__MNXC2" name="unknown" compartment="MNXC2">
    <notes>SOURCE: M_glc__D_e FORMULA: C6H12O6</notes></species>
    <species id="{other}__64__BOUNDARY" compartment="BOUNDARY"/></listOfSpecies>
    <listOfReactions><reaction id="opaque_exchange" sboTerm="{sbo}" reversible="true"><notes>{note}</notes>
    <listOfReactants><speciesReference species="MNXM900001__64__MNXC2"/></listOfReactants>
    <listOfProducts><speciesReference species="{other}__64__BOUNDARY"/></listOfProducts>
    </reaction></listOfReactions></model></sbml>'''


class BoundaryRecognitionTests(unittest.TestCase):
    def context(self,xml):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'fixture.xml';p.write_text(xml);source=load_model(p)
            return source,discovery_context(p,source)

    def test_boundary_sbo_and_source_hide_only_the_display_participant(self):
        source,context=self.context(boundary_xml())
        original=copy.deepcopy(source);strategy=BacterialStrategy(context)
        self.assertEqual(strategy.mapped['MNXM900001@MNXC2'],'glc__D')
        exchange=strategy.extraction['exchanges'][0]
        self.assertEqual(exchange['reaction_id'],'opaque_exchange')
        self.assertEqual(exchange['boundary_placeholder_id'],'MNXM900001@BOUNDARY')
        self.assertTrue(exchange['import_allowed'])
        classes=strategy.classify(source)
        projected,_,_=display_graph(source,{'opaque_exchange'},classes)
        self.assertEqual(projected['reactions'][0]['metabolites'],{'MNXM900001@MNXC2':-1})
        self.assertEqual(source,original)

    def test_sbo_alone_or_related_ex_annotation_is_sufficient(self):
        for use_sbo in (True,False):
            xml=boundary_xml(sbo='SBO:0000627' if use_sbo else '',note='')
            if not use_sbo:
                xml=xml.replace('<notes></notes>','''<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
                xmlns:bqbiol="http://biomodels.net/biology-qualifiers/"><rdf:Description rdf:about="#meta">
                <bqbiol:isRelatedTo><rdf:Bag><rdf:li rdf:resource="https://identifiers.org/bigg.reaction:EX_glc__D_e"/>
                </rdf:Bag></bqbiol:isRelatedTo></rdf:Description></rdf:RDF></annotation>''')
                xml=xml.replace('id="opaque_exchange"','id="opaque_exchange" metaid="meta"')
            _,context=self.context(xml);exchanges=BacterialStrategy(context).extraction['exchanges']
            self.assertEqual(len(exchanges),1)

    def test_different_boundary_chemical_and_unlabeled_transport_are_rejected(self):
        for xml in (boundary_xml(other='MNXM900002'),boundary_xml(sbo='',note='')):
            _,context=self.context(xml)
            self.assertFalse(any(m.get('_boundary_placeholder') for m in context['metabolites']))
            self.assertFalse(exchange_reactions(context,{}))

    def test_source_participant_count_and_demand_sink_subtypes(self):
        model=dict(metabolites=[dict(id='x',compartment='c'),dict(id='h',compartment='c')],
                   reactions=[dict(id=rid,metabolites={'x':-1},lower_bound=0,upper_bound=10)
                              for rid in ('unlabeled','DM_x','SK_x')])
        model['reactions'].append(dict(id='normal',metabolites={'x':-1,'h':1},lower_bound=0,upper_bound=10))
        kinds={r['reaction_id']:r['boundary_kind'] for r in exchange_reactions(model,{})}
        self.assertEqual(kinds,{'unlabeled':'exchange_candidate','DM_x':'demand','SK_x':'sink'})
        self.assertNotIn('normal',kinds)

    def test_related_sink_alias_does_not_replace_explicit_exchange(self):
        reaction=dict(id='opaque',notes='SOURCE: EX_x_e REACTION: MNXR1',
                      metabolites={'x':-1},lower_bound=-10,upper_bound=1000,
                      annotations=[dict(qualifier='{q}isRelatedTo',uri='https://identifiers.org/bigg.reaction:DM_x')])
        model=dict(metabolites=[dict(id='x',compartment='e')],reactions=[reaction])
        self.assertEqual(exchange_reactions(model,{})[0]['boundary_kind'],'exchange')

    @unittest.skipUnless((Path(__file__).parent/'data/misc_iCH360.sbml3.xml').exists(),'Local model is not distributed')
    def test_local_misc_model_has_all_central_pathways_and_original_exchanges(self):
        path=Path(__file__).parent/'data/misc_iCH360.sbml3.xml';source=load_model(path);original=copy.deepcopy(source)
        strategy=BacterialStrategy(discovery_context(path,source))
        self.assertEqual(strategy.extraction['status'],'complete_central_backbone')
        routes={r['id']:r for r in strategy.routes}
        self.assertEqual(len(routes['glycolysis']['reaction_ids']),11)
        self.assertTrue(routes['tca']['closed'])
        exchanges=strategy.extraction['exchanges']
        self.assertEqual(sum(any(v.startswith('EX_') for _,v in reaction_labels(
            next(r for r in strategy.context['reactions'] if r['id']==e['reaction_id']))) for e in exchanges),26)
        self.assertTrue(all(e['boundary_kind']=='exchange' for e in exchanges))
        self.assertEqual(source,original)


if __name__=='__main__':unittest.main()
