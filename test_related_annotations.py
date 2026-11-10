"""Related chemical annotations support discovery without replacing identities."""
import copy
import json
from pathlib import Path
import unittest

from auto_pathways import annotation_ids
from bacterial_pathways import REFERENCE, discover, map_roles
from model_input import load_model


def annotation(value, qualifier='isRelatedTo', database='bigg.metabolite'):
    return dict(qualifier='{http://biomodels.net/biology-qualifiers/}' + qualifier,
                uri='https://identifiers.org/' + database + ':' + value)


class RelatedAnnotationTests(unittest.TestCase):
    def mapping(self, annotations, overrides=None):
        model=dict(metabolites=[dict(id='opaque', name='', annotations=annotations)])
        return map_roles(model, json.loads(REFERENCE.read_text()), overrides)

    def test_related_identifiers_are_used_and_audited(self):
        mapped, records=self.mapping([annotation('glc__D')])
        self.assertEqual(mapped['opaque'], 'glc__D')
        self.assertEqual(records[0]['evidence'][0]['source'], 'annotation:bigg.metabolite:isRelatedTo')
        # The shared reader still requires explicit opt-in for other consumers.
        self.assertEqual(annotation_ids(dict(annotations=[annotation('glc__D')]), 'bigg.metabolite'), [])

    def test_identity_annotations_take_precedence_across_databases(self):
        for values in ([annotation('pyr', 'is'), annotation('accoa')],
                       [annotation('C00022', 'is', 'kegg.compound'), annotation('accoa')]):
            with self.subTest(values=values):
                mapped, records=self.mapping(values)
                self.assertEqual(mapped['opaque'], 'pyr')
                self.assertIn('accoa', records[0]['lower_tier_conflicts'])

    def test_conflicting_related_roles_stay_ambiguous_and_allow_override(self):
        values=[annotation('pyr'), annotation('accoa')]
        mapped, records=self.mapping(values)
        self.assertNotIn('opaque', mapped)
        self.assertEqual(records[0]['status'], 'ambiguous')
        self.assertEqual(self.mapping(values, {'opaque': 'pyr'})[0]['opaque'], 'pyr')

    def test_other_qualifiers_are_not_used(self):
        mapped, records=self.mapping([annotation('pyr', 'isPartOf')])
        self.assertEqual(mapped, {})
        self.assertEqual(records[0]['status'], 'unmapped')

    def test_pathways_with_only_related_annotations_and_opaque_ids(self):
        base=Path(__file__).resolve().parent
        path=base/'models/e_coli_core.xml' if (base/'models').exists() else base/'e_coli_core.xml'
        model=load_model(path)
        mids={m['id']: 'opaque_' + str(i) for i, m in enumerate(model['metabolites'])}
        for m in model['metabolites']:
            role=m['id'].removesuffix('_' + m['compartment'])
            m.update(id=mids[m['id']], name='', annotation={}, annotations=[annotation(role)])
        for r in model['reactions']:
            r['metabolites']={mids[mid]: value for mid, value in r['metabolites'].items()}
        original=copy.deepcopy(model)
        routes, report, _, _=discover(model)
        self.assertEqual(report['status'], 'complete_central_backbone')
        self.assertTrue(next(r for r in routes if r['id']=='tca')['closed'])
        self.assertEqual(model, original)


if __name__=='__main__':
    unittest.main()
