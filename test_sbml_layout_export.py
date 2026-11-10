"""Opt-in Layout export preserves the original full SBML and saved geometry."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from export_sbml_layout import export_sbml_layout,validate_export_request
from model_input import load_model
from sbml2map import parser,run

try:
    import libsbml
except ImportError:
    libsbml=None

BASE=Path(__file__).resolve().parent
MODEL=BASE/'models/e_coli_core.xml' if (BASE/'models').exists() else BASE/'e_coli_core.xml'
RESULTS=BASE/'results/e_coli_core' if (BASE/'results/e_coli_core').exists() else BASE/'public/results/e_coli_core'


class DefaultExportTests(unittest.TestCase):
    def test_default_no_export_even_without_optional_dependency(self):
        self.assertFalse(parser().parse_args([]).export_sbml_layout)
        with tempfile.TemporaryDirectory() as temporary:
            args=parser().parse_args(['--model',str(MODEL),'--no-sampling','--coverage','0',
                                      '--output-dir',temporary,'--time-limit','1','--max-milp-nodes','1'])
            with patch.dict('sys.modules',{'libsbml':None}),contextlib.redirect_stdout(io.StringIO()):run(args)
            self.assertFalse((Path(temporary)/'layout.sbml').exists())
            self.assertFalse((Path(temporary)/'sbml_layout_report.json').exists())

    def test_json_level2_and_invalid_layout_id_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'not JSON'):validate_export_request(RESULTS/'source_model.json')
        with self.assertRaisesRegex(ValueError,'valid SBML SId'):validate_export_request(MODEL,'bad id')
        with tempfile.TemporaryDirectory() as temporary:
            p=Path(temporary)/'level2.xml';p.write_text('<sbml level="2" version="4"><model id="m"/></sbml>')
            with self.assertRaisesRegex(ValueError,'Level 3'):validate_export_request(p)

    def test_optional_dependency_message_and_discovery_conflict(self):
        with patch.dict('sys.modules',{'libsbml':None}),self.assertRaisesRegex(RuntimeError,'requirements-sbml-layout'):
            validate_export_request(MODEL)
        args=parser().parse_args(['--model',str(MODEL),'--discover-only','--export-sbml-layout'])
        with self.assertRaisesRegex(ValueError,'discover-only'):run(args)


@unittest.skipUnless(libsbml is not None,'Optional python-libsbml is not installed')
class LayoutExportTests(unittest.TestCase):
    def test_full_model_coordinates_copies_and_fba_round_trip(self):
        original=MODEL.read_bytes();model=load_model(MODEL)
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'layout.sbml'
            report=export_sbml_layout(MODEL,RESULTS,path,verify_fba=True)
            self.assertEqual(report['status'],'valid');self.assertFalse(report['render_exported'])
            self.assertEqual(report['source_reactions'],95);self.assertEqual(report['reaction_glyphs'],76)
            self.assertEqual(report['species_glyphs'],130)
            self.assertEqual(report['line_segments'],sum(len(r['segments']) for r in json.loads((RESULTS/'escher_map.json').read_text())[1]['reactions'].values()))
            self.assertTrue(report['undisplayed_reaction_ids'])
            self.assertEqual(load_model(path),model)
            self.assertAlmostEqual(report['fba_verification']['after'],0.8739215069684304,places=8)
            document=libsbml.readSBMLFromFile(str(path));layout=document.getModel().getPlugin('layout').getLayout(report['layout_id'])
            coa=[g for g in layout.getListOfSpeciesGlyphs() if g.getSpeciesId()=='M_coa_c']
            self.assertGreater(len(coa),1)
            for reaction in layout.getListOfReactionGlyphs():
                for glyph in reaction.getListOfSpeciesReferenceGlyphs():
                    self.assertIsNotNone(layout.getSpeciesGlyph(glyph.getSpeciesGlyphId()))
                    species=layout.getSpeciesGlyph(glyph.getSpeciesGlyphId()).getSpeciesId()
                    original_reaction=document.getModel().getReaction(reaction.getReactionId())
                    is_substrate=any(ref.getSpecies()==species for ref in original_reaction.getListOfReactants())
                    self.assertEqual(glyph.getRole(),libsbml.SPECIES_ROLE_SUBSTRATE if is_substrate else libsbml.SPECIES_ROLE_PRODUCT)
                    for segment in glyph.getCurve().getListOfCurveSegments():
                        a,b=segment.getStart(),segment.getEnd()
                        self.assertTrue(abs(a.getXOffset()-b.getXOffset())<1e-6 or abs(a.getYOffset()-b.getYOffset())<1e-6)
        self.assertEqual(MODEL.read_bytes(),original)

    def test_existing_layout_and_annotations_are_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            first=Path(temporary)/'first.sbml';second=Path(temporary)/'second.sbml'
            first_report=export_sbml_layout(MODEL,RESULTS,first)
            first.write_text(first.read_text().replace('</model>','<!-- keep-export-test-comment --></model>'))
            second_report=export_sbml_layout(first,RESULTS,second)
            self.assertEqual(second_report['existing_layouts_preserved'],1)
            self.assertNotEqual(first_report['layout_id'],second_report['layout_id'])
            self.assertIn('keep-export-test-comment',second.read_text())
            doc=libsbml.readSBMLFromFile(str(second))
            self.assertEqual(doc.getModel().getPlugin('layout').getNumLayouts(),2)
            self.assertTrue(second_report['validation']['original_xml_semantics_preserved'])

    def test_original_escaped_reaction_ids_are_referenced(self):
        with tempfile.TemporaryDirectory() as temporary:
            source=Path(temporary)/'escaped.xml';output=Path(temporary)/'layout.sbml'
            source.write_text(MODEL.read_text().replace(' id="R_PFK"',' id="R_P__70__K"'))
            self.assertEqual(load_model(source),load_model(MODEL))
            report=export_sbml_layout(source,RESULTS,output)
            doc=libsbml.readSBMLFromFile(str(output));layout=doc.getModel().getPlugin('layout').getLayout(report['layout_id'])
            self.assertIn('R_P__70__K',[g.getReactionId() for g in layout.getListOfReactionGlyphs()])

    def test_mismatched_source_and_input_overwrite_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'overwrite'):export_sbml_layout(MODEL,RESULTS,MODEL)
        with tempfile.TemporaryDirectory() as temporary:
            source=Path(temporary)/'changed.xml';output=Path(temporary)/'layout.sbml'
            document=libsbml.readSBMLFromFile(str(MODEL));document.getModel().getParameter('R_EX_glc__D_e_lower_bound').setValue(-5)
            libsbml.writeSBMLToFile(document,str(source))
            with self.assertRaisesRegex(ValueError,'does not match'):export_sbml_layout(source,RESULTS,output)
            self.assertFalse(output.exists())

    def test_explicit_cli_flag_exports_empty_layout_without_replacing_model(self):
        with tempfile.TemporaryDirectory() as temporary:
            args=parser().parse_args(['--model',str(MODEL),'--no-sampling','--coverage','0',
                                      '--output-dir',temporary,'--time-limit','1','--max-milp-nodes','1','--export-sbml-layout'])
            with contextlib.redirect_stdout(io.StringIO()):report=run(args)
            self.assertTrue((Path(temporary)/'layout.sbml').exists())
            self.assertEqual(load_model(Path(temporary)/'layout.sbml'),load_model(MODEL))
            self.assertEqual(report['sbml_layout']['file'],'layout.sbml')


if __name__=='__main__':unittest.main()
