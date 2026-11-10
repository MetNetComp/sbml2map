"""Sampled percentages retain pathway priority, chemistry, and saved coordinates."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from export_escher import configure_reaction_coverage, write_pathway_preview
import test_escher_display_modes


class EscherCoverageTests(unittest.TestCase):
    def test_pathway_priority_then_absolute_sampled_median(self):
        data=[dict(metadata={}),dict(reactions={
            '1':dict(bigg_id='low'), '2':dict(bigg_id='high'),
            '3':dict(bigg_id='main'), '4':dict(bigg_id='tie')})]
        records=[dict(reaction_id='low',median=1),dict(reaction_id='high',median=-9),
                 dict(reaction_id='main',median=.1,is_major_pathway=True),dict(reaction_id='tie',median=9)]
        original=copy.deepcopy(data[1]);configure_reaction_coverage(data,records)
        self.assertEqual(data[0]['metadata']['reaction_coverage']['reaction_order'],['3','2','4','1'])
        self.assertEqual(data[1],original)

    def test_red_then_orange_override_larger_gray_fluxes(self):
        data=[dict(metadata=dict(segment_priorities={'s1':2,'s2':1,'s3':1,'s4':0})),
              dict(reactions={
                  '1':dict(bigg_id='red',segments={'s1':{}}),
                  '2':dict(bigg_id='orange_low',segments={'s2':{}}),
                  '3':dict(bigg_id='orange_high',segments={'s3':{}}),
                  '4':dict(bigg_id='gray',segments={'s4':{}})})]
        records=[dict(reaction_id='red',median=.01),dict(reaction_id='orange_low',median=.1),
                 dict(reaction_id='orange_high',median=-1),dict(reaction_id='gray',median=1000)]
        configure_reaction_coverage(data,records)
        self.assertEqual(data[0]['metadata']['reaction_coverage']['reaction_order'],['1','3','2','4'])

    def test_structural_or_incomplete_statistics_do_not_enable_sampled_filter(self):
        for records in ([],[dict(reaction_id='in',rank=1)], [dict(reaction_id='in',median=1)]):
            data=test_escher_display_modes.DisplayModeTests().fixture()
            configure_reaction_coverage(data,records)
            self.assertNotIn('reaction_coverage',data[0]['metadata'])

    def test_preview_embeds_sampled_order_without_removing_reactions(self):
        data=test_escher_display_modes.DisplayModeTests().fixture()
        original=copy.deepcopy(data[1])
        records=[dict(reaction_id='in',median=-1,is_major_pathway=True),dict(reaction_id='out',median=8)]
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'preview.html'
            write_pathway_preview(path,data,{'in'},records=records)
            page=path.read_text()
        self.assertEqual(data[1],original)
        embedded=json.loads(page.split('const mapData=',1)[1].split(';\nwindow.builder=',1)[0])
        ordered=embedded[0]['metadata']['reaction_coverage']['reaction_order']
        self.assertEqual([data[1]['reactions'][rid]['bigg_id'] for rid in ordered],['in','out'])
        self.assertIn('id="escher-coverage"',page)
        self.assertIn('data-escher-coverage-hidden',page)
        self.assertIn('window.syncEscherCoverage=sync',page)

    def test_nonfinite_sampled_median_is_rejected(self):
        data=[dict(metadata={}),dict(reactions={'1':dict(bigg_id='bad')})]
        with self.assertRaisesRegex(ValueError,'Nonfinite sampled median'):
            configure_reaction_coverage(data,[dict(reaction_id='bad',median=float('nan'))])


if __name__=='__main__':unittest.main()
