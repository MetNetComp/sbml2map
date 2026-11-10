"""Label visibility options preserve chemistry and structural priorities."""
import copy
import tempfile
from pathlib import Path
import unittest

from export_escher import configure_label_display, write_pathway_preview
import test_escher_display_modes
from sbml2map import parser


class LabelDisplayTests(unittest.TestCase):
    def test_cli_defaults_and_overrides(self):
        defaults = parser().parse_args([])
        self.assertFalse(defaults.hide_metabolite_labels)
        self.assertFalse(defaults.hide_reaction_labels)
        self.assertFalse(defaults.allow_label_overlap)
        args = parser().parse_args(['--hide-metabolite-labels', '--hide-reaction-labels', '--allow-label-overlap'])
        self.assertTrue(args.hide_metabolite_labels and args.hide_reaction_labels and args.allow_label_overlap)

    def test_priorities_and_local_cofactor_copies(self):
        nodes = {str(i): dict(node_type='metabolite', node_is_primary=i < 4) for i in range(6)}
        reactions = {}
        for i, identity in enumerate(['BOUNDARY', 'red', 'orange', 'gray']):
            reactions[str(i)] = dict(bigg_id=identity, segments={str(i): dict(from_node_id=str(i), to_node_id='4' if i == 1 else '5')})
        data = [dict(metadata=dict(exchange_reaction_ids=['BOUNDARY'], segment_priorities={'1': 2, '2': 1})),
                dict(nodes=nodes, reactions=reactions)]
        original = copy.deepcopy(data[1])
        options = configure_label_display(data)
        self.assertEqual(options['reaction_priorities'], {'0': 5, '1': 4, '2': 3, '3': 2})
        self.assertEqual(options['metabolite_priorities'], {'0': 5, '1': 4, '2': 3, '3': 2, '4': 1, '5': 5})
        self.assertEqual(data[1], original)

    def test_explicit_major_and_ex_prefix(self):
        data = [dict(metadata={}), dict(nodes={}, reactions={
            'a': dict(bigg_id='EX_glc', segments={}), 'b': dict(bigg_id='main', segments={})})]
        self.assertEqual(configure_label_display(data, ['main'])['reaction_priorities'], {'a': 5, 'b': 4})

    def test_preview_embeds_options_controls_and_keeps_labels(self):
        data = test_escher_display_modes.DisplayModeTests().fixture()
        data[0]['metadata']['label_display'] = dict(show_metabolites=False, show_reactions=False, allow_overlap=True)
        original = copy.deepcopy(data[1])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'preview.html'
            write_pathway_preview(path, data, {'in'})
            page = path.read_text()
        self.assertIn('id="allow-label-overlap"', page)
        self.assertIn('"show_metabolites": false', page)
        self.assertIn('"show_reactions": false', page)
        self.assertIn('"allow_overlap": true', page)
        self.assertEqual(data[1], original)


if __name__ == '__main__':
    unittest.main()
