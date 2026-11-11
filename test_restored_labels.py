"""Label controls and halos apply to all published viewer types."""
import json
from pathlib import Path
import re
import tempfile
import unittest
import xml.etree.ElementTree as ET
from escher_grid import render
from reaction_network import build_network
from test_escher_display_modes import builder_options

ROOT=Path(__file__).resolve().parent


class RestoredLabelTests(unittest.TestCase):
    def test_all_published_viewers_have_three_label_controls(self):
        for name in ('e_coli_core/index.html','e_coli_core/escher_preview.html','e_coli_core/escher_fba.html',
                     'eco_ich_360/index.html','eco_ich_360/escher_preview.html'):
            with self.subTest(view=name):
                page=(ROOT/'docs'/name).read_text()
                for control in ('show-metabolite-labels','show-reaction-labels','allow-label-overlap'):
                    self.assertIn(f'id="{control}"',page)

    def test_svg_render_honors_initial_flags_and_has_white_halos(self):
        model=dict(id='labels',metabolites=[dict(id='a'),dict(id='b')],reactions=[
            dict(id='r',metabolites={'a':-1,'b':1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={'met:a':[0,0],'rxn:r':[1,0],'met:b':[2,0]}
        records=[dict(reaction_id='r',selected=True,is_major_pathway=True,rank=1)]
        selection=dict(mode='structure',coverage_percent=100,eligible_reactions=1)
        backbones=dict(routes=[dict(metabolite_ids=['a','b'],reaction_ids=['r'],steps=[dict(reaction_id='r',source='a',target='b')])])
        options=dict(show_metabolites=False,show_reactions=False,allow_overlap=True)
        with tempfile.TemporaryDirectory() as temporary:
            folder=Path(temporary)
            render(folder,nodes,edges,xy,records,selection,model['id'],backbones,source_model=model,label_options=options)
            page=(folder/'index.html').read_text()
            actual=json.loads(re.search(r'window.gridLabelOptions=(.*?);</script>',page).group(1))
            self.assertEqual(actual,options)
            self.assertIn('label_kind',page)
            svg=ET.parse(folder/'layout.svg')
            labels=[e for e in svg.iter() if (e.get('id') or '').startswith('label-')]
            self.assertEqual(len(labels),3)
            for group in labels:
                self.assertTrue(any('stroke: #ffffff' in e.get('style','') for e in group.iter()),group.get('id'))

    def test_fba_display_controls_preserve_signed_fluxes_and_scale(self):
        folder=ROOT/'results/e_coli_core';page=(folder/'escher_fba.html').read_text()
        options=builder_options(page)
        self.assertEqual(options['reaction_data'],json.loads((folder/'escher_fluxes.json').read_text()))
        self.assertEqual(options['reaction_scale'][0]['color'],'#b8b8b8')
        self.assertEqual(options['reaction_scale'][1]['color'],'#2166ac')
        self.assertEqual(options['reaction_scale'][-1]['color'],'#d73027')
        self.assertEqual(options['reaction_scale'][-1]['size'],8)
        self.assertIn('square-root scale',page)
        data=json.loads(re.search(r'const mapData=(.*);\n',page).group(1))
        for key in ('show_metabolites','show_reactions','allow_overlap'):
            self.assertIn(key,data[0]['metadata']['label_display'])
        report=json.loads((folder/'fba_report.json').read_text())
        self.assertEqual(data[0]['metadata']['fba']['objective_value'],report['objective_value'])

