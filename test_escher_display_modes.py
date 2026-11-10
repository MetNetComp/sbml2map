"""Default previews show structure; simulation overlays own the data colors."""
import json
from pathlib import Path
import tempfile
import unittest

from export_escher import create_map,write_preview,write_pathway_preview
from reaction_network import build_network


def builder_options(page):
    start=page.index("escher.libs.d3_select('#map'),{")+len("escher.libs.d3_select('#map'),{")
    prefix=page[start:page.index('enable_editing:true',start)].strip().rstrip(',')
    return json.loads('{'+prefix+'}')


class DisplayModeTests(unittest.TestCase):
    def fixture(self):
        model=dict(id='sampled',metabolites=[dict(id='a')],reactions=[
            dict(id='in',metabolites={'a':1},lower_bound=0,upper_bound=10),
            dict(id='out',metabolites={'a':-1},lower_bound=0,upper_bound=10)])
        nodes,edges=build_network(model,[])
        xy={key:[float(i),float(i)] for i,key in enumerate(nodes)}
        data=create_map(model,nodes,edges,xy,layout_method='flux-median-ranked grid map')
        return data

    def test_default_has_no_reaction_data_or_numeric_flux_overlay(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'escher_preview.html'
            write_preview(path,self.fixture())
            page=path.read_text()
            self.assertNotIn('reaction_data',builder_options(page))
            self.assertNotIn('window.fbaFluxes=',page)
            self.assertIn('No flux data loaded.',page)
            self.assertFalse('<b>sampled.flux-median-ranked' in page.split('<script>',1)[0], 'Default title still advertises a flux overlay')

    def test_structure_colors_have_a_data_update_switch(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'escher_preview.html'
            write_pathway_preview(path,self.fixture(),{'in'},set())
            page=path.read_text()
            self.assertIn('<style id="major-pathway-colors">',page)
            self.assertIn('data-pathway-legend',page)
            self.assertNotIn('reaction_data',builder_options(page))
            self.assertIn("callback_manager.set('update_data.pathway_colors'",page)

    def test_fba_uses_flux_only_blue_scale_without_pathway_css(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'escher_fba.html'
            write_preview(path,self.fixture(),{'in':5.0,'out':5.0},{'objective_value':5.0})
            page=path.read_text();options=builder_options(page)
            self.assertEqual(options['reaction_data'],{'in':5.0,'out':5.0})
            colors=[x['color'] for x in options['reaction_scale']]
            self.assertEqual(colors,['#b8b8b8','#deebf7','#9ecae1','#4292c6','#08519c'])
            self.assertNotIn('<style id="major-pathway-colors">',page)
            self.assertIn('Fixed pathway colors are not applied.',page)

    def test_saved_e_coli_views_have_distinct_initial_modes(self):
        root=Path(__file__).resolve().parent;folder=root/'results/e_coli_core'
        if not folder.exists():folder=root/'public/results/e_coli_core'
        normal=(folder/'escher_preview.html').read_text()
        fba=(folder/'escher_fba.html').read_text()
        self.assertNotIn('reaction_data',builder_options(normal))
        self.assertNotIn('window.fbaFluxes=',normal)
        self.assertIn('reaction_data',builder_options(fba))
        self.assertNotIn('<style id="major-pathway-colors">',fba)


if __name__=='__main__':unittest.main()
