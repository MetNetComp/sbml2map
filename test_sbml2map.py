"""The renamed entry point keeps existing imports and saved maps usable."""
import copy
import json
from pathlib import Path
import unittest

import sbml2map
import sbml2escher
from escher_routing import separate_routes


class Sbml2mapTests(unittest.TestCase):
    def test_new_command_defaults(self):
        parser=sbml2map.parser();args=parser.parse_args([])
        self.assertEqual(parser.prog,'sbml2map.py')
        self.assertEqual(args.program_name,'sbml2map')
        self.assertEqual(args.output_dir.name,'e_coli_core_sbml2map_output')
        self.assertEqual(args.sbml_layout_id,'sbml2map_layout')
        self.assertEqual(args.orange_flux_weight,1)

    def test_legacy_entry_point_uses_the_same_implementation(self):
        self.assertIs(sbml2escher.run,sbml2map.run)
        self.assertIs(sbml2escher.BacterialStrategy,sbml2map.BacterialStrategy)
        self.assertEqual(sbml2escher.parser().parse_args([]).program_name,'sbml2map')

    def test_old_and_new_metadata_names_route_identically(self):
        path=Path(__file__).parent/'results/e_coli_core/escher_map.json'
        old=json.loads(path.read_text());new=copy.deepcopy(old)
        context=new[0]['metadata'].pop('sbml2escher')
        new[0]['metadata']['sbml2map']=context
        major={rid for route in context['major_backbones']['routes'] for rid in route['reaction_ids']}
        meta=old[0]['metadata']
        orange={frozenset(meta['segment_source_edges'][sid]) for sid,priority in meta['segment_priorities'].items() if priority==1}
        separate_routes(old,major,orange);separate_routes(new,major,orange)
        self.assertEqual(old[1],new[1])
        self.assertEqual(old[0]['metadata']['overlap_routing'],new[0]['metadata']['overlap_routing'])


if __name__=='__main__':unittest.main()
