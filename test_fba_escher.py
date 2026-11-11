import unittest
import contextlib
import io
import json
import shutil
import tempfile
from pathlib import Path

from model_input import load_model
from run_fba_escher import read_objective, solve_fba, run, read_fba_layout

BASE = Path(__file__).resolve().parent
MODEL = BASE/'models/e_coli_core.xml' if (BASE/'models').is_dir() else BASE/'e_coli_core.xml'
LAYOUT = BASE/'results/e_coli_core' if (BASE/'results/e_coli_core').is_dir() else BASE/'public/results/e_coli_core'


class FbaTests(unittest.TestCase):
    def test_sbml_growth_and_balance(self):
        path = MODEL
        model = load_model(path)
        objective, sense = read_objective(path, {r['id'] for r in model['reactions']})
        fluxes, report = solve_fba(model, objective, sense)
        self.assertAlmostEqual(report['objective_value'], 0.8739215069684304, places=8)
        self.assertLess(report['mass_balance_max_abs'], 1e-9)
        self.assertEqual(report['bound_violation_max'], 0)
        self.assertAlmostEqual(fluxes['EX_glc__D_e'], -10)

    def test_negative_flux_and_infeasible_model(self):
        model = dict(metabolites=[dict(id='a')], reactions=[
            dict(id='uptake', metabolites={'a': -1}, lower_bound=-3, upper_bound=0),
            dict(id='growth', metabolites={'a': -1}, lower_bound=0, upper_bound=10)])
        fluxes, report = solve_fba(model, {'growth': 1})
        self.assertAlmostEqual(fluxes['uptake'], -3)
        self.assertAlmostEqual(report['objective_value'], 3)
        model['reactions'][1]['lower_bound'] = 4
        with self.assertRaisesRegex(RuntimeError, 'FBA failed'):
            solve_fba(model, {'growth': 1})

    def test_full_model_fba_on_filtered_map_keeps_coordinates_and_sampling(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)/'layout'
            shutil.copytree(LAYOUT,directory,ignore=shutil.ignore_patterns('flux_samples.npz','run.log'))
            original_map=(directory/'escher_map.json').read_bytes()
            original_medians=(directory/'escher_flux_medians.json').read_bytes()
            with contextlib.redirect_stdout(io.StringIO()):fluxes,report=run(MODEL,directory)
            self.assertEqual(len(fluxes),95)
            self.assertAlmostEqual(report['objective_value'],0.8739215069684304,places=8)
            display=json.loads((directory/'escher_fluxes.json').read_text())
            data=json.loads((directory/'escher_fba_map.json').read_text())
            self.assertEqual(set(display),{r['bigg_id'] for r in data[1]['reactions'].values()})
            self.assertEqual(len(display),95)
            self.assertEqual(set(display),set(fluxes))
            self.assertTrue(report['map_validation']['schema_valid'])
            self.assertEqual(report['displayed_connections'],360)
            self.assertTrue(report['all_participants_drawn'])
            self.assertNotIn('reaction_coverage',data[0]['metadata'])
            self.assertEqual(report['complete_layout_validation']['existing_nodes_preserved'],206)
            self.assertEqual((directory/'escher_map.json').read_bytes(),original_map)
            self.assertEqual((directory/'escher_flux_medians.json').read_bytes(),original_medians)
            self.assertEqual(set(json.loads((directory/'fba_fluxes.json').read_text())),set(fluxes))
            self.assertTrue(report['source_model_used_for_fba'])
            self.assertFalse(report['coordinates_changed'])
            self.assertEqual(report['undisplayed_active_reactions'],{})
            self.assertEqual(report['undisplayed_reactions'],[])
            original_nodes=json.loads(original_map)[0]['metadata']['source_node_ids']
            new_nodes=data[0]['metadata']['source_node_ids']
            old_body=json.loads(original_map)[1]['nodes']
            new_by_key={key:data[1]['nodes'][nid] for nid,key in new_nodes.items()}
            for nid,key in original_nodes.items():
                # Original grid coordinates remain fixed; display origin may translate.
                old_transform=json.loads(original_map)[0]['metadata']['coordinate_transform']
                new_transform=data[0]['metadata']['coordinate_transform']
                for axis,origin,sign in [('x','x_min',1),('y','y_max',-1)]:
                    shift=sign*old_transform['scale']*(old_transform[origin]-new_transform[origin])
                    self.assertAlmostEqual(new_by_key[key][axis],old_body[nid][axis]+shift)
            self.assertIn('highlight-panel',(directory/'escher_fba.html').read_text())
            self.assertIn('FBA uses all 95 source reactions',(directory/'escher_fba.html').read_text())
            self.assertIn('fba-result-link',(directory/'index.html').read_text())
            for rid, value in display.items():
                multiplier=data[0]['metadata']['reactions'][rid]['coefficient_multiplier']
                self.assertAlmostEqual(value,multiplier*fluxes[rid],places=7)

    def test_modified_display_chemistry_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)/'layout'
            shutil.copytree(LAYOUT,directory,ignore=shutil.ignore_patterns('flux_samples.npz','run.log'))
            path=directory/'display_model.json';model=json.loads(path.read_text())
            r=model['reactions'][0];mid=next(iter(r['metabolites']));r['metabolites'][mid]*=2
            path.write_text(json.dumps(model))
            with self.assertRaisesRegex(ValueError,'chemistry or bounds'):
                read_fba_layout(directory,load_model(MODEL))

    def test_different_source_bounds_are_rejected(self):
        model=load_model(MODEL);model['reactions'][0]['upper_bound']+=1
        with self.assertRaisesRegex(ValueError,'does not match'):
            read_fba_layout(LAYOUT,model)


if __name__ == '__main__':
    unittest.main()
