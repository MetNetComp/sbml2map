"""Boundary reservoirs are excluded from steady-state balance, never from chemistry."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from escher_grid import model_hash, sampling_data
from flux_sampling2 import sample_fluxes, sampling_constraints
from model_input import load_model
from sbml2map import parser, run
from stoichiometric_network import stoichiometric_matrix, steady_state_matrix


def reservoir_model():
    return dict(id='reservoir', metabolites=[dict(id='external', boundaryCondition=True),
                dict(id='a_c', formula='C2H4'), dict(id='h_c', formula='H')], reactions=[
        dict(id='EX_in', metabolites={'external': -1, 'a_c': 1, 'h_c': -1}, lower_bound=2, upper_bound=5),
        dict(id='out', metabolites={'a_c': -1, 'h_c': 1}, lower_bound=0, upper_bound=5)])


class BoundarySamplingTests(unittest.TestCase):
    def test_loader_keeps_explicit_boundary_and_defaults_to_internal(self):
        xml='''<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core"><model id="test">
        <listOfSpecies><species id="M_external" compartment="BOUNDARY" boundaryCondition="true"/>
        <species id="M_a_c" compartment="c" boundaryCondition="false"/></listOfSpecies>
        <listOfReactions><reaction id="R_in"><listOfReactants><speciesReference species="M_external"/></listOfReactants>
        <listOfProducts><speciesReference species="M_a_c"/></listOfProducts></reaction></listOfReactions>
        </model></sbml>'''
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'input.xml';path.write_text(xml);model=load_model(path)
        self.assertIs(model['metabolites'][0]['boundaryCondition'], True)
        self.assertFalse(model['metabolites'][1].get('boundaryCondition', False))
        self.assertEqual(steady_state_matrix(model)[1], ['a_c'])

    def test_nonzero_exchange_sampling_keeps_internal_proton_balance(self):
        model=reservoir_model();original=copy.deepcopy(model)
        fluxes, stats, report=sample_fluxes(model, samples=1000, chains=4, thinning=2, burn_in=50,
                                          seed=4, progress=lambda *a, **k: None)
        np.testing.assert_allclose(fluxes[:,:,0], fluxes[:,:,1], atol=1e-8)
        self.assertGreaterEqual(fluxes.min(), 2-1e-8)
        self.assertLessEqual(fluxes.max(), 5+1e-8)
        self.assertAlmostEqual(stats[0]['median'], 3.5, delta=.15)
        self.assertEqual(report['balanced_metabolite_ids'], ['a_c','h_c'])
        self.assertEqual(report['boundary_metabolite_ids'], ['external'])
        self.assertTrue(report['performed'])
        self.assertTrue(report['mixing_checks_passed'])
        self.assertEqual(model, original)
        full,_,_=stoichiometric_matrix(model)
        self.assertGreater(np.max(np.abs(fluxes@np.asarray(full).T)), 2)

    def test_names_and_drawing_flags_do_not_exempt_internal_species(self):
        model=reservoir_model();model['metabolites'][0].pop('boundaryCondition')
        model['metabolites'][0].update(compartment='BOUNDARY', _boundary_placeholder=True)
        self.assertEqual(len(steady_state_matrix(model)[1]), 3)
        with self.assertRaisesRegex(ValueError, 'LP failed'):
            sample_fluxes(model, samples=20, chains=2, thinning=1, burn_in=1)
        model['metabolites'][0]['boundaryCondition']='false'
        with self.assertRaisesRegex(ValueError, 'must be boolean'):
            sampling_constraints(model)

    def test_saved_samples_require_matching_boundary_constraints(self):
        model=reservoir_model();digest=model_hash(model)
        args=SimpleNamespace(reuse_samples=None, samples=20, chains=2, thinning=1, burn_in=1, seed=0)
        samples, stats, report=sampling_data(args, model, digest)
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);args.reuse_samples=root
            np.savez_compressed(root/'flux_samples.npz', fluxes=samples, reaction_ids=[r['id'] for r in model['reactions']])
            report['samples_sha256']=hashlib.sha256((root/'flux_samples.npz').read_bytes()).hexdigest()
            (root/'flux_statistics.json').write_text(json.dumps(stats))
            def save(): (root/'sampling_report.json').write_text(json.dumps(report))
            save();reused,_,metadata=sampling_data(args,model,digest)
            np.testing.assert_array_equal(reused,samples)
            self.assertEqual(metadata['mode'],'reused')
            report['constraints_sha256']='wrong';save()
            with self.assertRaisesRegex(ValueError,'constraint hash differs'):
                sampling_data(args,model,digest)
            del report['constraints_sha256'];save()
            with self.assertRaisesRegex(ValueError,'require a constraint fingerprint'):
                sampling_data(args,model,digest)

    def test_require_sampling_stops_without_a_structural_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);path=root/'model.json';path.write_text(json.dumps(reservoir_model()))
            args=parser().parse_args(['--model',str(path),'--output-dir',str(root/'output'),'--require-sampling'])
            with patch('escher_grid.sample_fluxes', side_effect=ValueError('infeasible')):
                with self.assertRaisesRegex(RuntimeError,'Required flux sampling failed: infeasible'):
                    run(args)
            self.assertFalse((args.output_dir/'index.html').exists())
            args.discover_only=True
            with self.assertRaisesRegex(ValueError,'cannot be combined'):
                run(args)
            args.discover_only=False
            args.no_sampling=True
            with self.assertRaisesRegex(ValueError,'cannot be combined'):
                run(args)


if __name__=='__main__':unittest.main()
