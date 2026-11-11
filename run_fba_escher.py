#!/usr/bin/env python3
"""Solve SBML FBA with SciPy/HiGHS and overlay signed fluxes on an Escher map."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
import scipy
from scipy.optimize import linprog

from export_escher import read_layout, write_preview, validate_map, create_map
from reaction_network import build_network, validate_network
from highlighting import inject_controls, grid_elements
from full_fba_layout import complete_layout
from model_input import load_model
from stoichiometric_network import stoichiometric_matrix

BASE = Path(__file__).resolve().parent
TOLERANCE = 1e-7


def read_objective(path, reaction_ids):
    root = ET.parse(path).getroot()
    # Boundary species require a different balance matrix; reject rather than misinterpret.
    if any(e.get('boundaryCondition') in ('true', '1') for e in root.iter()
           if e.tag.split('}')[-1] == 'species'):
        raise ValueError('Boundary species are unsupported by this FBA runner')
    lists = [e for e in root.iter() if e.tag.split('}')[-1] == 'listOfObjectives']
    if len(lists) != 1:
        raise ValueError('Expected one FBC objective list')
    def attr(e, name):
        return next((v for k, v in e.attrib.items() if k.split('}')[-1] == name), None)
    active = attr(lists[0], 'activeObjective')
    objectives = [e for e in lists[0] if attr(e, 'id') == active]
    if len(objectives) != 1 or attr(objectives[0], 'type') not in ('maximize', 'minimize'):
        raise ValueError('Missing or invalid active FBC objective')
    weights = {}
    for e in objectives[0].iter():
        if e.tag.split('}')[-1] != 'fluxObjective':
            continue
        raw = attr(e, 'reaction')
        rid = re.sub(r'__(\d+)__', lambda m: chr(int(m[1])), raw.removeprefix('R_'))
        if rid not in reaction_ids:
            raise ValueError('Unknown objective reaction: ' + rid)
        value = float(attr(e, 'coefficient'))
        if not np.isfinite(value):
            raise ValueError('Nonfinite objective coefficient')
        weights[rid] = weights.get(rid, 0) + value
    if not any(weights.values()):
        raise ValueError('Empty objective')
    return weights, attr(objectives[0], 'type')


def solve_fba(model, objective, sense='maximize'):
    matrix, mids, rids = stoichiometric_matrix(model)
    matrix = np.asarray(matrix)
    weights = np.array([objective.get(rid, 0) for rid in rids], dtype=float)
    bounds = [(r['lower_bound'], r['upper_bound']) for r in model['reactions']]
    result = linprog(-weights if sense == 'maximize' else weights,
                     A_eq=matrix, b_eq=np.zeros(len(mids)), bounds=bounds, method='highs',
                     options={'primal_feasibility_tolerance': 1e-9,
                              'dual_feasibility_tolerance': 1e-9})
    if not result.success:
        raise RuntimeError('FBA failed: ' + result.message)
    fluxes = result.x
    balance = float(np.max(np.abs(matrix @ fluxes)))
    violation = max([0.0] + [max(lo - v, v - hi) for v, (lo, hi) in zip(fluxes, bounds)])
    if not np.all(np.isfinite(fluxes)) or balance > TOLERANCE or violation > TOLERANCE:
        raise RuntimeError('FBA balance or bound verification failed')
    report = dict(status='optimal', method='FBA (linear programming)', solver='SciPy linprog / HiGHS',
                  scipy_version=scipy.__version__, objective=objective, objective_sense=sense,
                  objective_value=float(weights @ fluxes), flux_units='mmol/gDW/h',
                  metabolites=len(mids), reactions=len(rids), mass_balance_max_abs=balance,
                  bound_violation_max=violation, verification_tolerance=TOLERANCE,
                  bounds_modified=False, active_reactions=sum(abs(v) > TOLERANCE for v in fluxes),
                  limitation='One optimal FBA solution; alternative optima may have different fluxes. Not pFBA.')
    report['active_reactions'] = int(report['active_reactions'])
    return dict(zip(rids, map(float, fluxes))), report


def read_fba_layout(directory, model):
    """Read either a full legacy layout or a filtered sbml2map display map."""
    if not (directory / 'display_model.json').exists():
        layout_model, nodes, edges, xy = read_layout(directory)
        if model != layout_model:
            raise ValueError('SBML model does not match the saved layout model')
        return layout_model, nodes, edges, xy
    saved = json.loads((directory / 'source_model.json').read_text())
    if model != saved:
        raise ValueError('SBML model does not match the saved layout model')
    layout_model = json.loads((directory / 'display_model.json').read_text())
    original = {r['id']: r for r in model['reactions']}
    if len({r['id'] for r in layout_model['reactions']}) != len(layout_model['reactions']):
        raise ValueError('Duplicate display reaction IDs')
    for reaction in layout_model['reactions']:
        source = original.get(reaction['id'])
        if source is None or {k:v for k,v in reaction.items() if k!='metabolites'} != {
                k:v for k,v in source.items() if k!='metabolites'} or any(
                source['metabolites'].get(mid)!=value for mid,value in reaction['metabolites'].items()):
            raise ValueError('Display model changes source reaction chemistry or bounds')
    with (directory / 'coordinates.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    currency = {row['original_id'] for row in rows if row['kind']=='currency'}
    nodes, edges = build_network(layout_model, currency)
    if len(rows)!=len(nodes) or {row['id'] for row in rows}!=set(nodes):
        raise ValueError('Display coordinates do not match the saved model')
    for row in rows:
        if any(row[k]!=nodes[row['id']][k] for k in ('kind','original_id')):
            raise ValueError('Display coordinate identity mismatch')
    xy = {row['id']: [float(row['x']),float(row['y'])] for row in rows}
    validate_network(layout_model,nodes,edges,xy)
    return layout_model,nodes,edges,xy


def run(model_path, directory):
    model_path, directory = Path(model_path), Path(directory)
    model = load_model(model_path)
    layout_model, nodes, edges, xy = read_fba_layout(directory, model)
    data = json.loads((directory / 'escher_map.json').read_text())
    validate_map(data, layout_model, nodes, xy)
    transform = data[0]['metadata']['coordinate_transform']
    label_options = data[0]['metadata'].get('label_display', {})
    nodes, edges, xy, layout_validation = complete_layout(model, nodes, xy)
    data = create_map(model, nodes, edges, xy, scale=transform['scale'],
                      padding=transform['padding'], layout_method='complete FBA layout')
    data[0]['metadata']['label_display'] = {
        key: value for key, value in label_options.items() if isinstance(value, bool)}
    routing_optimization = None
    saved_routing = directory / 'routing_geometry.json'
    if saved_routing.exists():
        from routing_geometry import build_routes
        from routing_optimizer import optimize
        from escher_routing import separate_routes
        saved = json.loads(saved_routing.read_text())
        settings = saved['optimization']
        # Full FBA adds reaction-local participants. Reuse all matching structural
        # incidences exactly, and optimize added branches without moving nodes.
        records = [dict(reaction_id=r['id'], selected=True) for r in model['reactions']]
        paths, _ = build_routes(nodes, edges, xy, records, {'mode': 'structure'}, route_offset=0)
        old = {(p['rid'], frozenset((p['source'], p['target']))): p for p in saved['routes']}
        for path in paths:
            parent = old.get((path['rid'], frozenset((path['source'], path['target']))))
            if parent:
                points = list(parent['points'])
                if parent['source'] != path['source']:
                    points.reverse()
                path.update(points=points, routing_locked=True)
        paths, _, routing_optimization = optimize(paths, xy, mode=settings['objective'],
            crossing_weight=settings['crossing_weight'], bend_weight=settings['bend_weight'],
            time_limit=settings['time_limit_seconds'], node_moves=False)
        separate_routes(data, enabled=False, route_geometry=paths,
                        routing_report={'optimization': routing_optimization, 'routing_preference': settings['objective']})
    map_validation = validate_map(data, model, nodes, xy)
    mapped = {r['bigg_id']: r for r in data[1]['reactions'].values()}
    # Balance and optimize the full source model, including invisible participants.
    objective, sense = read_objective(model_path, {r['id'] for r in model['reactions']})
    fluxes, report = solve_fba(model, objective, sense)
    display = {}
    for reaction in model['reactions']:
        rid = reaction['id']
        multiplier = data[0]['metadata']['reactions'][rid]['coefficient_multiplier']
        expected = -1 if reaction['upper_bound'] <= 0 and reaction['lower_bound'] < 0 else 1
        if multiplier != expected or {m['bigg_id']: m['coefficient'] for m in mapped[rid]['metabolites']} != {
                mid: multiplier * c for mid, c in reaction['metabolites'].items()}:
            raise ValueError('Map stoichiometry/sign mismatch: ' + rid)
        value = multiplier * fluxes[rid]
        display[rid] = 0.0 if abs(value) <= TOLERANCE else value
    if routing_optimization is not None:
        report['routing_optimization'] = routing_optimization
    omitted = sorted(set(fluxes)-set(display))
    report.update(source_sbml=str(model_path), source_sha256=hashlib.sha256(model_path.read_bytes()).hexdigest(),
                  map_id=data[0]['map_id'], display_zero_tolerance=TOLERANCE,
                  source_model_used_for_fba=True,
                  coordinates_changed=layout_validation['existing_coordinates_changed'],
                  complete_layout_validation=layout_validation, all_participants_drawn=True,
                  displayed_connections=len(edges), pathway_highlighting=False,
                  displayed_reactions=len(display), map_validation=map_validation,
                  undisplayed_reactions=omitted,
                  undisplayed_active_reactions={rid:fluxes[rid] for rid in omitted if abs(fluxes[rid])>TOLERANCE},
                  display_sign_multipliers={rid: d['coefficient_multiplier'] for rid, d in data[0]['metadata']['reactions'].items()
                                            if d['coefficient_multiplier'] != 1},
                  browser_validation='not run')
    with (directory / 'fba_fluxes.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(['reaction_id', 'flux', 'escher_flux', 'lower_bound', 'upper_bound'])
        for r in model['reactions']:
            writer.writerow([r['id'], fluxes[r['id']], display.get(r['id'], ''), r['lower_bound'], r['upper_bound']])
    for name, content in [('fba_fluxes.json', fluxes), ('escher_fluxes.json', display), ('fba_report.json', report)]:
        (directory / name).write_text(json.dumps(content, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    data[0]['map_name'] = model['id'] + ': FBA with all reactions and connections'
    data[0]['map_description'] = 'All source reactions and participants, including zero flux. Existing pathway coordinates retained.'
    data[0]['metadata']['fba'] = report
    (directory / 'escher_fba_map.json').write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    write_preview(directory / 'escher_fba.html', data, display, report)
    page = (directory / 'escher_fba.html').read_text()
    connections = sum(len(r['metabolites']) for r in model['reactions'])
    notice = (f'<p>FBA uses all {len(fluxes)} source reactions. '
              f'Displaying {len(display)}/{len(fluxes)} reactions and {len(edges)}/{connections} connections, including zero flux. '
              'Water, protons and cofactors are included as reaction-local copies. '
              '<a href="escher_fba_map.json">Complete FBA map (JSON)</a> · '
              '<a href="fba_fluxes.csv">All fluxes (CSV)</a> · <a href="fba_fluxes.json">All fluxes (JSON)</a> · '
              '<a href="escher_fluxes.json">Escher reaction data</a> · <a href="fba_report.json">FBA report</a> · '
              '<a href="index.html">Original drawing</a></p>')
    page = page.replace('<header>', '<header>'+notice, 1)
    (directory / 'escher_fba.html').write_text(inject_controls(page, model, grid_elements(nodes, edges), escher=True))
    for filename in ('index.html','escher_preview.html'):
        path=directory/filename
        if path.exists():
            page=path.read_text()
            if 'id="fba-result-link"' not in page:
                page=page.replace('<header>', '<header><p id="fba-result-link"><a href="escher_fba.html">View FBA result</a></p>',1)
                path.write_text(page)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return fluxes, report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, default=BASE / ('models/e_coli_core.xml' if (BASE/'models').is_dir() else 'e_coli_core.xml'))
    parser.add_argument('--layout-dir', type=Path, default=BASE / ('results/e_coli_core' if (BASE/'results').is_dir() else 'public/results/e_coli_core'))
    args = parser.parse_args()
    run(args.model, args.layout_dir)
