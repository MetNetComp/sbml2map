#!/usr/bin/env python3
"""Sampled or structural SBML maps with inorganic species hidden and grid placement."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from routing_simplification import POLICY,simplify_important,yield_gray

from red_route_layout import validate_red_spacing,plan_red_routes,spacing_report

from parallel_routing import separate_parallel, avoid_nodes, validate_route_offset, structural_node_priorities

from highlighting import HIGHLIGHT_COLOR, read_selection, grid_elements, catalog, inject_controls

from flux_sampling2 import sample_fluxes, sampling_constraints
from annotation_context import hidden_metabolite
from grid_layout2 import optimize_grid
from major_backbones2 import RED, backbone_plan, pathway_definitions, validate_backbones
from model_input import load_model
from reaction_network import build_network, validate_network
from stoichiometric_network import stoichiometric_matrix
from exchange_routes import ORANGE, shared_coenzymes, connector_plan, connector_length, connector_cost, orange_flux_costs

BASE = Path(__file__).resolve().parent
# Explicit biochemical display policy; acyl-CoA metabolites remain shared grid nodes.
COENZYMES = {'atp', 'adp', 'amp', 'gtp', 'gdp', 'gmp', 'utp', 'udp', 'ump',
             'ctp', 'cdp', 'cmp', 'nad', 'nadh', 'nadp', 'nadph', 'fad', 'fadh2',
             'fmn', 'fmnh2', 'coa', 'q8', 'q8h2', 'q10', 'q10h2', 'mqn8', 'mql8',
             'thf', 'dhf', '10fthf', '5mthf', 'mlthf', 'pydx5p', 'btn', 'thmpp'}
INORGANIC_IDS = {'h', 'h2o', 'o2', 'co2', 'co', 'hco3', 'co3', 'nh4', 'nh3',
                 'pi', 'ppi', 'so4', 'so3', 'h2s', 'no2', 'no3', 'h2', 'n2',
                 'na1', 'k', 'cl', 'ca2', 'mg2', 'fe2', 'fe3', 'cu2', 'zn2'}
INORGANIC_CARBON_FORMULAS = {'CO', 'CO2', 'CO3', 'CHO3', 'HCO3', 'CH2O3', 'H2CO3'}


def dump(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def model_hash(model):
    return hashlib.sha256(json.dumps(model, sort_keys=True).encode()).hexdigest()


def classify_metabolites(model, override=None):
    overrides = override or {}
    ids = {m['id'] for m in model['metabolites']}
    if set(overrides) - ids or set(overrides.values()) - {'inorganic', 'coenzyme', 'organic'}:
        raise ValueError('Classification overrides must map known metabolite IDs to inorganic/coenzyme/organic')
    result = {}
    for m in model['metabolites']:
        mid = m['id']
        base = re.sub(r'_[a-z][0-9]?$', '', mid).lower()
        formula = m.get('chemicalFormula', m.get('formula', '')) or ''
        elements = set(re.findall(r'[A-Z][a-z]?', formula))
        if mid in overrides:
            kind, reason = overrides[mid], 'explicit override'
        elif base in INORGANIC_IDS or formula in INORGANIC_CARBON_FORMULAS:
            kind, reason = 'inorganic', 'known inorganic ID or carbon oxide/carbonate formula'
        elif elements and 'C' not in elements:
            kind, reason = 'inorganic', 'chemical formula contains no carbon'
        elif base in COENZYMES:
            kind, reason = 'coenzyme', 'explicit coenzyme/nucleotide/electron-carrier ID list'
        else:
            kind, reason = 'organic', 'retained; no inorganic/coenzyme rule matched'
        result[mid] = dict(kind=kind, reason=reason, name=m.get('name', ''), formula=formula,
                           formula_missing=not bool(formula))
    return result


def select_reactions(model, stats, classification, coverage, zero_tolerance=1e-7, major_reactions=()):
    if not math.isfinite(coverage) or not 0 <= coverage <= 100:
        raise ValueError('coverage must be between 0 and 100')
    if not math.isfinite(zero_tolerance) or zero_tolerance <= 0:
        raise ValueError('zero-tolerance must be positive and finite')
    major_reactions = set(major_reactions)
    lookup = {r['id']: r for r in model['reactions']}
    records, eligible = [], []
    if {s['reaction_id'] for s in stats} != set(lookup) or len(stats) != len(lookup):
        raise ValueError('Statistics must contain each source reaction exactly once')
    for stat in stats:
        r = lookup[stat['reaction_id']]
        median = float(stat['median'])
        if not math.isfinite(median):
            raise ValueError('Nonfinite median')
        visible = [mid for mid in r['metabolites'] if not hidden_metabolite(classification[mid])]
        record = dict(stat, priority=abs(median), selected=False, rank=None, median_rank=None,
                      is_major_pathway=r['id'] in major_reactions,
                      visible_metabolites=visible, omitted_inorganic=[m for m in r['metabolites'] if m not in visible])
        record['exclusion_reason'] = ('zero_median' if abs(median) <= zero_tolerance else
                                      'inorganic_only' if not visible else 'below_coverage_cutoff')
        omitted_placeholders = [m for m in r['metabolites'] if classification[m].get('boundary_placeholder')]
        if omitted_placeholders:
            record['omitted_boundary_placeholders'] = omitted_placeholders
            record['omitted_inorganic'] = [m for m in record['omitted_inorganic'] if m not in omitted_placeholders]
        records.append(record)
        if abs(median) > zero_tolerance and visible:
            eligible.append(record)
    eligible.sort(key=lambda r: (-r['priority'], r['reaction_id']))
    for rank, record in enumerate(eligible, 1):
        record['median_rank'] = rank
    eligible.sort(key=lambda r: (not r['is_major_pathway'], -r['priority'], r['reaction_id']))
    count = math.ceil(coverage * len(eligible) / 100)
    for rank, record in enumerate(eligible, 1):
        record['rank'] = rank
        if rank <= count:
            record.update(selected=True, exclusion_reason='')
    selected = {r['reaction_id'] for r in eligible[:count]}
    return selected, records, dict(coverage_percent=coverage, eligible_reactions=len(eligible),
                                  selected_reactions=len(selected), zero_tolerance=zero_tolerance,
                                  ranking='Major pathways first, then descending abs(signed median flux) within each group; reaction ID breaks ties',
                                  major_eligible_reactions=sum(r['is_major_pathway'] for r in eligible),
                                  major_selected_reactions=sum(r['is_major_pathway'] for r in eligible[:count]),
                                  count_rule='ceil(coverage / 100 * eligible_count); 0 gives empty map')


def display_graph(model, selected, classification, required_metabolites=(), shared_metabolites=()):
    reactions = [dict(r, metabolites={m: c for m, c in r['metabolites'].items()
                                     if not hidden_metabolite(classification[m])})
                 for r in model['reactions'] if r['id'] in selected]
    used = {m for r in reactions for m in r['metabolites']}
    used.update(required_metabolites)
    projected = dict(model, metabolites=[m for m in model['metabolites'] if m['id'] in used], reactions=reactions)
    coenzymes = {m for m in used if classification[m]['kind'] == 'coenzyme'} - set(shared_metabolites)
    nodes, edges = build_network(projected, coenzymes)
    return projected, nodes, edges


def select_structural_reactions(model, classification, coverage, major_reactions=()):
    """Select by pathway membership and ID, without inventing flux statistics."""
    if not math.isfinite(coverage) or not 0 <= coverage <= 100:
        raise ValueError('coverage must be between 0 and 100')
    major_reactions = set(major_reactions)
    records, eligible = [], []
    for reaction in model['reactions']:
        visible = [mid for mid in reaction['metabolites'] if not hidden_metabolite(classification[mid])]
        record = dict(reaction_id=reaction['id'], selected=False, rank=None,
                      is_major_pathway=reaction['id'] in major_reactions,
                      visible_metabolites=visible,
                      omitted_inorganic=[mid for mid in reaction['metabolites'] if mid not in visible],
                      lower_bound=reaction['lower_bound'], upper_bound=reaction['upper_bound'],
                      exclusion_reason='below_coverage_cutoff' if visible else 'inorganic_only')
        omitted_placeholders = [m for m in reaction['metabolites'] if classification[m].get('boundary_placeholder')]
        if omitted_placeholders:
            record['omitted_boundary_placeholders'] = omitted_placeholders
            record['omitted_inorganic'] = [m for m in record['omitted_inorganic'] if m not in omitted_placeholders]
        records.append(record)
        if visible:
            eligible.append(record)
    eligible.sort(key=lambda record: (not record['is_major_pathway'], record['reaction_id']))
    count = math.ceil(coverage * len(eligible) / 100)
    for rank, record in enumerate(eligible, 1):
        record['rank'] = rank
        if rank <= count:
            record.update(selected=True, exclusion_reason='')
    selected = {record['reaction_id'] for record in eligible[:count]}
    return selected, records, dict(mode='structure', coverage_percent=coverage,
                                  eligible_reactions=len(eligible), selected_reactions=len(selected),
                                  major_eligible_reactions=sum(r['is_major_pathway'] for r in eligible),
                                  major_selected_reactions=sum(r['is_major_pathway'] for r in eligible[:count]),
                                  ranking='Major pathways first, then reaction ID; no flux or activity ranking',
                                  count_rule='ceil(coverage / 100 * eligible_count); 0 gives empty map')


def attach_coenzymes(nodes, xy, radius=.45, grid_shape=None, forbidden_region=None, boundary=None):
    groups = {}
    for key, node in nodes.items():
        if node['kind'] == 'currency':
            groups.setdefault(node['parent_id'], []).append(key)
    for parent, keys in sorted(groups.items()):
        if grid_shape is None and forbidden_region is None and boundary is None:
            angles = 2 * math.pi * (np.arange(len(keys)) + .25) / len(keys)
        else:
            angles = 2 * math.pi * (np.arange(4096) + .25) / 4096
        vectors = np.column_stack([np.cos(angles), np.sin(angles)])
        vectors *= radius / np.abs(vectors).sum(axis=1, keepdims=True)
        candidates = np.array(xy[parent]) + vectors
        if grid_shape is not None:
            cols, rows = grid_shape
            valid = candidates[(candidates[:, 0] >= 0) & (candidates[:, 0] <= cols - 1) &
                               (candidates[:, 1] >= 0) & (candidates[:, 1] <= rows - 1)]
            if len(valid) < len(keys):
                raise ValueError('Grid cannot contain off-grid coenzyme copies; use at least 2 rows and columns')
            candidates = valid[np.arange(len(keys)) * len(valid) // len(keys)]
        if boundary is not None:
            xmin,xmax,ymin,ymax=boundary
            candidates=np.asarray(xy[parent])+vectors
            mask=(candidates[:,0]>=xmin)&(candidates[:,0]<=xmax)&(candidates[:,1]>=ymin)&(candidates[:,1]<=ymax)
            if grid_shape is not None:mask &= (candidates[:,0]>=0)&(candidates[:,0]<=cols-1)&(candidates[:,1]>=0)&(candidates[:,1]<=rows-1)
            candidates=candidates[mask]
            if forbidden_region is None and len(candidates)>=len(keys):candidates=candidates[np.arange(len(keys))*len(candidates)//len(keys)]
        if forbidden_region is not None:
            from tca_geometry import allowed_points
            valid = candidates[allowed_points(candidates, forbidden_region)]
            if len(valid) < len(keys):
                # Recompute all outward candidates; keep the fixed L1 cofactor radius.
                angles = 2*math.pi*(np.arange(4096)+.25)/4096
                vectors = np.column_stack([np.cos(angles),np.sin(angles)])
                vectors *= radius/np.abs(vectors).sum(axis=1,keepdims=True)
                valid = np.asarray(xy[parent])+vectors
                mask = allowed_points(valid,forbidden_region)
                if grid_shape is not None:
                    cols,rows=grid_shape
                    mask &= (valid[:,0]>=0)&(valid[:,0]<=cols-1)&(valid[:,1]>=0)&(valid[:,1]<=rows-1)
                if boundary is not None:
                    mask &= (valid[:,0]>=xmin)&(valid[:,0]<=xmax)&(valid[:,1]>=ymin)&(valid[:,1]<=ymax)
                valid = valid[mask]
            if len(valid) < len(keys):
                raise ValueError('No cofactor position outside TCA; enlarge the grid')
            candidates = valid[np.arange(len(keys))*len(valid)//len(keys)]
        if len(candidates)<len(keys):raise ValueError('No cofactor space inside exchange boundary')
        for key, position in zip(sorted(keys), candidates):
            xy[key] = position.tolist()


def separate_static_labels(ax, renderer):
    """Pack labels in screen coordinates without moving scientific nodes."""
    from matplotlib.transforms import Bbox
    labels=sorted([a for a in ax.texts if (a.get_gid() or '').startswith('label-')],key=lambda a:-a.get_zorder())
    occupied=[];moved=0
    for text in labels:
        if not hasattr(text,'_preferred_label_offset'):text._preferred_label_offset=text.get_position()
        preferred=text._preferred_label_offset
        candidates=[preferred]
        for distance in range(8,97,8):
            candidates.extend([(distance,5),(-distance,5),(4,distance),(4,-distance),
                               (distance,distance),(-distance,distance),(distance,-distance),(-distance,-distance)])
        for offset in candidates:
            text.set_position(offset);box=text.get_window_extent(renderer)
            padded=Bbox.from_extents(box.x0-1,box.y0-1,box.x1+1,box.y1+1)
            if not any(padded.overlaps(other) for other in occupied):break
        occupied.append(padded);moved+=offset!=preferred
    overlaps=int(sum(a.overlaps(b) for i,a in enumerate(occupied) for b in occupied[i+1:]))
    return dict(labels=len(labels),moved=moved,overlaps=overlaps)


def render(output, nodes, edges, xy, records, selection, model_id, backbones=None, program_name='sbml2map', sampling_report=None, highlights=None, source_model=None, route_offset=.3, label_options=None):
    highlights = highlights or dict(metabolites=[], reactions=[])
    label_options={key:bool((label_options or {}).get(key,default)) for key,default in
        [('show_metabolites',True),('show_reactions',True),('allow_overlap',False)]}
    elements = grid_elements(nodes, edges)
    if source_model is None:
        source_model = {kind: [dict(id=key, name=next(n.get('name', '') for n in nodes.values() if n['original_id']==key))
                             for key in sorted({e['key'] for e in elements if e['kind']==kind})]
                        for kind in ('metabolites', 'reactions')}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as path_effects
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    stats = {r['reaction_id']: r for r in records}
    structural = selection.get('mode') == 'structure'
    major_mids = {mid for route in (backbones or {}).get('routes', []) for mid in route['metabolite_ids']}
    connectors = (backbones or {}).get('exchange_connectors')
    orange_edges = set(connectors['edge_indices']) if connectors else set()
    orange_nodes = {key for i in orange_edges for key in (edges[i]['source'], edges[i]['target'])}
    orange_pairs={frozenset((edges[i]['source'],edges[i]['target'])) for i in orange_edges}
    node_priorities=structural_node_priorities(nodes,
        {rid for rid,record in stats.items() if record.get('is_major_pathway',False)},major_mids,orange_pairs)
    obstacle_keys=sorted(key for key,rank in node_priorities.items() if rank>0)
    obstacles=[xy[key] for key in obstacle_keys]
    selected = [r for r in records if r['selected']]
    from display_ranking import reaction_display_order
    orange_reactions = {edges[i]['reaction_id']: 1 for i in orange_edges}
    display_ranks = {rid: rank for rank, rid in enumerate(
        reaction_display_order(selected, orange_reactions), 1)}
    priorities = [] if structural else [r['priority'] for r in selected if not r.get('is_major_pathway', False)]
    vmin, vmax = (min(priorities), max(priorities)) if priorities else (1, 10)
    vmax = max(vmax, vmin * 1.01)
    norm = LogNorm(vmin=max(vmin, 1e-12), vmax=vmax)
    cmap = plt.get_cmap('viridis')
    span = np.ptp(np.array(list(xy.values())), axis=0) if xy else np.array([8, 8])
    frame=(backbones or {}).get('render_frame')
    if frame:span=np.array([frame[1]-frame[0],frame[3]-frame[2]])
    size = (max(10, min(22, .8 * span[0] + 3)), max(10, min(30, .8 * span[1] + 3)))
    fig, ax = plt.subplots(figsize=size)
    ax.set_aspect('equal')
    if xy:
        arr = np.array(list(xy.values()))
        ax.set_xticks(np.arange(math.floor(arr[:, 0].min()), math.ceil(arr[:, 0].max()) + 1))
        ax.set_yticks(np.arange(math.floor(arr[:, 1].min()), math.ceil(arr[:, 1].max()) + 1))
    ax.grid(color='#d7dee5', alpha=.55, linewidth=.5, zorder=0)
    # Defer autoscaling until all artists exist; scatter otherwise rescans
    # previous collections for every added node on genome-scale maps.
    ax.set_autoscale_on(False)
    from tca_curves import curved_pairs, arc_controls, sample_curve
    tca_pairs=curved_pairs(backbones)
    tca_center=((backbones or {}).get('tca_exclusion',{}).get('region') or {}).get('center')
    region=((backbones or {}).get('tca_exclusion',{}).get('region') or {})
    if region.get('shape')=='circle' and region.get('node_grid_step',1)<1:
        # Show the declared half-step lattice locally, keeping other regions clear.
        step=region['node_grid_step'];cx,cy=region['center'];radius=region['radius']
        low=[math.floor((v-radius)/step)*step for v in (cx,cy)]
        high=[math.ceil((v+radius)/step)*step for v in (cx,cy)]
        xs=[v for v in np.arange(low[0],high[0]+step/2,step) if abs(v-round(v))>1e-8]
        ys=[v for v in np.arange(low[1],high[1]+step/2,step) if abs(v-round(v))>1e-8]
        ax.vlines(xs,low[1],high[1],color='#e6ebef',linewidth=.4,zorder=0)
        ax.hlines(ys,low[0],high[0],color='#e6ebef',linewidth=.4,zorder=0)

    backbone_pairs=set()
    for route in (backbones or {}).get('routes', []):
        for step in route.get('steps', []):
            for mid in (step['source'],step['target']):
                backbone_pairs.add(frozenset(('rxn:'+step['reaction_id'],'met:'+mid)))
    display_paths=[]
    for j, edge in enumerate(edges):
        rid=edge['reaction_id']
        a,b=tuple(xy[edge['source']]),tuple(xy[edge['target']])
        reverse=edge['direction']=='reverse' if structural else stats[rid]['median']<0
        if reverse:a,b=b,a
        elbow=(b[0],a[1])
        points=[a]+([elbow] if elbow!=a and elbow!=b else [])+[b]
        display_paths.append(dict(rid=rid,points=points,
            priority=2 if stats[rid].get('is_major_pathway',False) else 1 if j in orange_edges else 0,
            backbone=frozenset((edge['source'],edge['target'])) in backbone_pairs))
    original_route_length=sum(sum(abs(a[k]-b[k]) for k in (0,1)) for p in display_paths for a,b in zip(p['points'],p['points'][1:]))
    canonical=(backbones or {}).get('crossing_routing')
    red_plan=None if canonical else (backbones or {}).get('red_parallel_layout')
    if red_plan and red_plan.get('enabled') and red_plan.get('verified'):
        planned={(r['rid'],r['mid']):r for r in red_plan['routes']}
        for edge,path in zip(edges,display_paths):
            route=planned.get((edge['reaction_id'],edge['metabolite_id']))
            if route:
                points=[tuple(p) for p in route['points']]
                if path['points'][0]!=points[0]:points.reverse()
                path.update(points=points,routing_locked=True)
    straight_paths=[] if canonical else [p for edge,p in zip(edges,display_paths) if frozenset((edge['source'],edge['target'])) not in tca_pairs]
    parallel_report=separate_parallel(straight_paths,obstacles=obstacles)
    parallel_report['node_detours']=avoid_nodes(straight_paths,obstacles=obstacles,route_offset=route_offset)
    parallel_report['routing_preference']=POLICY
    requested_red_spacing=red_plan['requested_spacing_grid'] if red_plan and red_plan.get('enabled') else .3
    parallel_report['important_simplification']=simplify_important(straight_paths,obstacles=obstacles,
        route_offset=route_offset,red_spacing=requested_red_spacing)
    parallel_report['gray_yield']=yield_gray(straight_paths,obstacles=obstacles,route_offset=route_offset)

    if red_plan:
        parallel_report['red_parallel_layout']={k:v for k,v in red_plan.items() if k not in ('routes','allocations')}
        parallel_report['red_parallel_spacing']=spacing_report(straight_paths,red_plan['requested_spacing_grid'])
        if red_plan.get('enabled') and red_plan.get('verified') and not parallel_report['red_parallel_spacing']['verified']:
            raise RuntimeError('Preallocated red spacing changed during SVG rendering')
    if canonical:
        saved={(p['rid'],p['mid']):p for p in canonical['routes']}
        for edge,path in zip(edges,display_paths):
            record=saved[edge['reaction_id'],edge['metabolite_id']]
            points=list(record['points'])
            if math.dist(points[0],path['points'][0])>1e-7:points.reverse()
            path['points']=points
        parallel_report['crossing_constraints']=canonical['report']
    for edge,path in zip(edges,display_paths):
        if frozenset((edge['source'],edge['target'])) in tca_pairs:
            a,b=path['points'][0],path['points'][-1]
            path['curve_controls']=arc_controls(a,b,tca_center)
            path['points']=sample_curve(a,path['curve_controls'],b)
    parallel_report['tca_curves']=dict(style='circular_bezier',curved_incidences=sum('curve_controls' in p for p in display_paths),
                                      scope='supported TCA cycle incidences only')
    route_points=[p for path in display_paths for p in path['points']]
    if route_points:
        ax.set_xticks(np.arange(math.floor(min(p[0] for p in route_points)),math.ceil(max(p[0] for p in route_points))+1))
        ax.set_yticks(np.arange(math.floor(min(p[1] for p in route_points)),math.ceil(max(p[1] for p in route_points))+1))
    parallel_report['node_obstacle_policy']='red_orange_structural_nodes'
    parallel_report['obstacle_node_ids']=obstacle_keys
    parallel_report['structural_node_priorities']=node_priorities
    parallel_report['original_route_length']=original_route_length
    parallel_report['display_route_length']=sum(math.dist(a,b) for p in display_paths for a,b in zip(p['points'],p['points'][1:]))
    dump(output/'parallel_routing.json',parallel_report)
    for j, edge in enumerate(edges):
        rid = edge['reaction_id']
        a, b = np.array(xy[edge['source']]), np.array(xy[edge['target']])
        reverse = edge['direction'] == 'reverse' if structural else stats[rid]['median'] < 0
        if reverse:
            a, b = b, a
        major = stats[rid].get('is_major_pathway', False)
        orange = j in orange_edges
        if (backbones or {}).get('adaptive_layout'):
            major=canonical['routes'][j]['priority']==2
            orange=canonical['routes'][j]['priority']==1
        color = RED if major else ORANGE if orange else '#cbd3d8' if connectors or structural else cmap(norm(stats[rid]['priority']))
        coenzyme = nodes[edge['source']]['kind'] == 'currency' or nodes[edge['target']]['kind'] == 'currency'
        points=np.array(display_paths[j]['points'])
        style=dict(color=color, lw=.7 if coenzyme else 2.3 if major or orange else 1.0,
                   alpha=.95 if major or orange else .35 if connectors else .65 if coenzyme else .55,
                   zorder=2 if major or orange else 1)
        if 'curve_controls' in display_paths[j]:
            from matplotlib.path import Path as PlotPath
            from matplotlib.patches import PathPatch
            vertices=[points[0],*display_paths[j]['curve_controls'],points[-1]]
            path=PathPatch(PlotPath(vertices,[PlotPath.MOVETO,*([PlotPath.CURVE4]*3)]),fill=False,**style)
            ax.add_patch(path)
        else:
            path, = ax.plot(points[:,0], points[:,1],**style)
        path.set_gid(f'edge-{j}')
        direction = points[-1]-points[-2]
        length = np.linalg.norm(direction)
        if length and (not structural or edge['direction'] != 'blocked'):
            start = b - direction / length * min(.25, length * .6)
            arrow = ax.annotate('', xy=b, xytext=start, arrowprops=dict(arrowstyle='->', color=color, lw=.8), zorder=2)
            arrow.set_gid(f'arrow-{j}')
            if structural:
                arrow.arrow_patch.set_gid(f'arrow-{j}')
        if structural and edge['direction'] == 'reversible':
            initial = points[1]-points[0]
            initial_length = np.linalg.norm(initial)
            if initial_length:
                start = a + initial / initial_length * min(.25, initial_length * .6)
                arrow = ax.annotate('', xy=a, xytext=start, arrowprops=dict(arrowstyle='->', color=color, lw=.8), zorder=2)
                arrow.set_gid(f'reverse-arrow-{j}')
                arrow.arrow_patch.set_gid(f'reverse-arrow-{j}')
    for j, (key, node) in enumerate(nodes.items()):
        p = xy[key]; rid = node['original_id'] if node['kind'] == 'reaction' else node['parent_id'].removeprefix('rxn:')
        is_reaction = node['kind'] == 'reaction'; coenzyme = node['kind'] == 'currency'
        major = stats[rid].get('is_major_pathway', False) if is_reaction else node['original_id'] in major_mids
        color = RED if node_priorities[key]==2 else ORANGE if node_priorities[key]==1 else '#83949d' if connectors or structural else cmap(norm(stats[rid]['priority'])) if is_reaction else '#83949d' if coenzyme else '#173c53'
        port_spec=((backbones or {}).get('crossing_routing') or {}).get('report',{})
        if key in port_spec.get('node_port_nodes',[]) and not port_spec.get('uniform_node_sizes'):
            from matplotlib.patches import Circle,Rectangle
            radius=port_spec['node_port_glyph_radius_grid']
            artist=Rectangle((p[0]-radius,p[1]-radius),2*radius,2*radius,facecolor=color,edgecolor='white',linewidth=.4,zorder=3) if is_reaction else Circle(p,radius,facecolor=color,edgecolor='white',linewidth=.4,zorder=3)
            ax.add_patch(artist)
        else:
            artist = ax.scatter(*p, marker='s' if is_reaction else 'o', s=25 if is_reaction else 9 if coenzyme else 32,
                                c=[color], zorder=3, edgecolor='white', linewidth=.4)
        artist.set_gid(f'node-{j}')
        label = node['original_id']
        if is_reaction and not structural:
            label += '\n' + f"{stats[rid]['median']:.3g}"
        text = ax.annotate(label, p, xytext=(4, 5) if not coenzyme else (3, -8), textcoords='offset points',
                           fontsize=5 if coenzyme else 6 if is_reaction else 7,
                           color=color, fontweight='normal' if coenzyme else 'bold',
                           zorder=22 if node_priorities[key]==2 else 21 if node_priorities[key]==1 else 4)
        if node_priorities[key]>0:
            text.set_path_effects([path_effects.Stroke(linewidth=2.5,foreground='white'),path_effects.Normal()])
        text.set_gid(f'label-{j}')
    provisional = bool(not structural and sampling_report and not sampling_report['mixing_checks_passed'])
    title_prefix = 'PROVISIONAL: sampling convergence criteria not met\n' if provisional else ''
    ax.set_title(title_prefix + f"{model_id} | {program_name} | top {selection['coverage_percent']:g}%: "
                 f"{len(selected)}/{selection['eligible_reactions']} reactions\n"
                 + ('No flux sampling; arrows show allowed bounds' if structural else
                  'Inorganic species hidden; integer-grid nodes; arrows follow signed median flux'), fontsize=13)
    ax.set_xlabel('Grid x'); ax.set_ylabel('Grid y')
    if priorities and not connectors:
        fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=.6, label='Other reactions: |median flux| (model flux units)')
    if not nodes:
        ax.text(.5, .5, 'No reactions selected', transform=ax.transAxes, ha='center')
    ax.set_autoscale_on(True)
    if frame:
        ax.update_datalim([(frame[0],frame[2]),(frame[1],frame[3])])
        ax.set_xticks(np.arange(math.floor(frame[0]),math.ceil(frame[1])+1))
        ax.set_yticks(np.arange(math.floor(frame[2]),math.ceil(frame[3])+1))
    ax.margins(.08)
    if (backbones or {}).get('adaptive_layout') and xy:
        # Measure actual node, arrow and label extents after the final aspect
        # ratio is applied. Shrink right/bottom independently of routing bounds.
        points=[p for route in display_paths for p in route['points']]+list(xy.values())
        xmin=min(p[0] for p in points);xmax=max(p[0] for p in points)
        ymin=min(p[1] for p in points);ymax=max(p[1] for p in points)
        ax.set_xlim(xmin-.5,xmax+.5);ax.set_ylim(ymin-.5,ymax+.5)
        for _ in range(3):
            fig.canvas.draw();renderer=fig.canvas.get_renderer();inverse=ax.transData.inverted()
            label_report=separate_static_labels(ax,renderer)
            boxes=[artist.get_window_extent(renderer) for artist in ax.get_children()
                   if (artist.get_gid() or '').startswith(('label-','node-','arrow-','reverse-arrow-'))]
            bounds=[inverse.transform(box.get_points()) for box in boxes if np.isfinite(box.get_points()).all()]
            right=max([xmax]+[b[1,0] for b in bounds]);bottom=min([ymin]+[b[0,1] for b in bounds])
            ax.set_xlim(xmin-.5,right+.25);ax.set_ylim(bottom-.25,ymax+.5)
        backbones['render_content_frame']=[xmin-.5,right+.25,bottom-.25,ymax+.5]
        fig.canvas.draw()
        label_report=separate_static_labels(ax,fig.canvas.get_renderer())
        dump(output/'render_validation.json',dict(frame=backbones['render_content_frame'],labels=label_report,uniform_node_sizes=True))
    oxidation_label='pyruvate oxidation' if program_name in ('sbml2map', 'sbml2escher') else 'PDH'
    connector_label = (f'Orange: shared undirected EX connectors; unique orange L1 length = {connectors["total_length"]:g}.\n' if connectors else '')
    if connectors and connectors.get('flux_costs', {}).get('flux_weight', 0):
        connector_label += (f'Flux-weighted orange cost = {connectors["total_cost"]:.4g}; '
                            f'low-flux penalty strength = {connectors["flux_costs"]["flux_weight"]:g}.\n')
    red_optimization=(backbones or {}).get('red_optimization')
    if red_optimization:
        connector_label=f'Red total routed length = {red_optimization["preservation"]["final_length_grid"]:.3f} grid units.\n'+connector_label
    fig.text(.02, .018, f'Red: glucose exchange / glycolysis / {oxidation_label} / TCA. Gray nodes: local coenzymes.\n' + connector_label +
             ('Red length is optimized first; changes to the frozen red core are recorded. Connectivity does not establish steady-state activity.' if structural and red_optimization else
              'Layout objective uses original L1 routes; display detours add length. Connectivity and bounds do not establish steady-state activity.' if structural else
              'Layout objective uses original L1 routes; display detours add length. Median describes the sampled feasible space, not a measured flux.'), fontsize=9)
    # Keep the original SVG for reversible browser highlighting.
    import io
    original_svg = io.StringIO() if any(highlights.values()) else None
    if original_svg is not None:
        fig.savefig(original_svg, format='svg', bbox_inches='tight')
    selected_gids = {gid for element in elements if element['key'] in highlights[element['kind']] for gid in element['ids']}
    for artist in ax.get_children():
        if artist.get_gid() not in selected_gids:
            continue
        if hasattr(artist, 'set_color'):
            artist.set_color(HIGHLIGHT_COLOR)
        if getattr(artist, 'arrow_patch', None) is not None:
            artist.arrow_patch.set_color(HIGHLIGHT_COLOR)
        if hasattr(artist, 'set_alpha'):
            artist.set_alpha(1)
    if any(highlights.values()):
        fig.text(.02, .005, 'Purple: user selection', color=HIGHLIGHT_COLOR, fontsize=9)
    fig.savefig(output / 'layout.png', dpi=180, bbox_inches='tight')
    fig.savefig(output / 'layout.svg', bbox_inches='tight')
    if connectors:
        visible = set()
        for j, edge in enumerate(edges):
            if j in orange_edges or stats[edge['reaction_id']].get('is_major_pathway',False):
                visible.update((f'edge-{j}', f'arrow-{j}'))
                visible.add(f'reverse-arrow-{j}')
        for j,(key,node) in enumerate(nodes.items()):
            if key in orange_nodes or (node['kind']=='reaction' and stats[node['original_id']].get('is_major_pathway',False)) or (node['kind']=='metabolite' and node['original_id'] in major_mids):
                visible.update((f'node-{j}',f'label-{j}'))
        for artist in ax.get_children():
            if artist.get_gid() and artist.get_gid() not in visible:
                artist.set_visible(False)
        fig.savefig(output/'orange_routes.png',dpi=180,bbox_inches='tight')
        fig.savefig(output/'orange_routes.svg',bbox_inches='tight')
    plt.close(fig)
    # Preserve SVG IDs for filtering without changing exported scientific coordinates.
    svg = original_svg.getvalue() if original_svg is not None else (output / 'layout.svg').read_text()
    svg = svg[svg.index('<svg'):]
    registry = []
    for j, edge in enumerate(edges):
        rank = display_ranks[edge['reaction_id']]
        is_co = nodes[edge['source']]['kind'] == 'currency' or nodes[edge['target']]['kind'] == 'currency'
        highlight = j in orange_edges or stats[edge['reaction_id']].get('is_major_pathway',False)
        registry.extend([dict(id=f'edge-{j}', rank=rank, coenzyme=is_co, highlight=highlight), dict(id=f'arrow-{j}', rank=rank, coenzyme=is_co, highlight=highlight)])
        if structural:
            registry.append(dict(id=f'reverse-arrow-{j}', rank=rank, coenzyme=is_co, highlight=highlight))
    for j, (key, node) in enumerate(nodes.items()):
        ranks = [display_ranks[e['reaction_id']] for e in edges if key in (e['source'], e['target'])]
        if node['kind'] == 'reaction':
            ranks = [display_ranks[node['original_id']]]
        rank = min(ranks) if ranks else 1  # Always-visible glucose reference, except at 0%.
        highlight = key in orange_nodes or (node['kind']=='reaction' and stats[node['original_id']].get('is_major_pathway',False)) or (node['kind']=='metabolite' and node['original_id'] in major_mids)
        registry.extend([dict(id=f'node-{j}', rank=rank, coenzyme=node['kind'] == 'currency',highlight=highlight),
                         dict(id=f'label-{j}', rank=rank, coenzyme=node['kind'] == 'currency',highlight=highlight,
                              label_kind='reaction' if node['kind']=='reaction' else 'metabolite',label_priority=node_priorities[key])])
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>sbml2map</title>
<style>body{font:15px system-ui;color:#203744;margin:20px}#viewer{height:78vh;overflow:auto;border:1px solid #ccd;background:white}#viewer svg{width:1600px;height:auto}header{margin-bottom:12px}a{margin-right:14px}</style>
<header><h1>sbml2map — __MODEL__</h1><p>Inorganic compounds are hidden. Red: main pathways (glucose exchange, glycolysis, PDH, TCA). Red pathways come first, followed by orange connectors, then other reactions. Each group is sorted by absolute median flux. Squares: reactions. Circles: compounds. Small gray circles: cofactors. Arrows follow the sign of median flux.</p>
<label>Top <input id="coverage" type="range" min="0" max="__MAX__" step="1" value="__MAX__"><span id="value"></span>%</label>
<label><input id="coenzymes" type="checkbox" checked>Cofactors</label>
<p id="grid-label-controls"><label><input id="show-metabolite-labels" type="checkbox">Metabolite labels</label>
<label><input id="show-reaction-labels" type="checkbox">Reaction labels</label>
<label><input id="allow-label-overlap" type="checkbox">Allow label overlap</label></p>
<label>Zoom <input type="range" min="600" max="4200" value="1600" step="100" oninput="document.querySelector('#viewer svg').style.width=this.value+'px';window.syncGridLabels?.()"></label>
<p>The slider filters the saved layout. To calculate a new layout, run the CLI with --coverage.</p>
<a href="layout.png">PNG</a><a href="layout.svg">SVG</a><a href="escher_preview.html">Escher</a><a href="report.json">Validation</a><a href="reaction_ranking.csv">Ranking</a></header>
<div id="viewer">__SVG__</div>
<script>const registry=__REGISTRY__,eligible=__ELIGIBLE__;
function update(){const v=Number(document.getElementById('coverage').value),limit=Math.ceil(v*eligible/100),co=document.getElementById('coenzymes').checked;document.getElementById('value').textContent=v;for(const item of registry){const el=document.getElementById(item.id);if(el)el.style.display=(item.rank<=limit&&(co||!item.coenzyme))?'':'none';}window.syncGridLabels?.();}
document.getElementById('coverage').oninput=update;document.getElementById('coenzymes').onchange=update;update();</script></html>'''
    import html
    if structural:
        page = page.replace('Red pathways come first, followed by orange connectors, then other reactions. Each group is sorted by absolute median flux. ',
                            'No flux sampling. Red pathways come first, followed by orange connectors, then other reactions. Each group is sorted by reaction ID, not activity. ')
        page = page.replace('Arrows follow the sign of median flux.', 'Arrows show directions allowed by model bounds. Reversible reactions have arrows both ways. Reactions with both bounds at zero have no arrows.')
    if provisional:
        notice=(f'<p class="sampling-warning" style="background:#fff0cc;padding:12px">Provisional result: sampling convergence checks failed. Rankings are provisional. Maximum split R-hat={sampling_report["max_split_rhat"]:.3f}, minimum ESS={sampling_report["min_effective_sample_size"]:.1f}.</p>')
        page=page.replace('</h1>', '</h1>'+notice,1)
    if program_name in ('sbml2map', 'sbml2escher'):
        page=page.replace(', PDH, ',', pyruvate oxidation, ')
    if connectors:
        page=page.replace('</h1>', '</h1><p>Orange: undirected paths connecting main pathways to exchange (EX) reactions. The layout minimizes total path length, counting shared branches once. Gray: other branches. Only EX reactions use the outer grid boundary. Filters may hide parts of these paths.</p>',1)
        page=page.replace('<label>Zoom','<label><input id="focus" type="checkbox" checked>Main pathways and EX connections only</label><label>Zoom',1)
        page=page.replace('(item.rank<=limit&&(co||!item.coenzyme))', '(item.rank<=limit&&(co||!item.coenzyme)&&(!document.getElementById("focus").checked||item.highlight))')
        page=page.replace('update();</script>', 'document.getElementById("focus").onchange=update;update();</script>')
        if structural:
            page=page.replace('id="focus" type="checkbox" checked', 'id="focus" type="checkbox"')
            missing = len(connectors.get('disconnected_exchange_node_ids', []))
            page=page.replace('</h1>', f'</h1><p>EX reactions with no connection to main pathways: {missing}. These stay gray. No connections are added.</p>',1)
    if connectors and connectors.get('flux_costs', {}).get('flux_weight', 0):
        page = page.replace('The layout minimizes total path length, counting shared branches once.',
                            'The layout minimizes length with a low-flux penalty, counting shared branches once.')
    if (backbones or {}).get('orange_reposition',{}).get('enabled'):
        page=page.replace('The layout minimizes total path length, counting shared branches once.', 'Local relocation reduces overlapping orange lines from different reactions, then compares total length; shared branches within one reaction count once. Point crossings are allowed.')
    page = page.replace('sbml2map', html.escape(program_name)).replace('__MODEL__', html.escape(model_id)).replace('__MAX__', str(selection['coverage_percent'])).replace('__SVG__', svg).replace('__REGISTRY__', json.dumps(registry)).replace('__ELIGIBLE__', str(selection['eligible_reactions']))
    page=page.replace('</html>','<script>window.gridLabelOptions='+json.dumps(label_options)+';</script><script>'+BASE.joinpath('grid_labels.js').read_text()+'</script></html>')
    exclusion=(backbones or {}).get('tca_exclusion',{})
    alignment=(backbones or {}).get('exchange_alignment',{})
    if exclusion.get('enabled'):
        interior_text='Interior node placement is allowed; the TCA perimeter remains fixed.' if exclusion.get('interior_placement_allowed') else 'No nodes are placed inside the TCA perimeter.'
        description=f"TCA shape: {exclusion['shape']}. {interior_text}"
        if exclusion['shape']=='circle':
            region=exclusion['region']
            if region.get('anchor_mode')=='grid_line_intersection':
                description+=' TCA nodes lie at circle/grid-line intersections; each has at least one integer coordinate. TCA cycle arcs remain curved.'
            else:
                description+=f" Circle through {region.get('lattice_point_count',0)} grid intersections; TCA nodes use a {region.get('node_grid_step') or 1:g}-step auxiliary grid and cycle connections use curved Beziers."

        if alignment.get('mode')=='collinear':
            description+=f" {alignment['checked_triplets']} matched EX / extracellular / cytosolic triplets are collinear, ordered from outside inward."
        page=page.replace('</header>','<p id="tca-layout-policy">'+description+'</p></header>',1)
    page = inject_controls(page, source_model, elements, highlights, records)
    (output / 'index.html').write_text(page)


SAMPLING_FILES = ('flux_samples.npz', 'flux_statistics.json', 'sampling_report.json', 'escher_flux_medians.json')


def validate_sampling_options(args):
    """Protect sampling outputs before any structural-mode files are written."""
    if not math.isfinite(getattr(args, 'orange_flux_weight', 1)) or getattr(args, 'orange_flux_weight', 1) < 0:
        raise ValueError('orange-flux-weight must be finite and nonnegative')
    if getattr(args, 'require_sampling', False) and getattr(args, 'no_sampling', False):
        raise ValueError('--require-sampling cannot be combined with --no-sampling')
    if not getattr(args, 'no_sampling', False):
        return
    if args.reuse_samples:
        raise ValueError('--no-sampling cannot be combined with --reuse-samples')
    if any((args.output_dir / name).exists() for name in SAMPLING_FILES):
        raise ValueError('Use a separate output directory for --no-sampling; sampling files already exist')


def sampling_data(args, model, digest):
    if args.reuse_samples:
        source = args.reuse_samples
        sampling_report = json.loads((source / 'sampling_report.json').read_text())
        if sampling_report.get('model_sha256') != digest:
            raise ValueError('Saved sampling model hash differs from the current model')
        sample_path = source / 'flux_samples.npz'
        if sampling_report.get('samples_sha256') != hashlib.sha256(sample_path.read_bytes()).hexdigest():
            raise ValueError('Saved sample checksum mismatch')
        with np.load(sample_path, allow_pickle=False) as saved:
            samples = saved['fluxes']
            if saved['reaction_ids'].tolist() != [r['id'] for r in model['reactions']]:
                raise ValueError('Saved sample reaction order mismatch')
        if samples.ndim != 3 or samples.shape[2] != len(model['reactions']) or not np.isfinite(samples).all():
            raise ValueError('Invalid saved sample dimensions/values')
        S, _, _, constraints = sampling_constraints(model)
        saved_constraints = sampling_report.get('constraints_sha256')
        if saved_constraints is not None and saved_constraints != constraints['constraints_sha256']:
            raise ValueError('Saved sampling constraint hash differs from the current constraints')
        if saved_constraints is None and constraints['boundary_metabolite_ids']:
            raise ValueError('Boundary-aware saved samples require a constraint fingerprint')
        sampling_report.update(constraints, performed=True, mode='reused')
        lo = np.array([r['lower_bound'] for r in model['reactions']]); hi = np.array([r['upper_bound'] for r in model['reactions']])
        if np.max(np.abs(samples @ S.T), initial=0) > 1e-7 or max(0, (lo - samples).max(), (samples - hi).max()) > 1e-7:
            raise ValueError('Saved samples violate model constraints')
        stats = json.loads((source / 'flux_statistics.json').read_text())
        medians = np.median(samples.reshape(-1, samples.shape[-1]), axis=0)
        if [s['reaction_id'] for s in stats] != [r['id'] for r in model['reactions']] or not np.allclose([s['median'] for s in stats], medians, rtol=0, atol=1e-12):
            raise ValueError('Saved statistics differ from samples')
    else:
        samples, stats, sampling_report = sample_fluxes(model, args.samples, args.chains, args.thinning,
                                                       args.burn_in, args.seed)
        sampling_report['model_sha256'] = digest
    return samples, stats, sampling_report


def run(args):
    from export_escher import create_map, validate_map, write_pathway_preview
    from escher_routing import separate_routes
    if getattr(args,'export_sbml_layout',False):
        from export_sbml_layout import validate_export_request
        validate_export_request(args.model,getattr(args,'sbml_layout_id','sbml2map_layout'))
    validate_sampling_options(args)
    route_offset=validate_route_offset(getattr(args,'route_offset',.3))
    red_spacing=validate_red_spacing(getattr(args,'red_parallel_spacing',1.0))
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    model = load_model(args.model)
    highlights = read_selection(args, model)
    strategy = getattr(args, '_strategy', None)
    program_name = getattr(args, 'program_name', 'sbml2map')
    digest = model_hash(model)
    no_sampling = getattr(args, 'no_sampling', False)
    if no_sampling:
        sampling_report = dict(performed=False, mode='disabled', model_sha256=digest)
    else:
        try:
            reason = getattr(args, '_sampling_unavailable_reason', None)
            if reason:
                raise ValueError(reason)
            samples, stats, sampling_report = sampling_data(args, model, digest)
        except (ValueError, RuntimeError, OSError, ImportError, MemoryError, ArithmeticError, KeyError, TypeError, np.linalg.LinAlgError) as exc:
            # Only sampling acquisition is recoverable; layout/export failures propagate.
            if getattr(args, 'require_sampling', False):
                raise RuntimeError('Required flux sampling failed: ' + str(exc)) from exc
            reason = f'{type(exc).__name__}: {exc}'
            print(f'Flux sampling unavailable ({reason}); switching to no sample mode.', flush=True)
            no_sampling = True
            # Preserve existing sampled maps and statistics.
            if any((output / name).exists() for name in SAMPLING_FILES):
                output = output / 'no_sampling'
            structural_args = argparse.Namespace(**vars(args))
            structural_args.output_dir = output
            structural_args.no_sampling = True
            structural_args.reuse_samples = None
            validate_sampling_options(structural_args)
            output.mkdir(parents=True, exist_ok=True)
            sampling_report = dict(performed=False, mode='fallback', reason=reason, model_sha256=digest)
    if not no_sampling:
        if args.reuse_samples is None or args.reuse_samples.resolve() != output.resolve():
            np.savez_compressed(output / 'flux_samples.npz', fluxes=samples,
                                reaction_ids=[r['id'] for r in model['reactions']])
        sampling_report['samples_sha256'] = hashlib.sha256((output / 'flux_samples.npz').read_bytes()).hexdigest()
        dump(output / 'sampling_report.json', sampling_report); dump(output / 'flux_statistics.json', stats)
        if not no_sampling and not sampling_report['mixing_checks_passed']:
            print('Sampling diagnostic thresholds not met; see sampling_report.json. Increase samples/thinning/burn-in.', flush=True)
    override = json.loads(args.classification.read_text()) if args.classification else None
    classification = strategy.classify(model, override) if strategy else classify_metabolites(model, override)
    routes, skipped_routes = strategy.pathways(model) if strategy else pathway_definitions(model)
    major_reactions = {rid for route in routes for rid in route['reaction_ids']}
    if no_sampling:
        selected, records, selection = select_structural_reactions(model, classification, args.coverage, major_reactions)
    else:
        selected, records, selection = select_reactions(model, stats, classification, args.coverage, args.zero_tolerance, major_reactions)
    metabolite_ids = {m['id'] for m in model['metabolites']}
    glucose_id = strategy.glucose_id if strategy else args.glucose_id or next((m for m in ('glc__D_e', 'glc__D_c') if m in metabolite_ids), None)
    if glucose_id is not None and glucose_id not in metabolite_ids:
        raise ValueError('Unknown --glucose-id: ' + glucose_id)
    if glucose_id is not None and classification[glucose_id]['kind'] != 'organic':
        raise ValueError('Glucose must be classified as organic for its fixed grid position')
    references = {glucose_id} if selected and glucose_id else set()
    shared = set()
    roots = {'met:'+mid for route in routes for mid in route['metabolite_ids']}
    roots |= {'rxn:'+rid for rid in major_reactions if rid in selected}
    if strategy and selected:
        visible_exchanges = {'rxn:'+e['reaction_id'] for e in strategy.extraction['exchanges'] if e['reaction_id'] in selected}
        shared = shared_coenzymes(model, selected, classification, roots, visible_exchanges, allow_disconnected=no_sampling)
    projected, nodes, edges = display_graph(model, selected, classification, references, shared)
    grid_ids = sorted(key for key, node in nodes.items() if node['kind'] != 'currency')
    grid_pairs = [(e['source'], e['target']) for e in edges
                  if nodes[e['source']]['kind'] != 'currency' and nodes[e['target']]['kind'] != 'currency']
    glucose_node = 'met:' + glucose_id if references else None
    exchange_nodes = strategy.exchanges(grid_ids) if strategy else {key for key in grid_ids if nodes[key]['kind'] == 'reaction' and nodes[key]['original_id'].startswith('EX_')}
    fixed_offsets, backbones = strategy.plan(routes, grid_ids) if strategy else backbone_plan(routes, grid_ids)
    backbones['skipped_routes'] = skipped_routes
    from tca_geometry import configure_plan
    fixed_offsets, backbones = configure_plan(fixed_offsets,backbones,getattr(args,'tca_shape','circle'),getattr(args,'tca_interior','allow'))
    adaptive=getattr(args,'crossing_policy','legacy')=='adaptive-grid'
    if adaptive:
        from adaptive_geometry import compact_tca
        preview_connectors=connector_plan(nodes,edges,roots,exchange_nodes,major_reactions,allow_disconnected=no_sampling) if strategy and selected else {}
        degrees={}
        for index in preview_connectors.get('edge_indices',[]):
            for key in (edges[index]['source'],edges[index]['target']):degrees[key]=degrees.get(key,0)+1
        required_ports={key for key,degree in degrees.items() if degree>1}
        fixed_offsets,backbones=compact_tca(fixed_offsets,backbones,grid_point_keys=required_ports)
    for route in routes:
        if route['id']=='tca':route['shape']=getattr(args,'tca_shape','circle')
    grid_options = strategy.grid_options(nodes, backbones) if strategy else {}
    # Seed a feasible exterior embedding before testing inward relocations.
    seed_region=backbones.get('tca_geometry')
    grid_options['reserved_region'] = dict(seed_region,interior_policy='forbid') if seed_region else None
    from exchange_alignment import discover_alignment
    alignment_mode = getattr(args,'exchange_alignment','collinear')
    alignment_groups, alignment_skipped = discover_alignment(model,nodes,exchange_nodes) if alignment_mode == 'collinear' else ([],[])
    grid_options['alignment_groups'] = alignment_groups
    flux_costs = orange_flux_costs(nodes, edges, major_reactions, records,
                                    0 if no_sampling else getattr(args, 'orange_flux_weight', 1))
    flux_costs['requested_flux_weight'] = getattr(args, 'orange_flux_weight', 1)
    reaction_weights = flux_costs['reaction_weights']
    connectors = connector_plan(nodes, edges, roots, exchange_nodes, major_reactions, allow_disconnected=no_sampling, reaction_weights=reaction_weights) if strategy and selected else None
    record_map={record['reaction_id']:record for record in records}
    backbone_pairs={frozenset(('rxn:'+step['reaction_id'],'met:'+mid))
        for route in backbones.get('routes',[]) for step in route.get('steps',[])
        for mid in (step['source'],step['target'])}
    from tca_curves import curved_pairs
    tca_pairs=curved_pairs(backbones)
    red_edges=[dict(rid=e['reaction_id'],mid=e['metabolite_id'],source=e['source'],target=e['target'],
        reverse=(e['direction']=='reverse' if no_sampling else record_map[e['reaction_id']]['median']<0),
        backbone=frozenset((e['source'],e['target'])) in backbone_pairs)
        for e in edges if e['reaction_id'] in major_reactions and e['source'] in grid_ids and e['target'] in grid_ids
        and frozenset((e['source'],e['target'])) not in tca_pairs]
    red_objective_pairs=[(e['source'],e['target']) for e in red_edges]
    if adaptive:red_edges=[]
    major_mids={mid for route in routes for mid in route['metabolite_ids']}
    def routing_obstacles(plan):
        orange_pairs={frozenset((edges[i]['source'],edges[i]['target'])) for i in (plan or {}).get('edge_indices',[])}
        return sorted(key for key,rank in structural_node_priorities(nodes,major_reactions,major_mids,orange_pairs).items() if rank and key in grid_ids)
    history = []
    red_planning_seconds=0.0
    best = None
    previous_xy = None
    for iteration in range(1 if adaptive else 2 if connectors else 1):
        objective_pairs = red_objective_pairs if adaptive else [(edges[i]['source'], edges[i]['target']) for i in connectors['edge_indices']] if connectors else grid_pairs
        objective_weights = None if adaptive else [reaction_weights.get(edges[i]['reaction_id'], 1) for i in connectors['edge_indices']] if connectors else None
        candidate_xy, candidate_placement = optimize_grid(grid_ids, objective_pairs, args.grid_cols, args.grid_rows, args.time_limit, args.seed + iteration,
                                      top_left_node=glucose_node, boundary_nodes=exchange_nodes, fixed_offsets=fixed_offsets,
                                      initial_positions=previous_xy, edge_weights=objective_weights, red_edges=red_edges,
                                      red_parallel_spacing=red_spacing,red_obstacle_node_ids=routing_obstacles(connectors), **grid_options)
        red_planning_seconds+=candidate_placement.get('red_parallel_layout',{}).get('elapsed_seconds',0)
        if not connectors:
            xy, placement = candidate_xy, candidate_placement
            break
        candidate_plan = connector_plan(nodes, edges, roots, exchange_nodes, major_reactions, candidate_xy, allow_disconnected=no_sampling, reaction_weights=reaction_weights)
        before_reroute = connector_cost(connectors, edges, candidate_xy)
        if connector_cost(candidate_plan, edges, candidate_xy) > before_reroute:
            candidate_plan = connectors
        value = connector_cost(candidate_plan, edges, candidate_xy)
        candidate_red=candidate_placement.get('red_parallel_layout')
        if red_edges and (candidate_red is None or candidate_red.get('obstacle_node_ids')!=routing_obstacles(candidate_plan)):
            candidate_red=plan_red_routes(candidate_xy,red_edges,red_spacing,routing_obstacles(candidate_plan),
                time_limit=min(3,max(.05,args.time_limit)),use_milp=not bool(getattr(args,'max_milp_nodes',None) and len(grid_ids)>args.max_milp_nodes))
            candidate_placement['red_parallel_layout']=candidate_red
            red_planning_seconds+=candidate_red.get('elapsed_seconds',0)
        red_valid=not red_spacing or not red_edges or bool(candidate_red and candidate_red.get('verified'))
        history.append(dict(iteration=iteration+1, red_spacing_verified=red_valid, before_reroute=before_reroute, orange_cost=value,
                            orange_length=connector_length(candidate_plan, edges, candidate_xy),
                            optimized_fixed_network_cost=candidate_placement['objective'],
                            optimized_fixed_network_length=connector_length(connectors, edges, candidate_xy)))
        if best is None or (not red_valid,value)<(not best[4],best[0]):
            best = (value, candidate_xy, candidate_placement, candidate_plan,red_valid)
        connectors = candidate_plan
        previous_xy = candidate_xy
    if connectors:
        value, xy, placement, connectors, _ = best
        connectors.update(total_length=connector_length(connectors, edges, xy), total_cost=value,
                          flux_costs=flux_costs, shared_coenzyme_ids=sorted(shared), optimization_history=history)
        backbones['exchange_connectors'] = connectors
        placement['fixed_network_solver_objective'] = placement['objective']
        if adaptive:placement['red_seed_nonarc_L1_objective']=placement['objective']
        placement['objective'] = value
        placement['objective_metric'] = connectors['objective']
        placement['optimality_proven'] = False
        # Fixed-forest MILP bounds cannot certify the joint routing problem.
        placement['fixed_network_lower_bound'] = placement.pop('global_lower_bound', None)
        placement.pop('relative_gap', None)
    red_plan=placement.pop('red_parallel_layout',None)
    if red_edges and (red_plan is None or red_plan.get('obstacle_node_ids')!=routing_obstacles(connectors)):
        red_plan=plan_red_routes(xy,red_edges,red_spacing,routing_obstacles(connectors),
            time_limit=min(3,max(.05,args.time_limit)),use_milp=not bool(getattr(args,'max_milp_nodes',None) and len(grid_ids)>args.max_milp_nodes))
        red_planning_seconds+=red_plan.get('elapsed_seconds',0)
    placement['red_parallel_layout_total_seconds']=red_planning_seconds
    if red_plan:
        backbones['red_parallel_layout']=red_plan
        placement['red_parallel_layout']={k:v for k,v in red_plan.items() if k not in ('routes','allocations')}
    cols, rows = placement['grid_cols'], placement['grid_rows']
    from tca_geometry import absolute_region,validate_exclusion
    requested_tca_region = absolute_region(backbones.get('tca_geometry'),rows)
    tca_region = dict(requested_tca_region,interior_policy='forbid') if requested_tca_region else None
    if args.grid_cols is not None and args.grid_rows is not None:
        backbones['render_frame']=[0,cols-1,0,rows-1]
    if getattr(args,'crossing_policy','legacy')=='red-hard' and alignment_mode=='collinear':
        from exchange_corridors import place_corridors
        from crossing_routing import initial_paths,route as route_red_first
        attach_coenzymes(nodes,xy,grid_shape=(cols,rows),forbidden_region=tca_region)
        backbones['tca_exclusion']={'region':tca_region}
        red_paths=[p for p in initial_paths(nodes,edges,xy,backbones,major_reactions) if p['priority']==2]
        red_report=route_red_first(red_paths,xy,boundary=[0,cols-1,0,rows-1],attempts=0)
        if not red_report['verified']:raise ValueError('Red-first routing failed: '+str(red_report))
        backbones['red_first']={'routes':red_paths,'report':red_report}
        backbones['_placement_region']=tca_region
        specification=place_corridors(model,nodes,xy,backbones,exchange_nodes,alignment_groups,cols,rows)
        backbones.pop('_placement_region',None)
        backbones['exchange_corridors']=specification
        backbones['exchange_boundary_adjustments']=[]
        placement['objective_before_straight_exchange_adjustment']=placement['objective']
        if connectors:placement['objective']=connector_cost(connectors,edges,xy)
        placement['global_optimality_proven']=False
        placement.pop('relative_gap',None)
    attach_coenzymes(nodes, xy, grid_shape=(cols, rows),forbidden_region=tca_region)
    if getattr(args,'crossing_policy','legacy')=='red-hard' and not backbones.get('exchange_corridors'):
        from exchange_alignment import anchor_exchange_boundary
        adjustments=anchor_exchange_boundary(alignment_groups,exchange_nodes,xy,cols,rows,tca_region)
        backbones['exchange_boundary_adjustments']=adjustments
        if adjustments:
            placement['objective_before_exchange_boundary_adjustment']=placement['objective']
            if connectors:placement['objective']=connector_cost(connectors,edges,xy)
            placement['global_optimality_proven']=False
            placement.pop('relative_gap',None)
    if (getattr(args,'compact_boundary',True) and getattr(args,'crossing_policy','legacy')=='red-hard'
            and backbones.get('exchange_corridors')):
        from boundary_compaction import compact_boundary
        compact_boundary(model,nodes,edges,xy,backbones,exchange_nodes,alignment_groups,cols,rows,
            major_reactions,attach_coenzymes,attempts=getattr(args,'crossing_attempts',12),margin=getattr(args,'boundary_margin',1.))
        if backbones['boundary_compaction']['accepted']:
            placement['objective_before_boundary_compaction']=placement['objective']
            if connectors:placement['objective']=connector_cost(connectors,edges,xy)
            placement['global_optimality_proven']=False
            placement.pop('relative_gap',None)
    tca_region=requested_tca_region
    backbones['tca_exclusion']={'region':tca_region}
    backbones.setdefault('routing_boundary',[0,cols-1,0,rows-1])
    if not adaptive and tca_region and tca_region.get('interior_policy')=='allow':
        from orange_reposition import place_gray_inside_tca
        backbones['tca_interior_placement']=place_gray_inside_tca(nodes,edges,xy,backbones)
    if (getattr(args,'reposition_orange',True) and getattr(args,'crossing_policy','legacy')=='red-hard'
            and backbones.get('exchange_corridors')):
        from orange_reposition import reposition_orange
        reposition_orange(nodes,edges,xy,backbones,major_reactions,alignment_groups,cols,rows,
            attempts=getattr(args,'crossing_attempts',12))
        if backbones['orange_reposition']['accepted_moves']:
            placement['objective_before_orange_reposition']=placement['objective']
            if connectors:placement['objective']=connector_cost(connectors,edges,xy)
            placement['global_optimality_proven']=False
            placement.pop('relative_gap',None)
    if adaptive:
        from adaptive_routing import adaptive_layout
        adaptive_layout(model,nodes,edges,xy,backbones,exchange_nodes,alignment_groups,cols,rows,major_reactions,
                        red_time_limit=getattr(args,'red_time_limit',10.),red_candidates=getattr(args,'red_candidates',10))
        tca_region=backbones['tca_exclusion']['region']
        if glucose_node:placement['position_constraints']['top_left_position']=list(xy[glucose_node])
        if connectors:placement['objective']=connector_cost(connectors,edges,xy)
        placement['global_optimality_proven']=False
        if backbones.get('red_optimization'):
            placement['orange_objective']=placement['objective']
            placement['objective']=backbones['red_optimization']['preservation']['final_length_grid']
            placement['objective_metric']=backbones['red_optimization']['objective']
            placement['objective_order']=backbones['red_optimization']['objective_order']
    routing_boundary=backbones.get('routing_boundary',[0,cols-1,0,rows-1])
    xmin,xmax,ymin,ymax=routing_boundary
    tca_exclusion = validate_exclusion(nodes,xy,tca_region)
    backbones['tca_exclusion'] = tca_exclusion
    from exchange_alignment import validate_alignment
    exchange_alignment = validate_alignment(alignment_groups,xy,cols,rows,routing_boundary)
    exchange_alignment.update(mode=alignment_mode,skipped=alignment_skipped)
    backbones['exchange_alignment'] = exchange_alignment
    backbone_validation = validate_backbones(backbones, xy, rows)
    backbones['validation'] = backbone_validation
    if glucose_node and xy[glucose_node] != placement['position_constraints']['top_left_position']:
        raise RuntimeError('Glucose is not at its required top-left position')
    if any(xy[key][0] not in (xmin, xmax) and xy[key][1] not in (ymin, ymax) for key in exchange_nodes):
        raise RuntimeError('An EX reaction is not on the grid perimeter')
    if strategy and any(xy[key][0] in (xmin, xmax) or xy[key][1] in (ymin, ymax) for key in set(grid_ids)-exchange_nodes):
        raise RuntimeError('A non-EX grid node occupies the perimeter')
    if any(not (xmin <= p[0] <= xmax and ymin <= p[1] <= ymax) for p in xy.values()):
        raise RuntimeError('A displayed node lies outside the grid rectangle')
    network = validate_network(projected, nodes, edges, xy)
    fractional_allowed = set(backbones.get('exchange_corridors',{}).get('straight_node_ids',())) | set((tca_region or {}).get('node_ids',())) | {g[k] for g in alignment_groups for k in ('exchange','external','internal')}
    for key in grid_ids:
        if key not in fractional_allowed and not np.allclose(xy[key], np.rint(xy[key]), atol=1e-9, rtol=0):
            raise RuntimeError('A non-coenzyme node is off-grid')
    if len({tuple(xy[k]) for k in grid_ids}) != len(grid_ids):
        raise RuntimeError('Grid nodes overlap')
    if any(hidden_metabolite(classification[n['original_id']]) for n in nodes.values() if n['kind'] != 'reaction'):
        raise RuntimeError('Inorganic species in display')
    length = float(sum(np.abs(np.array(xy[e['source']]) - xy[e['target']]).sum() for e in edges))
    coenzyme_count = sum(n['kind'] == 'currency' for n in nodes.values())
    grid_length = sum(sum(abs(a-b) for a,b in zip(xy[a_key],xy[b_key])) for a_key,b_key in grid_pairs)
    if not adaptive and not math.isclose(length, grid_length + .45 * coenzyme_count, abs_tol=1e-7):
        raise RuntimeError('Edge-length objective differs from drawn orthogonal geometry')
    if connectors and not math.isclose(connector_cost(connectors,edges,xy),placement.get('orange_objective',placement['objective']),abs_tol=1e-7):
        raise RuntimeError('Orange objective differs from drawn orthogonal geometry')
    original = {r['id']: r for r in model['reactions']}
    for reaction in projected['reactions']:
        expected = {m: c for m, c in original[reaction['id']]['metabolites'].items() if not hidden_metabolite(classification[m])}
        if reaction['metabolites'] != expected:
            raise RuntimeError('Visible stoichiometry differs from source')
    if getattr(args,'crossing_policy','legacy')=='red-hard':
        backbones.setdefault('routing_boundary',[0,cols-1,0,rows-1])
        from crossing_routing import plan_routes
        if 'crossing_routing' not in backbones:
            backbones['crossing_routing']=plan_routes(nodes,edges,xy,backbones,major_reactions,
                attempts=getattr(args,'crossing_attempts',24))
    metadata = dict(source_model_sha256=digest, coverage=selection, sampling=sampling_report, major_backbones=backbones,
                    original_reactions={r['id']: r for r in model['reactions']}, classification=classification,
                    note='Display projection only: inorganic participants omitted; source_model.json retains the original source model.')
    if nodes:
        data = create_map(projected, nodes, edges, xy, scale=180, padding=80,
                          layout_method=program_name + (' structural grid map' if no_sampling else ' flux-median-ranked grid map'))
        data[0]['metadata'][program_name] = metadata
        data[0]['map_description'] = ('Structural display projection. No sampling or flux statistics; directions from source bounds.' if no_sampling else
                                      'Flux median ranked display projection. Inorganic participants hidden. Sampling uses the full original model.')
        orange_pairs = ({frozenset((edges[i]['source'],edges[i]['target']))
                         for i in connectors['edge_indices']} if connectors else None)
        separate_routes(data,major_reactions,orange_pairs,route_offset=route_offset)
        escher_check = validate_map(data, projected, nodes, xy)
        reaction_data = None if no_sampling else {r['reaction_id']: r['median'] * data[0]['metadata']['reactions'][r['reaction_id']]['coefficient_multiplier']
                         for r in records if r['selected']}
        orange_pairs = ({frozenset((edges[i]['source'],edges[i]['target']))
                         for i in connectors['edge_indices']} if connectors else None)
        data[0]['metadata']['exchange_reaction_ids'] = sorted(nodes[k]['original_id'] for k in exchange_nodes)
        data[0]['metadata']['label_display'] = dict(
            show_metabolites=not getattr(args, 'hide_metabolite_labels', False),
            show_reactions=not getattr(args, 'hide_reaction_labels', False),
            allow_overlap=getattr(args, 'allow_label_overlap', False))
        write_pathway_preview(output/'escher_preview.html',data,major_reactions,orange_pairs,records=records)

    else:
        data = [dict(schema='https://escher.github.io/escher/jsonschema/1-0-0#', homepage='https://escher.github.io',
                     map_id=digest[:24], map_name='sbml2map empty', map_description='No reactions selected'),
                dict(nodes={}, reactions={}, text_labels={}, canvas=dict(x=0, y=0, width=500, height=500))]
        from jsonschema import Draft4Validator
        from export_escher import SCHEMA
        Draft4Validator(json.loads(SCHEMA.read_text())).validate(data)
        escher_check = dict(schema_valid=True, reactions=0)
        reaction_data = {}
        (output / 'escher_preview.html').write_text('<!doctype html><meta charset="utf-8"><p>No reactions selected</p>')
    dump(output / 'escher_map.json', data)
    if not no_sampling:
        dump(output / 'escher_flux_medians.json', reaction_data)
    dump(output / 'source_model.json', model); dump(output / 'display_model.json', projected)
    dump(output / 'metabolite_classification.json', classification); dump(output / 'reaction_ranking.json', records)
    dump(output / 'escher_map_metadata.json', metadata)
    dump(output / 'major_backbones.json', backbones)
    if backbones.get('red_optimization'):dump(output/'red_optimization.json',backbones['red_optimization'])
    with (output / 'reaction_ranking.csv').open('w', newline='') as stream:
        fields = (['reaction_id', 'rank', 'is_major_pathway', 'selected', 'exclusion_reason', 'lower_bound', 'upper_bound'] if no_sampling else
                  ['reaction_id', 'rank', 'is_major_pathway', 'median_rank', 'median', 'priority', 'selected', 'exclusion_reason', 'split_rhat', 'ess', 'q05', 'q95'])
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore'); writer.writeheader(); writer.writerows(records)
    with (output / 'coordinates.csv').open('w', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['id', 'kind', 'original_id', 'x', 'y', 'on_grid'])
        for key, node in nodes.items():
            writer.writerow([key, node['kind'], node['original_id'], *xy[key], bool(np.allclose(xy[key],np.rint(xy[key]),atol=1e-9,rtol=0))])
    with (output / 'edges.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['source', 'target', 'reaction_id', 'metabolite_id', 'coefficient', 'direction'])
        writer.writeheader(); writer.writerows(edges)
    report = dict(source_model=str(args.model.resolve()), source_model_sha256=digest,
                  selection=selection, sampling=sampling_report, placement=placement, major_backbones=backbones,
                  source_reactions=len(model['reactions']), drawn_reactions=len(selected),
                  grid_nodes=len(grid_ids), coenzyme_copies=coenzyme_count, edges=len(edges),
                  glucose_metabolite_id=glucose_id, glucose_anchor_active=bool(glucose_node),
                  glucose_reference_only=bool(glucose_node and not any(e['metabolite_id'] == glucose_id for e in edges)),
                  exchange_reaction_ids=sorted(nodes[k]['original_id'] for k in exchange_nodes),
                  total_orthogonal_edge_length=length, coenzyme_edge_length=.45 * coenzyme_count,
                  inorganic_metabolites=[m for m, c in classification.items() if c['kind'] == 'inorganic'],
                  validation=dict(all_noncoenzyme_nodes_on_distinct_grid_points=all(np.allclose(xy[k],np.rint(xy[k]),atol=1e-9,rtol=0) for k in grid_ids),
                                  all_non_tca_nonaligned_nodes_on_grid=True,tca_interior=tca_exclusion,exchange_alignment=exchange_alignment,no_inorganic_nodes=True,
                                  glucose_top_left_verified=True, exchange_reactions_on_perimeter=True,
                                  all_nodes_within_grid_rectangle=True, major_backbones=backbone_validation,
                                  no_zero_median_reactions=True, visible_stoichiometry_preserved=True,
                                  full_model_used_for_sampling=True, network=network, escher=escher_check),
                  limitations=['Global length optimality is reported, not assumed.',
                               ('Orange overlaps outside common node ports and close parallel runs of at least two grids are forbidden. Red/orange point crossings require an audited reason.' if backbones.get('adaptive_layout') else 'Red crossings are forbidden; orange length is minimized over a bounded candidate search. Orange crossings are permitted.' if backbones.get('crossing_routing') else 'Crossings, coincident edge paths, and label overlap are not minimized.'),
                               'Display model omits inorganic participants and is not suitable for FBA.',
                               'TCA cycle connections may use circular Beziers; other connections use orthogonal routing. Biological coordinates and reaction chemistry are preserved.',
                               'Sampling uses original bounds without a growth floor or loopless constraints; internal cycles may dominate.'])
    if no_sampling:
        report['validation'].pop('no_zero_median_reactions')
        report['validation'].pop('full_model_used_for_sampling')
        report['validation']['no_sampling_or_flux_statistics_used'] = True
        report['validation']['arrows_follow_source_bounds'] = True
        report['limitations'][-1] = 'Structural connectivity includes inactive and zero-bound reactions; steady-state feasibility is not evaluated.'
    placeholders = [mid for mid, c in classification.items() if c.get('boundary_placeholder')]
    if placeholders:
        report['boundary_placeholder_metabolites'] = placeholders
        report['validation']['no_boundary_placeholder_nodes'] = True
        report['limitations'].append('CellDesigner source/sink placeholders are omitted only in the display; original stoichiometry and sampling constraints are unchanged.')
    report['red_parallel_spacing_grid'] = red_spacing
    report['route_offset_grid'] = route_offset
    report['program_name'] = program_name
    report['routing_boundary']=routing_boundary
    if backbones.get('boundary_compaction'):report['boundary_compaction']=backbones['boundary_compaction']
    if backbones.get('orange_reposition'):report['orange_reposition']=backbones['orange_reposition']
    if backbones.get('adaptive_layout'):report['adaptive_layout']=backbones['adaptive_layout']
    if backbones.get('red_optimization'):report['red_optimization']=backbones['red_optimization']
    if backbones.get('tca_interior_placement'):report['tca_interior_placement']=backbones['tca_interior_placement']
    report['tca_shape'] = getattr(args,'tca_shape','circle')
    report['tca_exclusion'] = tca_exclusion
    report['exchange_alignment'] = exchange_alignment
    if strategy:
        report['extraction_status'] = strategy.extraction['status']
        report['validation']['glucose_is_top_left_compound'] = True
        report['validation']['perimeter_exclusive_to_exchanges'] = True
        if connectors:
            report['exchange_connectors'] = connectors
            report['orange_orthogonal_edge_length'] = connectors['total_length']
            report['orange_weighted_cost'] = connectors['total_cost']
            report['validation']['all_exchanges_connected_to_major_pathways'] = connectors['all_exchanges_connected']
            dump(output/'exchange_connectors.json',connectors)
        strategy.save(output)
    elements = grid_elements(nodes, edges)
    highlight_catalog = catalog(model, elements, records)
    report['highlights'] = dict(selection=highlights, unavailable={kind: [item for item in highlight_catalog[kind]
                               if item['id'] in highlights[kind] and not item['drawn']] for kind in highlights})
    if any(report['highlights']['unavailable'].values()):
        print('Highlight targets not in saved map: ' + json.dumps(report['highlights']['unavailable']), flush=True)
    preview = output / 'escher_preview.html'
    preview.write_text(inject_controls(preview.read_text(), model, elements, highlights, records, escher=True))
    data[0].setdefault('metadata', {})['highlights'] = highlights
    dump(output / 'escher_map.json', data)
    dump(output / 'report.json', report)
    render(output, nodes, edges, xy, records, selection, model['id'], backbones, program_name, sampling_report, highlights, model, route_offset=route_offset,label_options=data[0].get('metadata',{}).get('label_display'))
    if getattr(args,'export_sbml_layout',False):
        from export_sbml_layout import export_sbml_layout
        layout_report=export_sbml_layout(args.model,output,layout_id=getattr(args,'sbml_layout_id','sbml2map_layout'))
        report['sbml_layout']=dict(file='layout.sbml',report='sbml_layout_report.json',layout_id=layout_report['layout_id'],
                                   validation=layout_report['validation'])
        dump(output/'report.json',report)
    print(json.dumps(dict(output=str(output.resolve()), selection=selection, sampling_passed=sampling_report.get('mixing_checks_passed'),
                          placement=placement, total_orthogonal_edge_length=length), indent=2), flush=True)
    return report


def validate_red_spacing_arg(value):
    try:return validate_red_spacing(float(value))
    except ValueError as exc:raise argparse.ArgumentTypeError(str(exc)) from exc


def validate_route_offset_arg(value):
    try:
        return validate_route_offset(float(value))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hide-metabolite-labels', action='store_true', help='Initially hide metabolite labels in interactive SVG and Escher views')
    p.add_argument('--hide-reaction-labels', action='store_true', help='Initially hide reaction labels in interactive SVG and Escher views')
    p.add_argument('--allow-label-overlap', action='store_true', help='Show all enabled Escher labels even when they overlap')
    p.add_argument('--highlight-metabolites', nargs='+', action='extend', default=[], help='Exact source metabolite IDs to highlight blue')
    p.add_argument('--highlight-reactions', nargs='+', action='extend', default=[], help='Exact source reaction IDs to highlight blue')
    p.add_argument('--highlight-file', type=Path, help='JSON with metabolites and reactions ID lists')
    p.add_argument('--export-sbml-layout', action='store_true', help='Write positions and connecting lines to layout.sbml (SBML Level 3 Layout; off by default)')
    p.add_argument('--sbml-layout-id', default='sbml2map_layout', help='ID for the added SBML Layout (existing layouts are preserved)')
    p.add_argument('--model', type=Path, default=BASE / 'e_coli_core.xml')
    p.add_argument('--output-dir', type=Path, default=BASE / 'e_coli_core_sbml2map_output')
    p.add_argument('--coverage', type=float, default=100, help='Top percentage of eligible reactions (0..100); --no-sampling uses pathway membership and IDs')
    p.add_argument('--crossing-policy',choices=('adaptive-grid','red-hard','legacy'),default='legacy',help='Forbid red/red and red/orange crossings; minimize orange length')
    p.add_argument('--routing-attempts','--crossing-attempts',dest='crossing_attempts',type=int,default=12,help='Orange length ordering candidates (default: 12)')
    p.add_argument('--no-reposition-orange',dest='reposition_orange',action='store_false',help='Disable local node relocation to reduce orange line overlap')
    p.add_argument('--no-compact-boundary',dest='compact_boundary',action='store_false',help='Keep the initial lower and right exchange boundaries')
    p.add_argument('--boundary-margin',type=float,default=1.,help='Minimum clearance around frozen red geometry during boundary compaction (default: 1 grid unit)')
    p.add_argument('--tca-interior',choices=('allow','forbid'),default='allow',help='Allow ordinary nodes and cofactor copies inside the fixed TCA perimeter (default: allow)')
    p.add_argument('--tca-shape',choices=('circle','square'),default='circle',help='TCA shape: circle (default) or square; only the TCA cycle follows clockwise order')
    p.add_argument('--exchange-alignment',choices=('collinear','off'),default='collinear',help='Align EX, extracellular and matching cytosolic nodes from outside inward (default: collinear)')
    p.add_argument('--red-parallel-spacing',type=validate_red_spacing_arg,default=1.0,help='Minimum red/red channel spacing for parallel runs of at least one grid unit (default: 1; 0 disables layout channel allocation)')
    p.add_argument('--route-offset', type=validate_route_offset_arg, default=.3, help='Whole-segment U-shaped node avoidance offset in grid units (default: 0.3; 0 disables node avoidance)')
    p.add_argument('--zero-tolerance', type=float, default=1e-7)
    p.add_argument('--orange-flux-weight', type=float, default=1, help='Penalty for low-flux orange routes (default: 1; 0: length only; larger values favor higher-flux reactions more; ignored without sampling)')
    p.add_argument('--require-sampling', action='store_true', help='Fail instead of falling back to a structural drawing if sampling is unavailable')
    p.add_argument('--samples', type=int, default=2000, help='Retained samples per chain')
    p.add_argument('--chains', type=int, default=4)
    p.add_argument('--thinning', type=int, default=100)
    p.add_argument('--burn-in', type=int, default=10000)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--time-limit', type=float, default=60, help='Grid MILP time limit in seconds')
    p.add_argument('--glucose-id', help='Glucose species to anchor; default glc__D_e, then glc__D_c if present')
    p.add_argument('--grid-cols', type=int)
    p.add_argument('--grid-rows', type=int)
    p.add_argument('--classification', type=Path, help='JSON metabolite ID -> inorganic/coenzyme/organic')
    sampling = p.add_mutually_exclusive_group()
    sampling.add_argument('--reuse-samples', type=Path, help='Reuse a validated sampling output directory; reoptimize layout')
    sampling.add_argument('--no-sampling', action='store_true', help='Draw source structure with no flux sampling/statistics; rank major pathways first, then IDs')
    return p
