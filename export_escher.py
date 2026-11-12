#!/usr/bin/env python3
"""Export a validated reaction-node layout as an Escher 1-0-0 map."""
import argparse
import base64
from collections import Counter
import csv
import copy
import hashlib
import html
import json
import math
from pathlib import Path
import re

from reaction_network import build_network, reaction_direction, validate_network
from escher_routing import separate_routes
from red_route_layout import plan_from_metadata
from routing_simplification import POLICY

BASE = Path(__file__).resolve().parent
SCHEMA = BASE / 'vendor/escher/map-schema-1-0-0.json'
VIEWER = BASE / 'vendor/escher/escher-1.8.1.min.js'


def fba_flux_style(value, maximum, first):
    """Square-root intensity on absolute flux; keep signed reaction data intact."""
    if value == 0:
        return dict(type='value', value=0, color='#b8b8b8', size=1)
    low = math.sqrt(first / maximum)
    intensity = max(0.0, min(1.0, (math.sqrt(abs(value) / maximum) - low) / (1 - low)))
    anchors = [(0.0, '#2166ac'), (.35, '#756bb1'), (.7, '#c05a88'), (1.0, '#d73027')]
    for (a, left), (b, right) in zip(anchors, anchors[1:]):
        if intensity <= b:
            weight = (intensity - a) / (b - a)
            channels = [round(int(left[i:i+2], 16) * (1-weight) + int(right[i:i+2], 16) * weight)
                        for i in (1, 3, 5)]
            return dict(type='value', value=value,
                        color='#' + ''.join(f'{channel:02x}' for channel in channels),
                        size=3 + 5 * intensity)


def fba_flux_scale(maximum, first):
    # Escher interpolates linearly in raw flux. Dense stops approximate the
    # square-root transfer without changing data, labels, tooltips or arrows.
    low = math.sqrt(first / maximum)
    values = [first] + [maximum * (low + (1-low)*i/64)**2 for i in range(1, 64)] + [maximum]
    return [fba_flux_style(0, maximum, first)] + [fba_flux_style(v, maximum, first) for v in values]


def read_layout(directory):
    model = json.loads((directory / 'model.json').read_text())
    with (directory / 'coordinates.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    if len({row['id'] for row in rows}) != len(rows):
        raise ValueError('Duplicate coordinate IDs')
    currency = sorted({row['original_id'] for row in rows if row['kind'] == 'currency'})
    nodes, expected_edges = build_network(model, currency)
    if {row['id'] for row in rows} != set(nodes):
        raise ValueError('Coordinates do not match the reaction graph')
    xy = {}
    for row in rows:
        node = nodes[row['id']]
        if any(row[field] != node[field] for field in ('kind', 'original_id', 'parent_id')):
            raise ValueError('Coordinate node identity mismatch')
        xy[row['id']] = [float(row['x']), float(row['y'])]
    with (directory / 'edges.csv').open() as stream:
        edges = list(csv.DictReader(stream))
    for edge in edges:
        edge['coefficient'] = float(edge['coefficient'])
    validate_network(model, nodes, edges, xy)
    # Validate input CSV reaction directions against actual model bounds as well.
    if len(edges) != len(expected_edges):
        raise ValueError('Missing reaction incidences')
    return model, nodes, edges, xy


def create_map(model, nodes, edges, xy, scale=4000.0, padding=180.0, layout_method='MILP reaction layout'):
    if not math.isfinite(scale) or scale <= 0 or not math.isfinite(padding) or padding < 0:
        raise ValueError('Scale must be positive and padding nonnegative')
    validate_network(model, nodes, edges, xy)
    x_min = min(p[0] for p in xy.values())
    y_max = max(p[1] for p in xy.values())
    def transform(p):
        return scale * (p[0] - x_min) + padding, scale * (y_max - p[1]) + padding
    ids = {key: str(i) for i, key in enumerate(sorted(nodes))}
    metabolites = {m['id']: m for m in model['metabolites']}
    map_nodes = {}
    for key in sorted(nodes):
        node = nodes[key]
        x, y = transform(xy[key])
        if node['kind'] == 'reaction':
            record = dict(node_type='midmarker', x=x, y=y)
        else:
            mid = node['original_id']
            record = dict(node_type='metabolite', x=x, y=y, bigg_id=mid,
                          name=metabolites[mid].get('name', mid), node_is_primary=node['kind'] == 'metabolite',
                          label_x=x, label_y=y - 24)
            if node['kind'] == 'currency':
                parent = xy[node['parent_id']]
                record.update(label_x=x + (10 if xy[key][0] >= parent[0] else -55),
                              label_y=y + (14 if xy[key][1] <= parent[1] else -10))
        map_nodes[ids[key]] = record
    incidence = {}
    for edge in edges:
        incidence.setdefault(edge['reaction_id'], []).append(edge)
    reactions, reaction_metadata, text_labels = {}, {}, {}
    gene_names = {g['id']: g.get('name', '') for g in model.get('genes', [])}
    segment_number = len(nodes) + len(model['reactions'])
    missing_gene_rules = []
    routing_node_ids = []
    segment_source_edges = {}
    for i, reaction in enumerate(model['reactions']):
        rid = reaction['id']
        reaction_id = str(len(nodes) + i)
        marker = ids['rxn:' + rid]
        x, y = transform(xy['rxn:' + rid])
        direction = reaction_direction(reaction['lower_bound'], reaction['upper_bound'])
        multiplier = -1 if direction == 'reverse' else 1
        rule = reaction.get('gene_reaction_rule', '')
        if 'gene_reaction_rule' not in reaction:
            missing_gene_rules.append(rid)
        genes = sorted(set(re.findall(r'[A-Za-z0-9_.:-]+', rule)) - {'and', 'or'})
        record = dict(name=reaction.get('name', rid), bigg_id=rid, reversibility=direction == 'reversible',
                      label_x=x + 18, label_y=y - 15, gene_reaction_rule=rule,
                      genes=[dict(bigg_id=g, name=gene_names.get(g, '')) for g in genes],
                      metabolites=[dict(bigg_id=mid, coefficient=multiplier * coefficient)
                                   for mid, coefficient in sorted(reaction['metabolites'].items())], segments={})
        for edge in incidence.get(rid, []):
            key = edge['source'] if nodes[edge['source']]['kind'] != 'reaction' else edge['target']
            target = ids[key]
            a, b = map_nodes[marker], map_nodes[target]
            path = [marker]
            if not math.isclose(a['x'], b['x'], abs_tol=1e-9, rel_tol=0) and not math.isclose(a['y'], b['y'], abs_tol=1e-9, rel_tol=0):
                # Substrates enter horizontally; products leave vertically after
                # an elbow. The terminal metabolite segment owns the arrowhead.
                substrate = multiplier * edge['coefficient'] < 0
                bend = str(len(map_nodes))
                map_nodes[bend] = dict(node_type='multimarker',
                                       x=a['x'] if substrate else b['x'],
                                       y=b['y'] if substrate else a['y'])
                routing_node_ids.append(bend)
                path.append(bend)
            path.append(target)
            for start, end in zip(path, path[1:]):
                sid = str(segment_number)
                record['segments'][sid] = dict(from_node_id=start, to_node_id=end, b1=None, b2=None)
                segment_source_edges[sid] = [edge['source'], edge['target']]
                segment_number += 1
        reactions[reaction_id] = record
        def bound(value):
            return value if math.isfinite(value) else str(value)
        reaction_metadata[rid] = dict(map_reaction_id=reaction_id, marker_node_id=marker,
                                     original_lower_bound=bound(reaction['lower_bound']),
                                     original_upper_bound=bound(reaction['upper_bound']),
                                     original_stoichiometry=reaction['metabolites'],
                                     coefficient_multiplier=multiplier, direction=direction)
        if direction == 'blocked':
            text_labels[str(len(text_labels))] = dict(x=x + 18, y=y + 35, text=rid + ': blocked (lb=ub=0)')
    coords = [(n['x'], n['y']) for n in map_nodes.values()]
    coords += [(n['label_x'], n['label_y']) for n in map_nodes.values() if n['node_type'] == 'metabolite']
    coords += [(r['label_x'], r['label_y']) for r in reactions.values()]
    coords += [(t['x'], t['y']) for t in text_labels.values()]
    lo_x, lo_y = min(p[0] for p in coords) - padding, min(p[1] for p in coords) - padding
    hi_x, hi_y = max(p[0] for p in coords) + padding, max(p[1] for p in coords) + padding
    # Extra right margin provides room for long reaction IDs and currency labels.
    canvas = dict(x=lo_x, y=lo_y, width=hi_x - lo_x + 350, height=hi_y - lo_y)
    metadata = dict(exporter='drawing/export_escher.py',
                    source_node_ids={ids[key]: key for key in ids},
                    routing='orthogonal', routing_node_ids=routing_node_ids, segment_source_edges=segment_source_edges,
                    coordinate_transform=dict(scale=scale, padding=padding, x_min=x_min, y_max=y_max, y_inverted=True),
                    reactions=reaction_metadata, currency_ids=sorted({n['original_id'] for n in nodes.values() if n['kind'] == 'currency'}),
                    missing_gene_rule_reactions=missing_gene_rules,
                    limitation='Bounds and blocked status are metadata; reverse-only coefficients are negated. Multiply corresponding fluxes by coefficient_multiplier.')
    signature = hashlib.sha256(json.dumps([model, xy, scale, padding, 'orthogonal-v1'], sort_keys=True).encode()).hexdigest()[:24]
    header = dict(schema='https://escher.github.io/escher/jsonschema/1-0-0#', homepage='https://escher.github.io',
                  map_id=signature, map_name=model['id'] + '.' + layout_method,
                  map_description=layout_method + '. Reaction-local currency copies, original bounds and node identities are retained in metadata.',
                  metadata=metadata)
    return [header, dict(nodes=map_nodes, reactions=reactions, text_labels=text_labels, canvas=canvas)]


def validate_map(data, model, source_nodes, source_xy, schema_path=SCHEMA, source_metadata=None):
    try:
        from jsonschema import Draft4Validator
    except ImportError as exc:
        raise RuntimeError('jsonschema is required: use a Python with jsonschema or install requirements-escher.txt') from exc
    Draft4Validator(json.loads(schema_path.read_text())).validate(data)
    report = validate_map_geometry(data, model, source_nodes, source_xy, source_metadata)
    report['schema_valid'] = True
    return report


def validate_map_geometry(data, model, source_nodes, source_xy, source_metadata=None):
    """Verify identity, chemistry and connected orthogonal paths independently of schema."""
    if len(data) != 2:
        raise ValueError('An Escher map must have exactly two items')
    header, body = data
    map_nodes, reactions = body['nodes'], body['reactions']
    metadata = header.get('metadata', source_metadata)
    if metadata is None:
        raise ValueError('Original export metadata is required to verify a saved map')
    source_ids = metadata['source_node_ids']
    routing_ids = set(metadata.get('routing_node_ids', []))
    if (set(source_ids) & routing_ids or set(source_ids) | routing_ids != set(map_nodes) or
            len(source_ids) != len(source_nodes) or set(source_ids.values()) != set(source_nodes)):
        raise ValueError('Map nodes were added or lost')
    orthogonal = metadata.get('routing') == 'orthogonal'
    if routing_ids and not orthogonal:
        raise ValueError('Routing nodes require orthogonal routing metadata')
    for node_id in routing_ids:
        node = map_nodes[node_id]
        if node['node_type'] != 'multimarker' or not all(math.isfinite(node[k]) for k in ('x','y')):
            raise ValueError('Invalid routing node')
    transform = metadata['coordinate_transform']
    from tca_curves import escher_curve_controls
    curved_count=0
    counts = Counter()
    for node_id, source_id in source_ids.items():
        node, source = map_nodes[node_id], source_nodes[source_id]
        expected_type = 'midmarker' if source['kind'] == 'reaction' else 'metabolite'
        if node['node_type'] != expected_type:
            raise ValueError('Node type mismatch')
        if expected_type == 'metabolite' and (node['bigg_id'] != source['original_id'] or
                node['node_is_primary'] != (source['kind'] == 'metabolite')):
            raise ValueError('Metabolite identity or primary flag mismatch')
        x, y = source_xy[source_id]
        expected = [transform['scale'] * (x - transform['x_min']) + transform['padding'],
                    transform['scale'] * (transform['y_max'] - y) + transform['padding']]
        if not all(math.isclose(node[axis], value, abs_tol=1e-6) for axis, value in zip(('x', 'y'), expected)):
            raise ValueError('Coordinate transform mismatch')
    expected_ids = {r['id'] for r in model['reactions']}
    if len(reactions) != len(expected_ids) or {r['bigg_id'] for r in reactions.values()} != expected_ids:
        raise ValueError('Reaction IDs missing or duplicated')
    original = {r['id']: r for r in model['reactions']}
    segment_ids = set()
    used_routing = set()
    for reaction_id, reaction in reactions.items():
        rid = reaction['bigg_id']
        source = original[rid]
        details = metadata['reactions'][rid]
        if details['map_reaction_id'] != reaction_id or source_ids.get(details['marker_node_id']) != 'rxn:' + rid:
            raise ValueError('Reaction ID metadata mismatch')
        for field in ('lower_bound', 'upper_bound'):
            value = source[field]
            if details['original_' + field] != (value if math.isfinite(value) else str(value)):
                raise ValueError('Flux bounds metadata mismatch')
        direction = reaction_direction(source['lower_bound'], source['upper_bound'])
        multiplier = -1 if direction == 'reverse' else 1
        if reaction['reversibility'] != (direction == 'reversible') or details['coefficient_multiplier'] != multiplier:
            raise ValueError('Direction was not preserved')
        coefficients = {m['bigg_id']: m['coefficient'] for m in reaction['metabolites']}
        if len(coefficients) != len(reaction['metabolites']) or coefficients != {m: multiplier * c for m, c in source['metabolites'].items()}:
            raise ValueError('Stoichiometry was not preserved')
        if details['original_stoichiometry'] != source['metabolites']:
            raise ValueError('Original stoichiometry metadata mismatch')
        participants = Counter()
        graph = {}
        for sid, segment in reaction['segments'].items():
            if sid in segment_ids:
                raise ValueError('Duplicate segment ID')
            segment_ids.add(sid)
            a, b = segment['from_node_id'], segment['to_node_id']
            if a not in map_nodes or b not in map_nodes or a == b:
                raise ValueError('Unknown or invalid segment endpoint')
            if orthogonal:
                pa, pb = map_nodes[a], map_nodes[b]
                # Escher may serialize straight lines as collinear Beziers.
                controls = [segment[k] for k in ('b1','b2') if segment[k] is not None]
                straight = any(math.isclose(pa[axis],pb[axis],abs_tol=1e-6,rel_tol=0) and
                               all(math.isfinite(c['x']) and math.isfinite(c['y']) and
                                   math.isclose(c[axis],pa[axis],abs_tol=1e-6,rel_tol=0) for c in controls)
                               for axis in ('x','y'))
                if not straight:
                    pair=metadata.get('segment_source_edges',{}).get(sid,[])
                    expected=escher_curve_controls(metadata,pair,[pa['x'],pa['y']],[pb['x'],pb['y']])
                    if (expected is None or a not in source_ids or b not in source_ids or
                            frozenset(pair)!=frozenset((source_ids[a],source_ids[b])) or
                            any(segment[k] is None or any(not math.isfinite(segment[k][axis]) or
                                not math.isclose(segment[k][axis],value,abs_tol=1e-5,rel_tol=0)
                                for axis,value in zip(('x','y'),point))
                                for k,point in zip(('b1','b2'),expected))):
                        raise ValueError('Segment is not horizontal or vertical, or a valid TCA circular Bezier')
                    curved_count+=1
            graph.setdefault(a,[]).append((b,sid))
            graph.setdefault(b,[]).append((a,sid))
            for endpoint in (a,b):
                if endpoint in source_ids:
                    counts[source_ids[endpoint]] += 1
        marker = details['marker_node_id']
        if marker not in graph:
            raise ValueError('Missing reaction marker connection')
        visited_segments = set()
        for node_id in graph:
            if node_id in routing_ids:
                if len(graph[node_id]) != 2:
                    raise ValueError('Routing node must join exactly two segments')
                continue
            node = source_nodes[source_ids[node_id]]
            if node['kind'] == 'reaction':
                if node_id != marker:
                    raise ValueError('Path attached to wrong reaction')
                continue
            if len(graph[node_id]) != 1:
                raise ValueError('Metabolite must terminate one reaction path')
            if node['kind'] == 'currency' and node['parent_id'] != 'rxn:' + rid:
                raise ValueError('Currency attached to wrong reaction')
            participants[node['original_id']] += 1
            previous, current = None, node_id
            seen = set()
            while current != marker:
                if current in seen:
                    raise ValueError('Cycle in reaction routing')
                seen.add(current)
                options = [(n,sid) for n,sid in graph[current] if n!=previous]
                if len(options) != 1:
                    raise ValueError('Branched or disconnected metabolite path')
                following, sid = options[0]
                visited_segments.add(sid)
                if orthogonal and set(metadata.get('segment_source_edges',{}).get(sid,[])) != {source_ids[marker],source_ids[node_id]}:
                    raise ValueError('Segment source incidence metadata mismatch')
                if current in routing_ids:
                    if current in used_routing:
                        raise ValueError('Routing node shared between reaction paths')
                    used_routing.add(current)
                if following != marker and following not in routing_ids:
                    raise ValueError('Path does not end at its reaction marker')
                previous,current = current,following
        if visited_segments != set(reaction['segments']):
            raise ValueError('Disconnected or unused reaction segments')
        if participants != Counter({m: 1 for m in source['metabolites']}):
            raise ValueError('Reaction participants lost or duplicated')
    if used_routing != routing_ids:
        raise ValueError('Routing nodes were lost or unused')
    for key, node in source_nodes.items():
        if node['kind'] == 'currency' and counts[key] != 1:
            raise ValueError('Currency copy must have one segment')
    crossing_check=None
    if metadata.get('crossing_routing',{}).get('verified'):
        from crossing_routing import validate_saved_geometry
        crossing_check=validate_saved_geometry(data,metadata)
    return dict(crossing_constraints=crossing_check,schema_valid=None, metadata_embedded='metadata' in header,
                reaction_ids_preserved=True, stoichiometry_preserved=True,
                coordinate_transform_verified=True, segment_endpoints_verified=True,
                nodes=len(map_nodes), source_nodes=len(source_nodes), routing_nodes=len(routing_ids),
                orthogonal_segments_verified=orthogonal and curved_count==0,
                non_tca_orthogonal_segments_verified=orthogonal,tca_curved_segments=curved_count, reactions=len(reactions), segments=len(segment_ids),
                primary_metabolites=sum(n['node_type'] == 'metabolite' and n['node_is_primary'] for n in map_nodes.values()),
                secondary_metabolites=sum(n['node_type'] == 'metabolite' and not n['node_is_primary'] for n in map_nodes.values()),
                currency_copy_counts=dict(Counter(n['original_id'] for n in source_nodes.values() if n['kind'] == 'currency')),
                reverse_only_reactions=[rid for rid, d in metadata['reactions'].items() if d['direction'] == 'reverse'],
                blocked_reactions=[rid for rid, d in metadata['reactions'].items() if d['direction'] == 'blocked'],
                missing_gene_rule_count=len(metadata['missing_gene_rule_reactions']))


def configure_label_display(data, major_reactions=()):
    """Keep structural label priorities independent of simulation colors."""
    meta = data[0].setdefault('metadata', {})
    options = meta.setdefault('label_display', {})
    for key, value in dict(show_metabolites=True, show_reactions=True, allow_overlap=False).items():
        options.setdefault(key, value)
    exchanges = set(meta.get('exchange_reaction_ids', ()))
    major = set(major_reactions)
    reactions = {}
    metabolites = {nid: 1 if not node.get('node_is_primary') else 2
                   for nid, node in data[1]['nodes'].items() if node['node_type'] == 'metabolite'}
    for rid, reaction in data[1]['reactions'].items():
        identity = reaction['bigg_id']
        rank = (5 if identity in exchanges or identity.upper().startswith('EX_') else
                4 if identity in major else
                2 + max((meta.get('segment_priorities', {}).get(sid, 0)
                         for sid in reaction['segments']), default=0))
        reactions[rid] = rank
        # Only actual terminal nodes inherit priority, so cofactor copies remain local.
        for segment in reaction['segments'].values():
            for nid in (segment['from_node_id'], segment['to_node_id']):
                if nid in metabolites and (data[1]['nodes'][nid].get('node_is_primary') or rank == 5):
                    metabolites[nid] = max(metabolites[nid], rank)
    options.update(reaction_priorities=reactions, metabolite_priorities=metabolites)
    return options


def configure_reaction_coverage(data, records):
    """Rank saved reactions red, orange, then other, using the sampled medians."""
    lookup = {record['reaction_id']: record for record in records}
    reactions = data[1]['reactions']
    if not reactions or not all(r['bigg_id'] in lookup and
                               lookup[r['bigg_id']].get('median') is not None for r in reactions.values()):
        return
    from display_ranking import reaction_display_order
    segment_priorities = data[0]['metadata'].get('segment_priorities', {})
    priorities = {reaction['bigg_id']: max((segment_priorities.get(sid, 0)
                   for sid in reaction.get('segments', {})), default=0)
                  for reaction in reactions.values()}
    ids = {reaction['bigg_id']: rid for rid, reaction in reactions.items()}
    ordered = reaction_display_order([lookup[identity] for identity in ids], priorities)
    data[0].setdefault('metadata', {})['reaction_coverage'] = dict(
        reaction_order=[ids[identity] for identity in ordered], initial_percent=100,
        ranking='Red first, then orange, then other reactions; descending absolute sampled median flux within each group; reaction ID breaks ties',
        denominator='Reactions present in this saved map')


def write_preview(path, data, reaction_data=None, fba_report=None):
    if 'reaction_priorities' not in data[0].get('metadata', {}).get('label_display', {}):
        configure_label_display(data)
    # Embedded library and map make this preview usable without a server or network.
    library = base64.b64encode(VIEWER.read_bytes()).decode('ascii')
    map_json = json.dumps(data, ensure_ascii=False, allow_nan=False).replace('<', '\\u003c')
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>Escher MILP layout</title>
<style>html,body{height:100%;margin:0;font:14px system-ui}header{padding:10px 18px;background:#f3f5f7}#map{height:calc(100vh - 100px);width:100%}#status{margin-left:16px}</style>
<header><b>e_coli_core: Escher layout</b><p>Reaction positions and CO2 copies are kept. This view uses Escher styles. Use the mouse wheel to zoom.</p>
<p id="escher-display-mode"></p>
<p id="escher-coverage-controls" hidden><label>Show reactions
<input id="escher-coverage" type="range" min="0" max="100" step="1" value="100">
<output id="escher-coverage-value">100%</output></label> <span id="escher-coverage-count"></span>
<br><small>Red first, then orange, then other reactions; within each group, descending absolute sampled median flux. Percentage of saved reactions.</small></p>
<p id="escher-label-controls"><label><input id="show-metabolite-labels" type="checkbox">Metabolite labels</label>
<label><input id="show-reaction-labels" type="checkbox">Reaction labels</label>
<label><input id="allow-label-overlap" type="checkbox">Allow label overlap</label></p><button id="save">Save map JSON</button><span id="status">Loading...</span></header><div id="map"></div>
<script>window.escherErrors=[];window.addEventListener('error',e=>window.escherErrors.push(e.message));</script>
<script src="data:application/javascript;base64,LIBRARY_PLACEHOLDER"></script><script>
const mapData=MAP_PLACEHOLDER;
window.builder=escher.Builder(mapData,null,null,escher.libs.d3_select('#map'),{
  enable_editing:true,enable_keys:true,enable_search:true,menu:'all',fill_screen:false,
  identifiers_on_map:'bigg_id',show_gene_reaction_rules:false,
  primary_metabolite_radius:10,secondary_metabolite_radius:3,marker_radius:3,
  reaction_no_data_size:3,hide_secondary_metabolites:false
});
function syncEscherDisplayMode(){
const loaded=!!window.builder.map?.has_data_on_reactions;
const colors=document.getElementById('major-pathway-colors');
if(colors?.sheet)colors.sheet.disabled=loaded;
for(const legend of document.querySelectorAll('[data-pathway-legend]'))legend.hidden=loaded;
document.getElementById('escher-display-mode').textContent=loaded?
  'Reaction data loaded: colors, width, values and arrows reflect data. Fixed pathway colors are off.':
  'No flux data loaded. Arrows show directions allowed by model bounds.';
window.syncEscherLayers?.();
}
const wait=setInterval(()=>{if(!window.builder.map)return;clearInterval(wait);window.escherReady=true;
window.builder.callback_manager.set('update_data.pathway_colors',syncEscherDisplayMode);syncEscherDisplayMode();
const m=window.builder.map;document.getElementById('status').textContent=Object.keys(m.reactions).length+' reactions, '+Object.keys(m.nodes).length+' nodes';},100);
document.getElementById('save').onclick=()=>{const data=window.builder.map.map_for_export();data[0].metadata=mapData[0].metadata;const url=URL.createObjectURL(new Blob([JSON.stringify(data,null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=url;a.download='escher_map.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
</script></html>'''
    from tca_curves import backbone_from_metadata
    routing=backbone_from_metadata(data[0].get('metadata',{})).get('crossing_routing',{}).get('report',{})
    port_keys=set(routing.get('node_port_nodes',[]))
    if port_keys and not routing.get('uniform_node_sizes'):
        radius=routing['node_port_glyph_radius_grid']*data[0]['metadata']['coordinate_transform']['scale']
        css=''.join('#n'+str(nid)+' circle{r:'+str(radius)+'px!important;}' for nid,key in data[0]['metadata']['source_node_ids'].items() if key in port_keys)
        page=page.replace('</style>',css+'</style>',1)
    if reaction_data is not None:
        if set(reaction_data) != {r['bigg_id'] for r in data[1]['reactions'].values()}:
            raise ValueError('Flux IDs must match all map reactions')
        if not all(math.isfinite(v) for v in reaction_data.values()):
            raise ValueError('Fluxes must be finite')
        maximum = max(abs(v) for v in reaction_data.values()) or 1.0
        active = [abs(v) for v in reaction_data.values() if v != 0]
        first = min(min(active), maximum * 0.02) if active else maximum * 0.02
        options = dict(reaction_data=reaction_data, reaction_styles=['color', 'size', 'text', 'abs'],
                       reaction_scale=[dict(type='value', value=0, color='#b8b8b8', size=1),
                                       dict(type='value', value=first, color='#2166ac' if fba_report else '#e6ab02', size=3),
                                       dict(type='value', value=maximum * 0.1, color='#67a9cf' if fba_report else '#ff7f00', size=5),
                                       dict(type='value', value=maximum * 0.5, color='#fdae61' if fba_report else '#e31a1c', size=9),
                                       dict(type='value', value=maximum, color='#d73027' if fba_report else '#9e0142', size=14)])
        if fba_report:
            options['reaction_scale'] = fba_flux_scale(maximum, first)
        page = page.replace('enable_editing:true,', json.dumps(options, allow_nan=False)[1:-1] + ',enable_editing:true,')
        page = page.replace('Escher layout', 'FBA flux in Escher')
        summary = ('Line width and color show absolute flux. Labels show signed values. '
                   'Zero is gray. Higher flux changes from yellow to orange, red, and purple. Thicker lines mean higher flux. Units: mmol/gDW/h. '
                   f'Scale: 0 to {maximum:.4g}. Arrows show the calculated flow direction. '
                   'Arrows on zero-flux reactions show allowed directions, not flow.')
        if fba_report:
            summary = summary.replace('Higher flux changes from yellow to orange, red, and purple.', 'Higher flux changes from blue through violet and magenta to red. Color and width use a square-root scale; legend ticks show actual absolute flux. Fixed pathway colors are not applied.')
            summary += (f' Biomass reaction: {fba_report["objective_value"]:.6f} h⁻¹. '
                        'This is one optimal solution that maximizes growth within the SBML bounds.')
        if fba_report:
            # The header legend uses readable physical ticks; the built-in
            # legend would label every interpolation stop.
            page = page.replace('</head>', '<style>#map .legend-container,#map .vertical-legend-container{display:none!important}</style></head>', 1)
            legend = '<p id="fba-color-legend" aria-label="Absolute FBA flux color scale">Absolute flux (mmol/gDW/h): '
            ticks = sorted({0, *[v for v in (1, 5, 10, 20) if first <= v < maximum], maximum})
            legend_stops = [fba_flux_style(v, maximum, first) for v in ticks]
            legend += ' &nbsp; '.join(
                f'<span><span aria-hidden="true" style="display:inline-block;width:1em;height:1em;background:{stop["color"]};border:1px solid #777"></span> {stop["value"]:.4g}</span>'
                for stop in legend_stops)
            legend += '</p>'
            page = page.replace('<p id="escher-display-mode">', legend + '<p id="escher-display-mode">', 1)
        page = page.replace('Reaction positions and CO2 copies are kept. This view uses Escher styles. Use the mouse wheel to zoom.', summary)
        page = page.replace('const mapData=MAP_PLACEHOLDER;',
                            'window.fbaFluxes=' + json.dumps(reaction_data, allow_nan=False) + ';const mapData=MAP_PLACEHOLDER;')
        page = page.replace('#map{height:calc(100vh - 100px)', '#map{height:calc(100vh - 140px)')
    if reaction_data is None:
        view_name = data[0]['map_name'].replace('flux-median-ranked grid map','grid map')
    else:
        view_name = data[0]['map_name']
    page = page.replace('e_coli_core: Escher layout', html.escape(view_name))
    page = page.replace('e_coli_core: FBA flux in Escher', html.escape(data[0]['map_name']) + ': FBA flux in Escher')
    page = page.replace('</html>','<script>'+BASE.joinpath('escher_layers.js').read_text()+'</script><script>'+BASE.joinpath('escher_labels.js').read_text()+'</script><script>'+BASE.joinpath('escher_coverage.js').read_text()+'</script></html>')
    path.write_text(page.replace('LIBRARY_PLACEHOLDER', library).replace('MAP_PLACEHOLDER', map_json))



def write_pathway_preview(path, data, major_reactions=(), connector_pairs=None,
                          model=None, elements=(), initial_highlights=None, records=()):
    """Start with pathway structure; disable fixed colors when reaction data loads."""
    if data[0]['metadata'].get('routing')=='orthogonal' and (not data[0]['metadata'].get('overlap_routing',{}).get('enabled') or
            data[0]['metadata'].get('overlap_routing',{}).get('node_obstacle_policy')!='red_orange_structural_nodes' or
            data[0]['metadata'].get('overlap_routing',{}).get('trunk_grouping')!='positive_length_overlap' or
            data[0]['metadata'].get('overlap_routing',{}).get('red_layout_signature')!=(plan_from_metadata(data[0]['metadata']) or {}).get('signature') or
            data[0]['metadata'].get('overlap_routing',{}).get('routing_preference')!=POLICY):
        data = separate_routes(copy.deepcopy(data),major_reactions,connector_pairs)
    configure_label_display(data, major_reactions)
    configure_reaction_coverage(data, records)
    write_preview(path, data)
    metadata = data[0]['metadata']
    source_ids = metadata['source_node_ids']
    major_reactions = set(major_reactions)
    css = ['#reactions .segment {stroke-width:3 !important;}']
    for rid, reaction in data[1]['reactions'].items():
        if connector_pairs is not None:
            css.extend([f'#r{rid} .segment {{stroke:#cbd3d8 !important;}}',
                        f'#r{rid} .arrowheads path {{stroke:#cbd3d8 !important;fill:#cbd3d8 !important;}}'])
            for sid, segment in reaction['segments'].items():
                pair = metadata.get('segment_source_edges',{}).get(sid)
                if pair is None:
                    pair = [source_ids[segment['from_node_id']],source_ids[segment['to_node_id']]]
                if reaction['bigg_id'] not in major_reactions and frozenset(pair) in connector_pairs:
                    css.extend([f'#r{rid} #s{sid} .segment {{stroke:#f28e2b !important;stroke-width:5 !important;}}',
                                f'#r{rid} #s{sid} .arrowheads path {{stroke:#f28e2b !important;fill:#f28e2b !important;}}'])
        if reaction['bigg_id'] in major_reactions:
            css.extend([f'#r{rid} .segment {{stroke:#d62728 !important;stroke-width:6 !important;}}',
                        f'#r{rid} .arrowheads path {{stroke:#d62728 !important;fill:#d62728 !important;}}',
                        f'#r{rid} text {{fill:#d62728 !important;}}'])
    blocked = '\n'.join(f'#r{rid} .arrowheads {{display:none !important;}}'
                        for rid, r in data[1]['reactions'].items()
                        if metadata['reactions'][r['bigg_id']]['direction']=='blocked')
    style = '<style id="major-pathway-colors">'+'\n'.join(css)+'</style>'
    style += '<style id="blocked-reaction-arrows">'+blocked+'</style>'
    legend = '<p data-pathway-legend>Red: main pathways.'
    if connector_pairs is not None:
        legend += ' Orange: connections to exchange reactions. Gray: other branches.'
    connector_details = metadata.get('sbml2map', metadata.get('sbml2escher', {})).get('major_backbones', {}).get('exchange_connectors', {})
    flux_weight = connector_details.get('flux_costs', {}).get('flux_weight', 0)
    if flux_weight:
        legend += f' Orange route selection favors larger absolute sampled median fluxes (penalty strength {flux_weight:g}).'
    legend += ' These colors describe structure, not flux.</p>'
    page = path.read_text().replace('<header>',style+'<header>'+legend,1)
    if (path.parent/'escher_fba.html').exists():
        page = page.replace('<header>','<header><p id="fba-result-link"><a href="escher_fba.html">View FBA result</a></p>',1)
    if model is not None:
        from highlighting import inject_controls
        page = inject_controls(page,model,elements,initial_highlights,records,escher=True)
    path.write_text(page)


def export_layout(directory, output, scale=4000.0, padding=180.0, preview=None):
    model, nodes, edges, xy = read_layout(directory)
    report_path = directory / 'report.json'
    layout_method = json.loads(report_path.read_text()).get('layout_method', 'MILP reaction layout') if report_path.exists() else 'MILP reaction layout'
    data = create_map(model, nodes, edges, xy, scale, padding, layout_method=layout_method)
    backbone_file = directory / 'major_pathway.json'
    if backbone_file.exists():
        backbone = json.loads(backbone_file.read_text())
        if 'pathways' in backbone:
            data[0]['metadata']['selected_pathways'] = backbone['pathways']
            data[0]['metadata']['pathway_emphasis'] = 'Representative pathway selection, not flux magnitude; original reversibility is retained.'
            transform = data[0]['metadata']['coordinate_transform']
            labels = data[1]['text_labels']
            for group in backbone['pathways']:
                x, y = group['label_position']
                labels[str(len(labels))] = dict(x=scale * (x - transform['x_min']) + padding,
                                               y=scale * (transform['y_max'] - y) + padding,
                                               text=group['name'])
            labels[str(len(labels))] = dict(x=data[1]['canvas']['x'] + padding,
                                           y=data[1]['canvas']['y'] + data[1]['canvas']['height'] - padding,
                                           text='Selected pathways are layout backbones, not flux magnitudes. ETC: electron transfer; oxygen is reduced to water. ATP synthase is separate.')
    if backbone_file.exists() and 'compartment_regions' in backbone:
        data[0]['metadata']['compartment_layout'] = dict(regions=backbone['compartment_regions'],
                                                        species_symbols=backbone['species_compartments'],
                                                        audit=backbone['compartment_audit'])
        transform = data[0]['metadata']['coordinate_transform']
        labels = data[1]['text_labels']
        for symbol, region in backbone['compartment_regions'].items():
            if not any(backbone['species_compartments'][m['id']] == symbol for m in model['metabolites']):
                continue
            labels[str(len(labels))] = dict(x=scale*(region['left']-transform['x_min'])+padding,
                                           y=scale*(transform['y_max']-region['top'])+padding,
                                           text=f"[{symbol}] {region['name']}")
    verification = validate_map(data, model, nodes, xy)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    output.with_name(output.stem + '_metadata.json').write_text(json.dumps(
        dict(map_id=data[0]['map_id'], metadata=data[0]['metadata']), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    report = dict(source_layout=str(directory), map_file=output.name, schema='Escher 1-0-0',
                  validation=verification, browser_validation='not run',
                  limitations=['Escher uses its own styles; original square markers and pathway colors are not reproduced.',
                               'Bounds and blocked status are metadata; gene rules absent from the loaded model are empty.',
                               'Native Escher 1.8.1 saves omit additional metadata; retain the separate metadata JSON.'])
    output.with_name(output.stem + '_report.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    if preview:
        preview.parent.mkdir(parents=True, exist_ok=True)
        write_preview(preview, data)
    print(json.dumps(verification, indent=2, ensure_ascii=False))
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--layout-dir', type=Path, default=BASE / 'layout_reaction_sbml_output')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--scale', type=float, default=4000)
    parser.add_argument('--padding', type=float, default=180)
    parser.add_argument('--preview', type=Path, help='Write a self-contained Escher HTML preview')
    args = parser.parse_args()
    output = args.output or args.layout_dir / 'escher_map.json'
    try:
        export_layout(args.layout_dir, output, args.scale, args.padding, args.preview)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
