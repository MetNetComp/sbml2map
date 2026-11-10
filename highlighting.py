"""Shared, offline highlight controls for the grid and Escher HTML views."""
import json
from pathlib import Path

BLUE = '#1677ff'


def validate_selection(value, model):
    if not isinstance(value, dict) or set(value) - {'metabolites', 'reactions'}:
        raise ValueError('Highlights must contain metabolites and/or reactions')
    result = {}
    for kind in ('metabolites', 'reactions'):
        ids = value.get(kind, [])
        if not isinstance(ids, list) or any(not isinstance(i, str) for i in ids):
            raise ValueError(f'Highlight {kind} must be a list of IDs')
        unknown = set(ids) - {item['id'] for item in model[kind]}
        if unknown:
            raise ValueError(f'Unknown highlight {kind}: {", ".join(sorted(unknown))}')
        result[kind] = sorted(set(ids))
    return result


def read_selection(args, model):
    value = json.loads(args.highlight_file.read_text()) if getattr(args, 'highlight_file', None) else {}
    value = validate_selection(value, model)
    for kind in ('metabolites', 'reactions'):
        value[kind] += getattr(args, 'highlight_' + kind, None) or []
    return validate_selection(value, model)


def grid_elements(nodes, edges):
    result = []
    keys = {key: j for j, key in enumerate(nodes)}
    for j, node in enumerate(nodes.values()):
        result.append(dict(kind='reactions' if node['kind'] == 'reaction' else 'metabolites',
                           key=node['original_id'], ids=[f'node-{j}', f'label-{j}'], role='node'))
    for j, edge in enumerate(edges):
        result.append(dict(kind='reactions', key=edge['reaction_id'], role='edge',
                           ids=[f'edge-{j}', f'arrow-{j}', f'reverse-arrow-{j}'],
                           endpoints=[f'{prefix}-{keys[k]}' for k in (edge['source'], edge['target'])
                                      for prefix in ('node', 'label')]))
    return result


def catalog(model, elements, records=()):
    drawn = {(e['kind'], e['key']) for e in elements}
    reasons = {r['reaction_id']: r.get('exclusion_reason', '') for r in records}
    return {kind: [dict(id=item['id'], name=item.get('name', ''),
                        compartment=item.get('compartment', ''),
                        drawn=(kind, item['id']) in drawn,
                        reason=reasons.get(item['id'], 'not_in_saved_map') if kind == 'reactions' else 'not_in_saved_map')
                   for item in model[kind]] for kind in ('metabolites', 'reactions')}


def inject_controls(page, model, elements=(), initial=None, records=(), escher=False):
    """Append controls without replacing existing layout, styles or viewer code."""
    if 'id="highlight-panel"' in page:
        return page
    if escher:
        # Escher's live map supplies element IDs, including nodes created by editing.
        elements = list(elements)
    payload = dict(catalog=catalog(model, elements, records), elements=elements,
                   initial=initial or dict(metabolites=[], reactions=[]), escher=escher)
    encoded = json.dumps(payload, ensure_ascii=False).replace('<', '\\u003c')
    panel = '''<section id="highlight-panel" aria-label="Highlight compounds and reactions">
<p><b>Blue: user selection</b></p>
<label>Type <select id="hl-kind"><option value="metabolites">Compounds</option><option value="reactions">Reactions</option></select></label>
<label>Search ID or name <input id="hl-search" type="search"></label>
<div id="hl-results" aria-label="Search results"></div>
<label>IDs (comma or newline separated) <textarea id="hl-bulk" rows="2"></textarea></label>
<button id="hl-add" type="button">Add</button>
<label><input id="hl-force" type="checkbox">Always show blue targets</label>
<div id="hl-selected" aria-label="Selected targets"></div>
<button id="hl-clear" type="button">Clear</button>
<button id="hl-save" type="button">Save selection JSON</button>
<label>Load selection JSON <input id="hl-load" type="file" accept=".json,application/json"></label>
<button id="hl-svg" type="button">Save current SVG</button>
<p id="hl-message" role="status" aria-live="polite"></p>
<p>Selection changes this HTML view. PNG files are updated by running the CLI.</p>
</section>'''
    style = '''<style id="highlight-style">
#highlight-panel{padding:10px;border:1px solid #ccd;margin:10px 0;font:14px system-ui}
#highlight-panel label{display:inline-block;margin:4px}#highlight-panel button{margin:3px}
#hl-results{max-height:130px;overflow:auto}#hl-results button{display:block;text-align:left}
#hl-selected{max-height:130px;overflow:auto}#hl-bulk{vertical-align:middle}
body #viewer .user-blue path,body #viewer .user-blue use,body #viewer .user-blue text{stroke:#1677ff!important;fill:#1677ff!important}
body #viewer .user-blue[data-hl-role="edge"] path{fill:none!important;stroke:#1677ff!important;stroke-opacity:1!important}
body #viewer .user-blue[data-hl-role="edge"]{opacity:1!important}
body #map#map .user-blue .segment,body #map#map .user-blue .arrowheads path{stroke:#1677ff!important;stroke-opacity:1!important}
body #map#map .user-blue .arrowheads path,body #map#map .user-blue circle,body #map#map .user-blue text{fill:#1677ff!important}
</style>'''
    js = Path(__file__).with_name('highlight_controls.js').read_text()
    injection = style + panel + '<script>\n(()=>{const config=' + encoded + ';\n' + js + '\n})();</script>'
    return page.replace('</html>', injection + '</html>') if '</html>' in page else page + injection
