"""Complete FBA incidence layout, retaining saved primary pathway coordinates."""
import math

from escher_grid import classify_metabolites
from reaction_network import build_network, validate_network


def complete_layout(model, saved_nodes, saved_xy, exclusion_region=None, exchange_alignment='off', grid_shape=None):
    classification = classify_metabolites(model)
    shared = {n['original_id'] for n in saved_nodes.values() if n['kind'] == 'metabolite'}
    currency = {mid for mid, record in classification.items()
                if record['kind'] in ('inorganic', 'coenzyme')} - shared
    currency.update(n['original_id'] for n in saved_nodes.values() if n['kind'] == 'currency')
    nodes, edges = build_network(model, currency)
    xy = {key: list(value) for key, value in saved_xy.items() if key in nodes}
    occupied = {tuple(value) for key, value in xy.items() if nodes[key]['kind'] != 'currency'}
    cols = max(2, math.ceil(max(p[0] for p in xy.values())) + 1)
    rows = max(2, math.ceil(max(p[1] for p in xy.values())) + 1)
    if grid_shape is not None:
        cols,rows=grid_shape
    if exclusion_region is not None or exchange_alignment == 'collinear':
        from reserved_grid import optimize_reserved_grid
        from exchange_alignment import discover_alignment
        grid_ids=[k for k,node in nodes.items() if node['kind']!='currency']
        boundary={k for k in grid_ids if nodes[k]['kind']=='reaction' and nodes[k]['original_id'].startswith('EX_')}
        groups,skipped=discover_alignment(model,nodes,boundary) if exchange_alignment=='collinear' else ([],[])
        fixed={k:[p[0],p[1]-(rows-1)] for k,p in xy.items() if k in grid_ids}
        region=dict(exclusion_region) if exclusion_region is not None else None
        if region:
            if region['shape']=='circle':region['center_offset']=[region['center'][0],region['center'][1]-(rows-1)]
            else:region['top_offset']=region['top']-(rows-1)
        pairs=[(e['source'],e['target']) for e in edges if e['source'] in grid_ids and e['target'] in grid_ids]
        placed,placement=optimize_reserved_grid(grid_ids,pairs,cols,rows,time_limit=0,
            fixed_offsets=fixed,boundary_nodes=boundary,max_milp_nodes=1,
            reserved_region=region,alignment_groups=groups)
        xy.update(placed)
        occupied={tuple(xy[k]) for k in grid_ids}
    adjacency = {key: [] for key in nodes}
    for edge in edges:
        a, b = edge['source'], edge['target']
        adjacency[a].append(b)
        adjacency[b].append(a)
    # Reactions with known participants first; newly visible organic species follow.
    pending = {key for key in nodes if key not in xy and nodes[key]['kind'] != 'currency'}
    while pending:
        key = min(pending, key=lambda k: (-sum(n in xy for n in adjacency[k]), k))
        neighbors = [xy[n] for n in adjacency[key] if n in xy]
        center = ([sum(p[i] for p in neighbors) / len(neighbors) for i in (0, 1)]
                  if neighbors else [(cols - 1) / 2, (rows - 1) / 2])
        exchange = nodes[key]['kind'] == 'reaction' and nodes[key]['original_id'].startswith('EX_')
        def cells():
            return [(x, y) for x in range(cols) for y in range(rows)
                    if (x, y) not in occupied and
                    (not exchange or x in (0, cols - 1) or y in (0, rows - 1))]
        candidates = cells()
        while not candidates:
            cols += 1
            candidates = cells()
        position = min(candidates, key=lambda p: (abs(p[0]-center[0])+abs(p[1]-center[1]), p))
        xy[key] = list(position)
        occupied.add(position)
        pending.remove(key)
    # Local copies preserve chemical IDs while avoiding long cofactor hubs.
    for parent in sorted(key for key in nodes if nodes[key]['kind'] == 'reaction'):
        copies = sorted(key for key, node in nodes.items()
                        if node['kind'] == 'currency' and node['parent_id'] == parent and key not in xy)
        for key in copies:
            candidates = []
            for radius in (.45, .65, .85):
                for index in range(128):
                    angle = 2 * math.pi * (index + .25) / 128
                    dx, dy = math.cos(angle), math.sin(angle)
                    norm = abs(dx) + abs(dy)
                    p = [xy[parent][0] + radius*dx/norm, xy[parent][1] + radius*dy/norm]
                    from tca_geometry import allowed_points
                    if 0 <= p[0] <= cols-1 and 0 <= p[1] <= rows-1 and allowed_points([p],exclusion_region)[0]:
                        distance = min(math.dist(p, q) for q in xy.values())
                        candidates.append((distance, -radius, -index, p))
            if not candidates:
                raise ValueError('No space for reaction-local participant: ' + key)
            xy[key] = max(candidates)[3]
    validation = validate_network(model, nodes, edges, xy)
    if exclusion_region is not None:
        from tca_geometry import validate_exclusion
        validation['tca_interior']=validate_exclusion(nodes,xy,exclusion_region)
    if exchange_alignment=='collinear':
        from exchange_alignment import validate_alignment
        validation['exchange_alignment']=validate_alignment(groups,xy,cols,rows)
        validation['exchange_alignment']['skipped']=skipped
    validation.update(existing_nodes_preserved=len(set(saved_xy) & set(xy)),
                      added_nodes=len(set(xy)-set(saved_xy)),
                      existing_coordinates_changed=any(xy.get(k) != list(p) for k, p in saved_xy.items()),
                      all_participants_drawn=True, grid_cols=cols, grid_rows=rows)
    return nodes, edges, xy, validation
