"""Complete FBA incidence layout, retaining saved primary pathway coordinates."""
import math

from escher_grid import classify_metabolites
from reaction_network import build_network, validate_network


def complete_layout(model, saved_nodes, saved_xy):
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
                    if 0 <= p[0] <= cols-1 and 0 <= p[1] <= rows-1:
                        distance = min(math.dist(p, q) for q in xy.values())
                        candidates.append((distance, -radius, -index, p))
            if not candidates:
                raise ValueError('No space for reaction-local participant: ' + key)
            xy[key] = max(candidates)[3]
    validation = validate_network(model, nodes, edges, xy)
    validation.update(existing_nodes_preserved=len(set(saved_xy) & set(xy)),
                      added_nodes=len(set(xy)-set(saved_xy)),
                      existing_coordinates_changed=any(xy.get(k) != list(p) for k, p in saved_xy.items()),
                      all_participants_drawn=True, grid_cols=cols, grid_rows=rows)
    return nodes, edges, xy, validation
