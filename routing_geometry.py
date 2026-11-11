"""One canonical grid-coordinate route set for SVG, Escher and SBML Layout."""
import copy
from parallel_routing import separate_parallel, avoid_nodes, structural_node_priorities
from red_route_layout import spacing_report
from routing_simplification import simplify_important, yield_gray


def build_routes(nodes, edges, xy, records, selection, backbones=None, route_offset=.3):
    backbones = backbones or {}
    stats = {r['reaction_id']: r for r in records}
    structural = selection.get('mode') == 'structure'
    connectors = backbones.get('exchange_connectors') or {}
    orange_edges = set(connectors.get('edge_indices', ()))
    orange_pairs = {frozenset((edges[i]['source'], edges[i]['target'])) for i in orange_edges}
    major = {rid for rid, record in stats.items() if record.get('is_major_pathway', False)}
    major_mids = {mid for route in backbones.get('routes', ()) for mid in route['metabolite_ids']}
    priorities = structural_node_priorities(nodes, major, major_mids, orange_pairs)
    obstacle_keys = sorted(key for key, rank in priorities.items() if rank > 0)
    obstacles = [xy[key] for key in obstacle_keys]
    backbone_pairs = {frozenset(('rxn:'+step['reaction_id'], 'met:'+mid))
                      for route in backbones.get('routes', ()) for step in route.get('steps', ())
                      for mid in (step['source'], step['target'])}
    paths = []
    for edge in edges:
        rid = edge['reaction_id']
        source, target = edge['source'], edge['target']
        reverse = edge['direction'] == 'reverse' if structural else stats[rid]['median'] < 0
        if reverse:
            source, target = target, source
        a, b = tuple(xy[source]), tuple(xy[target])
        elbow = (b[0], a[1])
        points = [a]+([elbow] if elbow != a and elbow != b else [])+[b]
        paths.append(dict(rid=rid, mid=edge['metabolite_id'], source=source, target=target,
                          points=points, priority=2 if rid in major else 1 if frozenset((source, target)) in orange_pairs else 0,
                          backbone=frozenset((source, target)) in backbone_pairs))
    original_length = sum(sum(abs(a[k]-b[k]) for k in (0, 1))
                          for path in paths for a, b in zip(path['points'], path['points'][1:]))
    red_plan = backbones.get('red_parallel_layout')
    if red_plan and red_plan.get('enabled') and red_plan.get('verified'):
        planned = {(r['rid'], r['mid']): r for r in red_plan['routes']}
        for path in paths:
            reserved = planned.get((path['rid'], path['mid']))
            if reserved:
                points = [tuple(q) for q in reserved['points']]
                if points[0] != tuple(xy[path['source']]):
                    points.reverse()
                path.update(points=points, routing_locked=True)
    report = separate_parallel(paths, obstacles=obstacles)
    report['node_detours'] = avoid_nodes(paths, obstacles=obstacles, route_offset=route_offset)
    spacing = red_plan['requested_spacing_grid'] if red_plan and red_plan.get('enabled') else .3
    report['important_simplification'] = simplify_important(paths, obstacles=obstacles, route_offset=route_offset, red_spacing=spacing)
    report['gray_yield'] = yield_gray(paths, obstacles=obstacles, route_offset=route_offset)
    report.update(node_obstacle_policy='red_orange_structural_nodes', obstacle_node_ids=obstacle_keys,
                  structural_node_priorities=priorities, original_route_length=original_length)
    if red_plan:
        report['red_parallel_layout'] = {k: v for k, v in red_plan.items() if k not in ('routes', 'allocations')}
        report['red_parallel_spacing'] = spacing_report(paths, red_plan['requested_spacing_grid'])
        if red_plan.get('enabled') and red_plan.get('verified') and not report['red_parallel_spacing']['verified']:
            raise RuntimeError('Preallocated red spacing changed during canonical routing')
    return paths, report


def orient_for_edges(paths, edges, xy):
    lookup = {(p['rid'], frozenset((p['source'], p['target']))): p for p in paths}
    ordered = []
    for edge in edges:
        path = copy.deepcopy(lookup[(edge['reaction_id'], frozenset((edge['source'], edge['target'])))])
        ordered.append(path)
    return ordered


def transformed_paths(paths, transform):
    scale = transform['scale']
    return [dict(path, points=[(scale*(x-transform['x_min'])+transform['padding'],
                              scale*(transform['y_max']-y)+transform['padding'])
                             for x, y in path['points']]) for path in paths]
