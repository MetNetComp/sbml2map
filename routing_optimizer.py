"""Bounded local routing search with explicit geometric and biological constraints."""
from collections import defaultdict
import copy
import math
import time
import numpy as np
from parallel_routing import EPS, _legs, _distance, _overlap
from routing_simplification import clean, candidates
from routing_metrics import MetricIndex, objective, route_key, validate_options


class ObstacleIndex(dict):
    def __init__(self, positions):
        super().__init__(positions)
        self.keys_array = list(positions)
        self.points_array = np.asarray(list(positions.values()), dtype=float).reshape((-1, 2))


def node_hits(path, obstacles):
    if not obstacles:
        return set()
    index = obstacles if isinstance(obstacles, ObstacleIndex) else ObstacleIndex(obstacles)
    points = clean(path['points'])
    positions = index.points_array
    eligible = np.ones(len(positions), dtype=bool)
    for endpoint in (points[0], points[-1]):
        eligible &= np.any(np.abs(positions-endpoint) > EPS, axis=1)
    hits = np.zeros(len(positions), dtype=bool)
    for leg in _legs(points):
        axis = leg['axis']
        along = np.maximum(np.maximum(leg['lo']-positions[:, axis], 0), positions[:, axis]-leg['hi'])
        side = positions[:, 1-axis]-leg['level']
        hits |= along*along+side*side < (.12-EPS)**2
    return {index.keys_array[i] for i in np.flatnonzero(hits & eligible)}


def new_overlaps(events, index):
    """Existing overlap may shrink, but no candidate may hide a crossing in a new overlap."""
    for key, (_, intervals) in events.items():
        previous = index.pairs.get(key, (set(), set()))[1]
        for axis, level, lo, hi in intervals:
            if not any(axis == a and abs(level-c) < EPS and lo >= l-EPS and hi <= h+EPS
                       for a, c, l, h in previous):
                return True
    return False


def red_valid(paths, others, spacing):
    for path in paths:
        if path.get('priority', 0) != 2:
            continue
        for other in others:
            if path['rid'] == other['rid'] or other.get('priority', 0) != 2:
                continue
            for a in _legs(path['points']):
                for b in _legs(other['points']):
                    if (a['axis'] == b['axis'] and _overlap(a, b) >= 1-EPS and
                            abs(a['level']-b['level']) < spacing-EPS):
                        return False
    return True


def shares(left, right):
    return any(a['axis'] == b['axis'] and abs(a['level']-b['level']) < EPS and _overlap(a, b) > EPS
               for a in _legs(left['points']) for b in _legs(right['points']))


def simple(points):
    points = clean(points)
    try:
        legs = list(_legs(points))
    except ValueError:
        return False
    for i, a in enumerate(legs):
        for b in legs[i+1:]:
            if a['axis'] == b['axis'] and abs(a['level']-b['level']) < EPS and _overlap(a, b) > EPS:
                return False
            if a['axis'] != b['axis']:
                h, v = (a, b) if a['axis'] == 0 else (b, a)
                if h['lo']+EPS < v['level'] < h['hi']-EPS and v['lo']+EPS < h['level'] < v['hi']-EPS:
                    return False
    return True


def preserve_trunks(original, proposed):
    for i, left in enumerate(original):
        for j, right in enumerate(original[i+1:], i+1):
            if shares(left, right) and not shares(proposed[i], proposed[j]):
                return False
    return True


def group_candidates(group, offset, blockers=()):
    # Change a positively shared trunk as a whole, never one branch in isolation.
    shared = {i for i, a in enumerate(group) if any(i != j and shares(a, b) for j, b in enumerate(group))}
    for i, path in enumerate(group):
        if path.get('routing_locked') or path.get('backbone') or i in shared:
            continue
        a, b = path['points'][0], path['points'][-1]
        points_to_avoid = [point for other in blockers for point in other['points']]
        xs = sorted({point[0]+sign*max(.15,offset) for point in points_to_avoid for sign in (-1,1)}, key=lambda x: (min(abs(x-a[0]),abs(x-b[0])),x))[:12]
        ys = sorted({point[1]+sign*max(.15,offset) for point in points_to_avoid for sign in (-1,1)}, key=lambda y: (min(abs(y-a[1]),abs(y-b[1])),y))[:12]
        extra = [[a,(x,a[1]),(x,b[1]),b] for x in xs]+[[a,(a[0],y),(b[0],y),b] for y in ys]
        for points in list(candidates(path, max(.15, offset), 1))+[clean(points) for points in extra]:
            if points == clean(path['points']):
                continue
            new = [dict(p) for p in group]
            new[i]['points'] = points
            yield new
    if any(p.get('routing_locked') or p.get('backbone') for p in group):
        return
    for vertical_first in (False, True):
        new = []
        for path in group:
            a, b = path['points'][0], path['points'][-1]
            elbow = (a[0], b[1]) if vertical_first else (b[0], a[1])
            new.append(dict(path, points=clean([a, elbow, b])))
        if preserve_trunks(group, new):
            yield new


def optimize(paths, xy, nodes=None, movable=(), cols=None, rows=None, valid_position=None,
             mode='crossings-first', crossing_weight=10., bend_weight=1., time_limit=20.,
             route_offset=.3, red_spacing=1., node_moves=True):
    """Never worsen the requested objective; retain incumbent on timeout.

    Constraints do not silently disappear when the incumbent is imperfect:
    no new node hits or unrelated overlaps are permitted. Existing violations
    are reported explicitly. Red channels and biochemical backbone shapes stay
    fixed. A bounded heuristic does not establish global optimality.
    """
    validate_options(mode, crossing_weight, bend_weight, time_limit)
    started = time.monotonic()
    deadline = started + time_limit
    paths = copy.deepcopy(paths)
    xy = {k: list(p) for k, p in xy.items()}
    index = MetricIndex(paths)
    before = index.metrics()
    initial_xy = copy.deepcopy(xy)
    obstacle_index = ObstacleIndex(xy)
    initial_node_hits = sum(len(node_hits(p, obstacle_index)) for p in paths)
    hit_cache = {}
    accepted_routes, accepted_nodes, evaluated, passes = 0, 0, 0, 0
    score = lambda m: objective(m, mode, crossing_weight, bend_weight)

    def assess(replacements, positions=None):
        nonlocal evaluated
        evaluated += 1
        obstacles = ObstacleIndex(positions) if positions is not None else obstacle_index
        for rid, group in replacements.items():
            old = index.groups[rid]
            if len(old) != len(group) or not preserve_trunks(old, group):
                return None
            for a, b in zip(old, group):
                if b['points'] != a['points'] and not simple(b['points']):
                    return None
                if positions is None and b['points'] == a['points']:
                    continue
                cache_key = (route_key(a), tuple(map(tuple, a['points'])))
                if cache_key not in hit_cache:
                    hit_cache[cache_key] = node_hits(a, obstacle_index)
                if node_hits(b, obstacles) - hit_cache[cache_key]:
                    return None
        changed = [p for group in replacements.values() for p in group]
        others = [p for rid, group in index.groups.items() if rid not in replacements for p in group]
        if not red_valid(changed, others+changed, red_spacing):
            return None
        metrics, events, geometries = index.evaluate(replacements)
        if new_overlaps(events, index) or score(metrics) >= score(index.metrics()):
            return None
        return metrics, events, geometries

    if route_offset > 0 and time_limit > 0:
        # Alternate route choices and legal one-cell moves/swaps; keep the best
        # incumbent after every accepted operation, including on timeout.
        for _ in range(3):
            if time.monotonic() >= deadline:
                break
            passes += 1
            improved = False
            order = sorted(index.groups, key=lambda rid: (-sum(len(v[0]) for pair, v in index.pairs.items() if rid in pair)/max(1,sum(len(p['points'])-1 for p in index.groups[rid]))**2, rid))
            route_deadline = min(deadline, time.monotonic()+max(0, deadline-time.monotonic())*(.6 if node_moves and movable else 1.))
            for rid in order:
                if time.monotonic() >= route_deadline:
                    break
                best = None
                crossing_rids = {b if a == rid else a for (a,b), events in index.pairs.items() if rid in (a,b) and events[0]}
                blockers = [p for r in sorted(crossing_rids) for p in index.groups[r]]
                for group in group_candidates(index.groups[rid], route_offset, blockers):
                    if time.monotonic() >= route_deadline:
                        break
                    result = assess({rid: group})
                    if result and (best is None or score(result[0]) < score(best[1][0])):
                        best = (group, result)
                if best:
                    group, (_, events, geometries) = best
                    index.commit({rid: group}, events, geometries)
                    accepted_routes += 1
                    improved = True
            if node_moves and nodes and time.monotonic() < deadline:
                occupied = {tuple(p): key for key, p in xy.items() if nodes[key]['kind'] != 'currency'}
                eligible = set(movable)
                # Nodes with many crossing incidences first.
                incidence = defaultdict(set)
                for rid, group in index.groups.items():
                    for path in group:
                        incidence[path.get('source')].add(rid)
                        incidence[path.get('target')].add(rid)
                for key in sorted(eligible, key=lambda k: (-len(incidence[k]), k)):
                    if time.monotonic() >= deadline:
                        break
                    x, y = xy[key]
                    for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                        if time.monotonic() >= deadline:
                            break
                        target = (x+dx, y+dy)
                        if cols is not None and not (0 <= target[0] < cols and 0 <= target[1] < rows):
                            continue
                        partner = occupied.get(target)
                        if partner is not None and partner not in eligible:
                            continue
                        shifts = {key: [target[0], target[1]]}
                        if partner:
                            shifts[partner] = [x, y]
                        if valid_position and any(not valid_position(k, pos) for k, pos in shifts.items()):
                            continue
                        # Reaction-local cofactor copies move rigidly with their parent.
                        for k, node in nodes.items():
                            parent = node.get('parent_id')
                            if node['kind'] == 'currency' and parent in shifts:
                                delta = [shifts[parent][i]-xy[parent][i] for i in (0, 1)]
                                shifts[k] = [xy[k][i]+delta[i] for i in (0, 1)]
                        positions = dict(xy, **shifts)
                        if cols is not None and any(not (0 <= p[0] <= cols-1 and 0 <= p[1] <= rows-1) for p in shifts.values()):
                            continue
                        if any(math.dist(p, positions[other]) < .12-EPS for k, p in shifts.items() for other in positions if other != k):
                            continue
                        affected = set().union(*(incidence[k] for k in shifts))
                        replacements = {}
                        for rid in affected:
                            group = []
                            for path in index.groups[rid]:
                                if path.get('source') in shifts or path.get('target') in shifts:
                                    if path.get('routing_locked') or path.get('backbone'):
                                        group = None
                                        break
                                    a, b = positions[path['source']], positions[path['target']]
                                    group.append(dict(path, points=clean([a, (b[0], a[1]), b])))
                                else:
                                    group.append(dict(path))
                            if group is None:
                                break
                            replacements[rid] = group
                        if len(replacements) != len(affected):
                            continue
                        # A moved node must not newly intersect any unchanged route.
                        blocked = False
                        for rid, group in index.groups.items():
                            if rid not in affected:
                                for path in group:
                                    if node_hits(path, shifts)-node_hits(path, {k: xy[k] for k in shifts}):
                                        blocked = True
                                        break
                            if blocked:
                                break
                        if blocked:
                            continue
                        result = assess(replacements, positions)
                        if result:
                            _, events, geometries = result
                            for k in shifts:
                                if nodes[k]['kind'] != 'currency':
                                    occupied.pop(tuple(xy[k]), None)
                            xy.update(shifts)
                            obstacle_index = ObstacleIndex(xy)
                            hit_cache.clear()
                            for k in shifts:
                                if nodes[k]['kind'] != 'currency':
                                    occupied[tuple(xy[k])] = k
                            index.commit(replacements, events, geometries)
                            accepted_nodes += 1
                            improved = True
                            break
            if not improved:
                break
    result_paths = [p for rid in sorted(index.groups) for p in index.groups[rid]]
    after = index.metrics()
    assert score(after) <= score(before)
    hits = [{'route': list(route_key(p)), 'nodes': sorted(node_hits(p, obstacle_index))}
            for p in result_paths if node_hits(p, obstacle_index)]
    report = dict(objective=mode, crossing_weight=crossing_weight, bend_weight=bend_weight,
                  before=before, after=after, nonworsening_verified=True,
                  global_optimality_proven=False, method='time-bounded alternating route and grid-node local search',
                  elapsed_seconds=time.monotonic()-started, time_limit_seconds=time_limit,
                  timed_out=time.monotonic() >= deadline and time_limit > 0,
                  accepted_route_changes=accepted_routes, accepted_node_moves=accepted_nodes,
                  moved_node_ids=sorted(k for k in xy if xy[k] != initial_xy[k]),
                  search_enabled=bool(route_offset > 0 and time_limit > 0), node_moves_enabled=bool(node_moves),
                  evaluated_candidates=evaluated, passes=passes,
                  initial_node_hits=initial_node_hits, remaining_node_hit_count=sum(len(item['nodes']) for item in hits),
                  remaining_node_hits=hits, remaining_coincident_intervals=after['coincident_intervals'],
                  constraints_satisfied=not hits and after['coincident_intervals'] == 0,
                  metric_definition='All colors; unique reaction-pair intersections; common biological endpoints excluded; same-reaction shared corners counted once',
                  limitation='Finite candidates; fixed backbone, EX and structural red/orange nodes; existing node hits or overlaps may remain and are reported')
    return result_paths, xy, report
