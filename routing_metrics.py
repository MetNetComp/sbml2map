"""Geometry-based routing metrics, independent of renderer and segment splitting."""
from collections import defaultdict
import math
import numpy as np
from parallel_routing import EPS, _legs
from routing_simplification import clean


def point_key(point):
    return tuple(round(float(v), 7) for v in point)


def route_key(path):
    return (path['rid'], path.get('source', ''), path.get('target', ''), path.get('mid', ''))


def geometry(paths):
    """Union collinear same-reaction trunks; count a shared corner once."""
    buckets = defaultdict(list)
    corners = set()
    endpoints = defaultdict(set)
    for path in paths:
        points = clean(path['points'])
        for point, name in ((points[0], path.get('source')), (points[-1], path.get('target'))):
            endpoints[point_key(point)].add(name if name is not None else point_key(point))
        for a, b, c in zip(points, points[1:], points[2:]):
            if (abs(a[0]-b[0]) < EPS) != (abs(b[0]-c[0]) < EPS):
                corners.add(point_key(b))
        for leg in _legs(points):
            buckets[(leg['axis'], round(float(leg['level']), 7))].append((leg['lo'], leg['hi']))
    segments = []
    for (axis, level), intervals in sorted(buckets.items()):
        lo, hi = None, None
        for start, end in sorted(intervals):
            if lo is not None and start > hi + EPS:
                segments.append((axis, level, lo, hi))
                lo, hi = start, end
            elif lo is None:
                lo, hi = start, end
            else:
                hi = max(hi, end)
        if lo is not None:
            segments.append((axis, level, lo, hi))
    return dict(segments=segments, array=np.asarray(segments, dtype=float), corners=corners, endpoints=endpoints,
                length=sum(hi-lo for _, _, lo, hi in segments))


def interactions(left, right):
    """Count inclusive intersections, except genuine shared biological endpoints.

    A crossing at a corner or at an unrelated line endpoint still counts. Each
    reaction-pair/coordinate counts once, independent of polyline segmentation.
    """
    crossings, overlaps = set(), set()
    if not left['segments'] or not right['segments']:
        return crossings, overlaps
    other = right['array']
    axis, level, lo, hi = other.T
    for a, c, start, end in left['segments']:
        mask = ((axis != a) & (lo <= c+EPS) & (hi >= c-EPS) &
                (level >= start-EPS) & (level <= end+EPS))
        for v in level[mask]:
            point = point_key((v, c) if a == 0 else (c, v))
            if left['endpoints'].get(point, set()) & right['endpoints'].get(point, set()):
                continue
            crossings.add(point)
        mask = ((axis == a) & (np.abs(level-c) < EPS) &
                (np.minimum(hi, end)-np.maximum(lo, start) > EPS))
        for l, h in zip(lo[mask], hi[mask]):
            overlaps.add((a, round(c, 7), round(max(l, start), 7), round(min(h, end), 7)))
    return crossings, overlaps


class MetricIndex:
    """Reaction geometry and pair events cached for local route changes."""
    def __init__(self, paths):
        groups = defaultdict(list)
        for path in paths:
            groups[path['rid']].append(path)
        self.groups = dict(groups)
        self.geometries = {rid: geometry(group) for rid, group in groups.items()}
        self.pairs = {}
        keys = sorted(groups)
        for i, a in enumerate(keys):
            for b in keys[i+1:]:
                events = interactions(self.geometries[a], self.geometries[b])
                if any(events):
                    self.pairs[(a, b)] = events

    def metrics(self):
        crossings = sum(len(v[0]) for v in self.pairs.values())
        bends = sum(len(g['corners']) for g in self.geometries.values())
        return dict(crossings=crossings, bends=bends, crossings_plus_bends=crossings+bends,
                    total_length=float(sum(g['length'] for g in self.geometries.values())),
                    coincident_intervals=sum(len(v[1]) for v in self.pairs.values()))

    def evaluate(self, replacements):
        """Return metrics and events for replacements keyed by reaction ID."""
        changed = {rid: geometry(paths) for rid, paths in replacements.items()}
        events = {}
        for a in sorted(self.geometries):
            for b in sorted(changed):
                if a == b or (a in changed and a > b):
                    continue
                key = tuple(sorted((a, b)))
                events[key] = interactions(changed.get(a, self.geometries[a]), changed[b])
        result = self.metrics()
        for key, value in events.items():
            previous = self.pairs.get(key, (set(), set()))
            result['crossings'] += len(value[0])-len(previous[0])
            result['coincident_intervals'] += len(value[1])-len(previous[1])
        result['bends'] += sum(len(g['corners'])-len(self.geometries[rid]['corners']) for rid, g in changed.items())
        result['total_length'] += sum(g['length']-self.geometries[rid]['length'] for rid, g in changed.items())
        result['crossings_plus_bends'] = result['crossings']+result['bends']
        return result, events, changed

    def commit(self, replacements, events, changed):
        self.groups.update(replacements)
        self.geometries.update(changed)
        for key, value in events.items():
            if any(value):
                self.pairs[key] = value
            else:
                self.pairs.pop(key, None)


def measure(paths):
    return MetricIndex(paths).metrics()


def objective(metrics, mode='crossings-first', crossing_weight=10.0, bend_weight=1.0):
    if mode == 'crossings-first':
        primary = (metrics['crossings'], metrics['bends'])
    elif mode == 'sum':
        primary = (metrics['crossings_plus_bends'],)
    elif mode == 'weighted':
        primary = (crossing_weight*metrics['crossings']+bend_weight*metrics['bends'],)
    else:
        raise ValueError('Unknown routing objective: '+mode)
    return (*primary, round(metrics['total_length'], 7))


def validate_options(mode, crossing_weight, bend_weight, time_limit):
    if mode not in ('crossings-first', 'sum', 'weighted'):
        raise ValueError('Invalid routing objective')
    if any(not math.isfinite(v) or v < 0 for v in (crossing_weight, bend_weight, time_limit)):
        raise ValueError('Routing weights and time limit must be finite and nonnegative')
    if crossing_weight == bend_weight == 0:
        raise ValueError('At least one routing weight must be positive')
