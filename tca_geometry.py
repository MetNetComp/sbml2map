"""TCA perimeter geometry and renderer-independent node exclusion checks."""
import math
from collections import defaultdict
from functools import lru_cache
import numpy as np

EPS = 1e-7
NODE_CLEARANCE = .16
TCA_GRID_STEP = .5


@lru_cache(maxsize=64)
def lattice_circle(required_slots, minimum_radius=2.):
    """Smallest circle with enough integer/half-integer centered lattice points.

    Search doubled integer coordinates, grouped by parity and squared radius.
    Translation then puts every selected point on the ordinary integer grid.
    """
    limit = max(8, required_slots)
    while True:
        groups = defaultdict(list)
        for x in range(-limit, limit+1):
            for y in range(-limit, limit+1):
                radius2 = x*x+y*y
                if 4*minimum_radius**2-EPS <= radius2 <= limit*limit:
                    groups[(radius2, x % 2, y % 2)].append((x/2, y/2))
        candidates = [(key, points) for key, points in groups.items() if len(points) >= required_slots]
        if candidates:
            (radius2, _, _), points = min(candidates, key=lambda item:item[0])
            # Begin near the upper-left entry; walk clockwise in y-up space.
            points.sort(key=lambda p:(3*math.pi/4-math.atan2(p[1],p[0])) % (2*math.pi))
            return tuple(points), math.sqrt(radius2)/2
        limit *= 2


def configure_plan(fixed, plan, shape='circle', interior_policy='forbid'):
    if interior_policy not in ('allow','forbid'):
        raise ValueError('tca-interior must be allow or forbid')
    if shape not in ('circle', 'square'):
        raise ValueError('tca-shape must be circle or square')
    route = next((r for r in plan['routes'] if r['id'] == 'tca'), None)
    square = plan.get('square')
    if route is None or square is None:
        plan['tca_geometry'] = None
        return fixed, plan
    offsets = route['node_offsets']
    # Missing cycle steps keep their original angular gaps; no closure is added.
    count = max(1, 2*len(route.get('pattern', ())[:-1])) if not route.get('closed', True) and route.get('pattern') else len(offsets)
    side = square['side_length']
    left, top = square['left'], square['top_offset']
    node_ids = list(route['displayed_node_ids'])
    region = dict(shape=shape, node_ids=node_ids, clearance_grid=NODE_CLEARANCE, interior_policy=interior_policy,
                  orientation='clockwise', closed=route.get('closed', True))
    if shape == 'circle':
        # A genuine circle through integer grid intersections, rather than
        # rounded samples of a smaller circle. Keep enough slots for absent steps.
        lattice, radius = lattice_circle(max(count, 4*side), max(2., side/2))
        # Half-step auxiliary grid halves the diameter without losing slots.
        lattice = tuple((x*TCA_GRID_STEP,y*TCA_GRID_STEP) for x,y in lattice)
        radius *= TCA_GRID_STEP
        xmin = min(x for x, y in lattice)
        ymax = max(y for x, y in lattice)
        shift = [left-xmin, top-ymax]
        center = [shift[0], shift[1]]
        for key, (x, y) in list(offsets.items()):
            if abs(y-top) < EPS:
                progress = x-left
            elif abs(x-left-side) < EPS:
                progress = side+top-y
            elif abs(y-top+side) < EPS:
                progress = 2*side+left+side-x
            else:
                progress = 3*side+y-(top-side)
            slot = round(progress*len(lattice)/(4*side)) % len(lattice)
            position = [round((lattice[slot][i]+shift[i])/TCA_GRID_STEP)*TCA_GRID_STEP for i in (0, 1)]
            offsets[key] = position
            if key in fixed:
                fixed[key] = position
                plan['fixed_node_offsets'][key] = position
        region.update(center_offset=center, radius=radius,
                      lattice_point_count=len(lattice), required_slots=max(count, 4*side),
                      all_tca_nodes_on_grid=True, node_grid_step=TCA_GRID_STEP,
                      diameter_scale=.5, connection_style='circular_bezier')
        plan['circle'] = dict(center_offset=center, radius=radius, orientation='clockwise',
                              lattice_point_count=len(lattice))
        plan['square'] = None
    else:
        region.update(left=left, top_offset=top, side_length=side)
        plan['circle'] = None
    route['shape'] = shape
    plan['tca_geometry'] = region
    return fixed, plan


def absolute_region(region, rows):
    if region is None:
        return None
    result = dict(region)
    if region['shape'] == 'circle':
        x, y = region['center_offset']
        result['center'] = [x, rows-1+y]
    else:
        result['top'] = rows-1+region['top_offset']
    return result


def allowed_points(points, region, clearance=None):
    points = np.asarray(points, dtype=float).reshape((-1, 2))
    if region is None:
        return np.ones(len(points), dtype=bool)
    margin = region.get('clearance_grid', NODE_CLEARANCE) if clearance is None else clearance
    if region['shape'] == 'circle':
        distance=np.linalg.norm(points-np.asarray(region['center']), axis=1)
        if region.get('interior_policy')=='allow':
            return np.abs(distance-region['radius']) >= margin-EPS
        return distance >= region['radius']+margin-EPS
    left, top, side = region['left'], region['top'], region['side_length']
    outside = ((points[:,0] <= left-margin+EPS) | (points[:,0] >= left+side+margin-EPS) |
            (points[:,1] <= top-side-margin+EPS) | (points[:,1] >= top+margin-EPS))
    if region.get('interior_policy')=='allow':
        inside=((points[:,0]>=left+margin-EPS)&(points[:,0]<=left+side-margin+EPS)&
                (points[:,1]>=top-side+margin-EPS)&(points[:,1]<=top-margin+EPS))
        return outside|inside
    return outside


def on_perimeter(point, region):
    x, y = point
    if region['shape'] == 'circle':
        return abs(math.dist(point, region['center'])-region['radius']) < EPS
    left, top, side = region['left'], region['top'], region['side_length']
    return (left-EPS <= x <= left+side+EPS and top-side-EPS <= y <= top+EPS and
            min(abs(x-left),abs(x-left-side),abs(y-top),abs(y-top+side)) < EPS)


def validate_exclusion(nodes, xy, region):
    if region is None:
        return dict(enabled=False, interior_node_count=0, clearance_violation_count=0, verified=True)
    ring = set(region['node_ids'])
    permit=region.get('interior_policy')=='allow'
    outside_region=dict(region,interior_policy='forbid')
    interior, clearance, misplaced = [], [], []
    for key, point in xy.items():
        if key in ring:
            if not on_perimeter(point, region):
                misplaced.append(key)
        else:
            if not allowed_points([point], outside_region, 0)[0]:
                interior.append(key)
            if not allowed_points([point], region)[0]:
                clearance.append(key)
    report = dict(enabled=True, shape=region['shape'], region=region, interior_policy='allow' if permit else 'forbid',
                  interior_placement_allowed=permit,
                  interior_node_count=len(interior), interior_node_ids=interior,
                  clearance_violation_count=len(clearance), clearance_violation_node_ids=clearance,
                  misplaced_tca_node_ids=misplaced, checked_nodes=len(xy),
                  applies_to='all metabolites, reactions and reaction-local cofactor copies',
                  verified=(permit or not interior) and not clearance and not misplaced)
    if not report['verified']:
        raise RuntimeError('TCA perimeter/interior exclusion violated: '+str(interior+clearance+misplaced))
    return report
