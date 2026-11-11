"""TCA perimeter geometry and renderer-independent node exclusion checks."""
import math
import numpy as np

EPS = 1e-7
NODE_CLEARANCE = .16


def configure_plan(fixed, plan, shape='circle'):
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
    region = dict(shape=shape, node_ids=node_ids, clearance_grid=NODE_CLEARANCE,
                  orientation='clockwise', closed=route.get('closed', True))
    if shape == 'circle':
        radius = max(2., side/2, count/(2*math.pi))
        center = (left+radius, top-radius)
        # Use square perimeter progress to preserve missing-step gaps and order.
        for key, (x, y) in list(offsets.items()):
            if abs(y-top) < EPS:
                progress = x-left
            elif abs(x-left-side) < EPS:
                progress = side+top-y
            elif abs(y-top+side) < EPS:
                progress = 2*side+left+side-x
            else:
                progress = 3*side+y-(top-side)
            angle = 3*math.pi/4-2*math.pi*progress/(4*side)
            position = [center[0]+radius*math.cos(angle), center[1]+radius*math.sin(angle)]
            offsets[key] = position
            if key in fixed:
                fixed[key] = position
                plan['fixed_node_offsets'][key] = position
        region.update(center_offset=list(center), radius=radius)
        plan['circle'] = dict(center_offset=list(center), radius=radius, orientation='clockwise')
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
        return np.linalg.norm(points-np.asarray(region['center']), axis=1) >= region['radius']+margin-EPS
    left, top, side = region['left'], region['top'], region['side_length']
    return ((points[:,0] <= left-margin+EPS) | (points[:,0] >= left+side+margin-EPS) |
            (points[:,1] <= top-side-margin+EPS) | (points[:,1] >= top+margin-EPS))


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
    interior, clearance, misplaced = [], [], []
    for key, point in xy.items():
        if key in ring:
            if not on_perimeter(point, region):
                misplaced.append(key)
        else:
            if not allowed_points([point], region, 0)[0]:
                interior.append(key)
            if not allowed_points([point], region)[0]:
                clearance.append(key)
    report = dict(enabled=True, shape=region['shape'], region=region,
                  interior_node_count=len(interior), interior_node_ids=interior,
                  clearance_violation_count=len(clearance), clearance_violation_node_ids=clearance,
                  misplaced_tca_node_ids=misplaced, checked_nodes=len(xy),
                  applies_to='all metabolites, reactions and reaction-local cofactor copies',
                  verified=not interior and not clearance and not misplaced)
    if not report['verified']:
        raise RuntimeError('TCA perimeter/interior exclusion violated: '+str(interior+clearance+misplaced))
    return report
