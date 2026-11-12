"""Shared TCA circular connection geometry for SVG, Escher and SBML export."""
import math


def backbone_from_metadata(metadata):
    return next((metadata[k].get('major_backbones', {})
                 for k in ('sbml2map', 'sbml2escher', 'escher_grid') if k in metadata), {})


def curved_pairs(backbones):
    region = (backbones or {}).get('tca_geometry')
    if not region or region['shape'] != 'circle':
        return set()
    return {frozenset(('rxn:'+step['reaction_id'], 'met:'+mid))
            for route in backbones.get('routes', []) if route['id'] == 'tca'
            for step in route.get('steps', []) for mid in (step['source'], step['target'])}


def arc_controls(start, end, center):
    """Cubic Bezier tangential to the circle at both biological endpoints."""
    a = math.atan2(start[1]-center[1], start[0]-center[0])
    b = math.atan2(end[1]-center[1], end[0]-center[0])
    delta = (b-a+math.pi) % (2*math.pi)-math.pi
    radius = math.dist(start, center)
    if not math.isclose(radius, math.dist(end, center), abs_tol=1e-6):
        raise ValueError('TCA arc endpoints must lie on the same circle')
    k = 4/3*math.tan(delta/4)*radius
    return [[start[0]-k*math.sin(a), start[1]+k*math.cos(a)],
            [end[0]+k*math.sin(b), end[1]-k*math.cos(b)]]


def sample_curve(start, controls, end, intervals=32):
    points = [start, *controls, end]
    return [[sum(w*p[axis] for w,p in zip(((1-t)**3,3*(1-t)**2*t,3*(1-t)*t*t,t**3), points))
             for axis in (0,1)] for t in (i/intervals for i in range(intervals+1))]


def transform_point(point, transform):
    return [transform['scale']*(point[0]-transform['x_min'])+transform['padding'],
            transform['scale']*(transform['y_max']-point[1])+transform['padding']]


def escher_curve_controls(metadata, pair, start, end):
    """Allow curves only on supported TCA incidences in a circular layout."""
    backbones = backbone_from_metadata(metadata)
    if frozenset(pair) not in curved_pairs(backbones):
        return None
    region = backbones['tca_exclusion']['region']
    center = transform_point(region['center'], metadata['coordinate_transform'])
    return arc_controls(start, end, center)
