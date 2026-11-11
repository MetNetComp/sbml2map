"""Display-only rectangular detours for long adjacent highlighted segments."""
from collections import defaultdict
import math

EPS = 1e-7



def structural_node_priorities(nodes, major_reactions=(), major_metabolites=(), connector_pairs=()):
    """Use the structural red/orange node policy, independent of UI selections.

    Cofactor copies on red reactions remain gray unless an orange connector
    explicitly uses that copy. Red metabolites are backbone participants only.
    """
    major_reactions=set(major_reactions);major_metabolites=set(major_metabolites)
    orange_nodes={key for pair in (connector_pairs or ()) for key in pair}
    return {key:(2 if ((node['kind']=='reaction' and node['original_id'] in major_reactions) or
                       (node['kind']=='metabolite' and node['original_id'] in major_metabolites)) else
                 1 if key in orange_nodes else 0)
            for key,node in nodes.items()}


def _legs(points):
    for a, b in zip(points, points[1:]):
        if math.dist(a, b) <= EPS:
            continue
        if abs(a[1] - b[1]) < EPS:
            axis = 0
        elif abs(a[0] - b[0]) < EPS:
            axis = 1
        else:
            raise ValueError('Parallel routing requires orthogonal segments')
        yield dict(a=tuple(a), b=tuple(b), axis=axis, level=a[1-axis],
                   lo=min(a[axis], b[axis]), hi=max(a[axis], b[axis]))


def _overlap(a, b):
    return min(a['hi'], b['hi']) - max(a['lo'], b['lo'])


def _distance(point, leg):
    axis = leg['axis']
    return math.hypot(max(leg['lo']-point[axis], 0, point[axis]-leg['hi']),
                      point[1-axis]-leg['level'])


def _detour(leg, offset, grid):
    a, b = leg['a'], leg['b']
    if abs(offset) < EPS:
        return [a, b]
    axis = leg['axis']; side = 1-axis
    stub = min(.08*grid, (leg['hi']-leg['lo'])/4)
    p, q = list(a), list(b)
    sign = 1 if b[axis] > a[axis] else -1
    p[axis] += sign*stub; q[axis] -= sign*stub
    u, v = list(p), list(q)
    u[side] += offset; v[side] += offset
    return [a, tuple(p), tuple(u), tuple(v), tuple(q), b]


def separate_parallel(paths, grid=1.0, obstacles=(), clearance=.3, min_overlap=1.0):
    """Mutate path points, keeping endpoints and fixed backbone segments intact.

    Paths provide rid, priority (2 red / 1 orange / 0 other), points and optional
    backbone. Same-reaction, same-color trunks share an offset only when they overlap
    with positive length; touching at an endpoint does not bundle them.
    Obstacles are biological node coordinates; crossings of other lines are
    scored when choosing a side. Unplaceable groups remain unchanged and are
    recorded. Geometry, never fluxes or chemistry, is changed.
    """
    gap, threshold = clearance*grid, min_overlap*grid
    obstacles = [tuple(p) for p in obstacles]
    entries=[]; buckets=defaultdict(list)
    for index, path in enumerate(paths):
        for part, leg in enumerate(_legs(path['points'])):
            leg.update(path=index, part=part, rid=path['rid'], priority=path['priority'],
                       fixed=bool(path.get('backbone') or path.get('routing_locked')), offset=0.0)
            entries.append(leg)
            if leg['priority'] > 0:
                buckets[(leg['rid'], leg['priority'], leg['axis'], round(leg['level']/grid, 8))].append(leg)
    groups=[]
    for key, legs in sorted(buckets.items()):
        group=[]; end=-math.inf
        for leg in sorted(legs, key=lambda x:(x['lo'],x['hi'],x['path'],x['part'])):
            # Sharing only an endpoint must not bundle a collision-free branch.
            if group and leg['lo'] >= end-EPS:
                groups.append(group); group=[]
            group.append(leg); end=max(end,leg['hi'])
        if group: groups.append(group)
    groups.sort(key=lambda g:(not any(x['fixed'] for x in g), -g[0]['priority'],
                              g[0]['rid'],g[0]['axis'],g[0]['level'],min(x['lo'] for x in g)))
    reserved=[]; unresolved=[]; shifted=0
    def conflicts(group, offset):
        return [other for other in reserved if any(
            leg['rid'] != other['rid'] and leg['axis'] == other['axis'] and
            _overlap(leg,other) >= threshold-EPS and
            abs(leg['level']+offset-other['level']) < gap-EPS for leg in group)]
    def collision(group, offset):
        endpoints={p for leg in group for p in (leg['a'],leg['b'])}
        segments=[s for leg in group for s in _legs(_detour(leg,offset,grid))]
        return any(_distance(point, s) < .06*grid-EPS
                   for point in obstacles if point not in endpoints for s in segments)
    def crossings(group, offset):
        count=0
        for leg in group:
            for part in _legs(_detour(leg,offset,grid)):
                for other in entries:
                    if other['rid']==leg['rid'] or part['axis']==other['axis']:continue
                    h,v=(part,other) if part['axis']==0 else (other,part)
                    if h['lo']+EPS < v['level'] < h['hi']-EPS and v['lo']+EPS < h['level'] < v['hi']-EPS:
                        count+=1
        return count
    for group in groups:
        clashes=conflicts(group,0)
        offset=0.0
        if clashes:
            candidates=set()
            if not any(leg['fixed'] for leg in group):
                for leg in group:
                    for other in reserved:
                        if leg['axis']==other['axis'] and _overlap(leg,other)>=threshold-EPS:
                            for side in (-1,1):
                                candidate=other['level']+side*gap-leg['level']
                                if EPS < abs(candidate) <= 1.2*grid+EPS:candidates.add(candidate)
            viable=[value for value in candidates if not conflicts(group,value) and not collision(group,value)]
            if viable:
                # Nearest valid lane first; choose the less crowded side on ties.
                offset=min(viable,key=lambda value:(round(abs(value)/grid,8),crossings(group,value),value<0))
                shifted+=len(group)
            else:
                unresolved.append(dict(reaction_id=group[0]['rid'],axis=group[0]['axis'],
                                       level=group[0]['level'],
                                       reason='fixed backbone' if any(x['fixed'] for x in group) else 'no free lane within 1.2 grids',
                                       conflicting_reactions=sorted({x['rid'] for x in clashes})))
        for leg in group:
            leg['offset']=offset
            # Reserve the long shaft. The original endpoints remain connected by short stubs.
            reserved.append(dict(leg,level=leg['level']+offset))
    by_path=defaultdict(list)
    for leg in entries:by_path[leg['path']].append(leg)
    for index,path in enumerate(paths):
        if not any(leg['offset'] for leg in by_path[index]):continue
        points=[]
        for leg in by_path[index]:
            detour=_detour(leg,leg['offset'],grid)
            points.extend(detour if not points else detour[1:])
        path['points']=points
    return dict(trunk_grouping='positive_length_overlap',clearance_grid=clearance,min_parallel_length_grid=min_overlap,
                shifted_segments=shifted,unresolved=unresolved,
                max_offset_grid=max((abs(x['offset'])/grid for x in entries),default=0),
                fixed_endpoints=True,shared_same_reaction_trunks=True)


def validate_route_offset(value):
    """A display displacement in grid units; zero disables node avoidance."""
    if not math.isfinite(value) or value < 0:
        raise ValueError('route-offset must be finite and nonnegative')
    return value


def _u_detour(leg, offset):
    """Exactly two right-angle bends at the original segment endpoints."""
    a,b=leg['a'],leg['b']
    if abs(offset)<EPS:return [a,b]
    side=1-leg['axis'];p=list(a);q=list(b)
    p[side]+=offset;q[side]+=offset
    return [a,tuple(p),tuple(q),b]



def _node_detour(leg, offset, grid, inset_start=False, inset_end=False):
    """Step along an endpoint before taking two right-angle turns, if needed."""
    if not inset_start and not inset_end:return _u_detour(leg,offset)
    axis=leg['axis'];side=1-axis;a,b=leg['a'],leg['b']
    setback=min(.15*grid,(leg['hi']-leg['lo'])/4)
    sign=1 if b[axis]>a[axis] else -1
    start,end=list(a),list(b)
    if inset_start:start[axis]+=sign*setback
    if inset_end:end[axis]-=sign*setback
    shifted_start,shifted_end=list(start),list(end)
    shifted_start[side]+=offset;shifted_end[side]+=offset
    return ([a]+([tuple(start)] if inset_start else [])+
            [tuple(shifted_start),tuple(shifted_end)]+
            ([tuple(end)] if inset_end else [])+[b])


def avoid_nodes(paths, grid=1.0, obstacles=(), route_offset=.3):
    """Offset a whole red/orange segment when it crosses an unrelated node.

    Biological endpoints stay fixed, including on primary backbones. Bends are
    display geometry only. Same-reaction overlapping trunks move together. Only
    the requested positive/negative displacement is tried, never a larger one.
    If a direct U is blocked, a 0.15-grid step is tried at either or both ends.
    """
    validate_route_offset(route_offset)
    report=dict(trunk_grouping='positive_length_overlap',route_offset_grid=route_offset,enabled=route_offset>0,
                shifted_segments=0,unresolved=[],shape='rectangle with optional endpoint steps',
                endpoint_step_grid=.15,endpoint_step_fallback_segments=0,detours=[],
                fixed_endpoints=True,shared_same_reaction_trunks=True,
                backbone_shifted_segments=0,max_offset_grid=0,backbone_max_offset_grid=0)
    if route_offset==0:return report
    amount=route_offset*grid;radius=.12*grid
    if not math.isfinite(amount):
        raise ValueError('route-offset is too large for the drawing scale')
    obstacles=sorted(set(tuple(p) for p in obstacles))
    entries=[];buckets=defaultdict(list)
    for index,path in enumerate(paths):
        endpoints={tuple(path['points'][0]),tuple(path['points'][-1])}
        for part,leg in enumerate(_legs(path['points'])):
            leg.update(path=index,part=part,rid=path['rid'],priority=path['priority'],
                       endpoints=endpoints,backbone=bool(path.get('backbone')),routing_locked=bool(path.get('routing_locked')),points=[leg['a'],leg['b']])
            entries.append(leg)
            if leg['priority']>0:
                buckets[(leg['rid'],leg['priority'],leg['axis'],round(leg['level']/grid,8))].append(leg)
    groups=[]
    for key,legs in sorted(buckets.items()):
        group=[];end=-math.inf
        for leg in sorted(legs,key=lambda x:(x['lo'],x['hi'],x['path'],x['part'])):
            # Sharing only an endpoint must not bundle a collision-free branch.
            if group and leg['lo']>=end-EPS:
                groups.append(group);group=[]
            group.append(leg);end=max(end,leg['hi'])
        if group:groups.append(group)
    groups.sort(key=lambda g:(-g[0]['priority'],g[0]['rid'],g[0]['axis'],g[0]['level'],min(x['lo'] for x in g)))
    def hits(leg,points):
        return [node for node in obstacles if node not in leg['endpoints'] and
                any(_distance(node,part)<radius-EPS for part in _legs(points))]
    def retained_endpoint_stem(leg,segment,current,length):
        # A step may keep an already coincident short stem at a biological
        # endpoint. It must not add overlap or preserve a long shared path.
        if segment['axis']!=leg['axis'] or abs(segment['level']-leg['level'])>EPS:return False
        if length>.15*grid+EPS or length>_overlap(leg,current)+EPS:return False
        return any(_distance(point,segment)<EPS and _distance(point,current)<EPS and
                   (point==leg['a'] or point==leg['b']) for point in leg['endpoints'])
    def interference(group,shapes):
        crossings=0;retained=0;gray_penalty=0
        for leg,points in zip(group,shapes):
            for segment in _legs(points):
                for other in entries:
                    if leg['rid']==other['rid']:continue
                    if other['priority']==0:
                        for current in _legs(other['points']):
                            if segment['axis']==current['axis'] and abs(segment['level']-current['level'])<EPS:
                                gray_penalty+=max(0,_overlap(segment,current))/grid
                        continue
                    for current in _legs(other['points']):
                        if segment['axis']==current['axis']:
                            length=_overlap(segment,current);distance=abs(segment['level']-current['level'])
                            # Avoid new coincident paths and long close red/orange lanes.
                            if length>EPS and distance<EPS:
                                if retained_endpoint_stem(leg,segment,current,length):retained+=1
                                else:return None
                            if other['priority']>0 and length>=grid-EPS and distance<.3*grid-EPS:return None
                        else:
                            h,v=(segment,current) if segment['axis']==0 else (current,segment)
                            if h['lo']+EPS<v['level']<h['hi']-EPS and v['lo']+EPS<h['level']<v['hi']-EPS:
                                crossings+=1
        return crossings,retained,gray_penalty
    for group in groups:
        collisions=[(leg,hits(leg,leg['points'])) for leg in group]
        if not any(nodes for _,nodes in collisions):continue
        if any(leg['routing_locked'] for leg in group):
            report['unresolved'].append(dict(reaction_id=group[0]['rid'],reason='preallocated red lane requires protected-node review'))
            continue
        candidates=[]
        # Preserve the simple U when possible. Add steps at one end, then both.
        for variants in [[(False,False)],[(True,False),(False,True)],[(True,True)]]:
            for mode,(inset_start,inset_end) in enumerate(variants):
                for offset in (amount,-amount):
                    shapes=[_node_detour(leg,offset,grid,inset_start,inset_end) for leg in group]
                    if any(hits(leg,points) for leg,points in zip(group,shapes)):continue
                    score=interference(group,shapes)
                    if score is not None:
                        candidates.append((score,mode,offset<0,offset,inset_start,inset_end,shapes))
            if candidates:break
        if candidates:
            score,_,_,offset,inset_start,inset_end,shapes=min(candidates,key=lambda c:c[:3])
            for leg,points in zip(group,shapes):leg['points']=points
            if inset_start or inset_end:report['endpoint_step_fallback_segments']+=len(group)
            report['detours'].append(dict(reaction_id=group[0]['rid'],axis=group[0]['axis'],
                offset_grid=offset/grid,inset_start=inset_start,inset_end=inset_end,
                segments=len(group),existing_endpoint_stem_overlaps=score[1]))
            report['shifted_segments']+=len(group)
            backbone_count=sum(leg['backbone'] for leg in group)
            report['backbone_shifted_segments']+=backbone_count
            report['max_offset_grid']=route_offset
            if backbone_count:report['backbone_max_offset_grid']=route_offset
        else:
            report['unresolved'].append(dict(reaction_id=group[0]['rid'],axis=group[0]['axis'],
                level=group[0]['level'],reason='requested U-shaped lanes and endpoint-step fallbacks are obstructed',
                nodes=[list(node) for node in sorted({node for _,nodes in collisions for node in nodes})]))
    by_path=defaultdict(list)
    for leg in entries:by_path[leg['path']].append(leg)
    for index,path in enumerate(paths):
        points=[]
        for leg in by_path[index]:points.extend(leg['points'] if not points else leg['points'][1:])
        if points:path['points']=points
    return report
