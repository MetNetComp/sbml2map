"""Separate coincident orthogonal incidences without moving biological nodes."""
from collections import defaultdict
import hashlib
import json
import math

from routing_simplification import POLICY,simplify_important,yield_gray

from red_route_layout import plan_from_metadata,spacing_report

from parallel_routing import separate_parallel, avoid_nodes, validate_route_offset, structural_node_priorities


def _simplify(points):
    out=[]
    for point in points:
        point=tuple(point)
        if out and all(abs(a-b)<1e-8 for a,b in zip(out[-1],point)):
            continue
        while len(out)>1:
            a,b=out[-2:]
            aligned=(abs(a[0]-b[0])<1e-8 and abs(b[0]-point[0])<1e-8 or
                     abs(a[1]-b[1])<1e-8 and abs(b[1]-point[1])<1e-8)
            forward=sum((b[k]-a[k])*(point[k]-b[k]) for k in (0,1))>=0
            if not (aligned and forward):break
            out.pop()
        out.append(point)
    return out


def separate_routes(data,major_reactions=(),connector_pairs=None,enabled=True,route_offset=None):
    """Use small deterministic lanes, red before orange before gray.

    A reaction's same-color shared trunks stay bundled. Short terminal stubs
    join offset paths to the original nodes. Disabling restores plain L routes
    for simulation views, with no structural priority.
    """
    header,body=data;meta=header['metadata']
    route_offset=meta.get('route_offset_grid',.3) if route_offset is None else route_offset
    validate_route_offset(route_offset)
    meta['route_offset_grid']=route_offset
    if meta.get('routing')!='orthogonal':return data
    major=set(major_reactions);orange=set(connector_pairs or ())
    sources=meta['source_node_ids'];old_nodes=body['nodes']
    base_id=meta.get('overlap_routing',{}).get('base_map_id',header['map_id'])
    backbone=set();major_metabolites=set()
    for namespace in ('sbml2map','sbml2escher','escher_grid'):
        for route in meta.get(namespace,{}).get('major_backbones',{}).get('routes',[]):
            major_metabolites.update(route.get('metabolite_ids',()))
            for step in route.get('steps',[]):
                reaction='rxn:'+step['reaction_id']
                backbone.add(frozenset((reaction,'met:'+step['source'])))
                backbone.add(frozenset((reaction,'met:'+step['target'])))
    node_records={key:dict(kind='reaction' if old_nodes[nid]['node_type']=='midmarker' else
        'metabolite' if old_nodes[nid].get('node_is_primary') else 'currency',
        original_id=key.removeprefix('rxn:') if old_nodes[nid]['node_type']=='midmarker' else old_nodes[nid]['bigg_id'])
        for nid,key in sources.items()}
    node_priorities=structural_node_priorities(node_records,major,major_metabolites,orange)
    obstacle_keys=sorted(key for key,rank in node_priorities.items() if rank>0)
    obstacles=[(old_nodes[nid]['x'],old_nodes[nid]['y']) for nid,key in sources.items() if key in obstacle_keys]
    meta['structural_node_priorities']=node_priorities
    routes=[]
    for rid,reaction in body['reactions'].items():
        marker=meta['reactions'][reaction['bigg_id']]['marker_node_id']
        targets={s[k] for s in reaction['segments'].values() for k in ('from_node_id','to_node_id')
                 if s[k] in sources and s[k]!=marker}
        details=meta['reactions'][reaction['bigg_id']]
        for target in sorted(targets,key=lambda n:sources[n]):
            mid=old_nodes[target]['bigg_id'];coefficient=details['original_stoichiometry'][mid]
            pair=[sources[target],sources[marker]] if coefficient<0 else [sources[marker],sources[target]]
            priority=2 if reaction['bigg_id'] in major else 1 if frozenset(pair) in orange else 0
            a=tuple(old_nodes[marker][k] for k in ('x','y'));b=tuple(old_nodes[target][k] for k in ('x','y'))
            points=[a]
            if abs(a[0]-b[0])>1e-9 and abs(a[1]-b[1])>1e-9:
                points.append((a[0],b[1]) if coefficient*details['coefficient_multiplier']<0 else (b[0],a[1]))
            points.append(b)
            routes.append(dict(reaction=rid,rid=reaction['bigg_id'],mid=mid,marker=marker,target=target,
                               pair=pair,priority=priority if enabled else 0,backbone=frozenset(pair) in backbone,
                               points=points,legs=[]))
    step=max(10.0,meta['coordinate_transform']['scale']*.06)
    legs=[]
    for route in routes:
        for a,b in zip(route['points'],route['points'][1:]):
            axis=0 if abs(a[1]-b[1])<1e-8 else 1
            leg=dict(route=route,axis=axis,level=a[1-axis],lo=min(a[axis],b[axis]),hi=max(a[axis],b[axis]),offset=0.0)
            route['legs'].append(leg);legs.append(leg)
    reservations=defaultdict(list)
    def conflicts(leg,level):
        bucket=math.floor(level/step);out=[]
        for i in (bucket-1,bucket,bucket+1):
            for other in reservations[(leg['axis'],i)]:
                if (other['route']['rid']==leg['route']['rid'] and
                        other['route']['priority']==leg['route']['priority']):continue
                length=min(leg['hi'],other['hi'])-max(leg['lo'],other['lo'])
                if length>step and abs(level-other['level']-other['offset'])<(step*.4 if leg['route']['priority']==other['route']['priority']==0 else step)-1e-7:
                    out.append((other,length))
        return out
    ordered=sorted(legs,key=lambda l:(-l['route']['priority'],not l['route']['backbone'],
                                      l['route']['rid'],l['route']['mid'],l['axis']))
    if enabled:
        for leg in ordered:
            clashes=conflicts(leg,leg['level'])
            if clashes and leg['hi']-leg['lo']>3*step:
                start=4 if leg['route']['priority']==0 and any(o['route']['priority']>0 for o,_ in clashes) else 1
                unit=step/2 if leg['route']['priority']==0 else step
                for lane in range(start,2*len(ordered)+2):
                    candidates=(lane*unit,-lane*unit)
                    chosen=next((offset for offset in candidates if not conflicts(leg,leg['level']+offset)),None)
                    if chosen is not None:
                        leg['offset']=chosen;break
            reservations[(leg['axis'],math.floor((leg['level']+leg['offset'])/step))].append(leg)
    for route in routes:
        if not any(l['offset'] for l in route['legs']):continue
        points=route['points'];first,last=route['legs'][0],route['legs'][-1]
        start,end=points[0],points[-1]
        port=max(28.0,2*step)
        a=list(start);z=list(end)
        first_axis,last_axis=first['axis'],last['axis']
        a[first_axis]+=math.copysign(min(port,(first['hi']-first['lo'])/3),points[1][first_axis]-start[first_axis])
        z[last_axis]-=math.copysign(min(port,(last['hi']-last['lo'])/3),end[last_axis]-points[-2][last_axis])
        shifted_a=list(a);shifted_a[1-first_axis]+=first['offset']
        shifted_z=list(z);shifted_z[1-last_axis]+=last['offset']
        joints=[]
        for previous,following in zip(route['legs'],route['legs'][1:]):
            horizontal=previous if previous['axis']==0 else following
            vertical=previous if previous['axis']==1 else following
            joints.append((vertical['level']+vertical['offset'],horizontal['level']+horizontal['offset']))
        route['points']=_simplify([start,a,shifted_a,*joints,shifted_z,z,end])
    red_plan=plan_from_metadata(meta) if enabled else None
    if red_plan and red_plan.get('enabled') and red_plan.get('verified'):
        planned={(r['rid'],r['mid']):r for r in red_plan['routes']}
        transform=meta['coordinate_transform'];scale=transform['scale']
        for route in routes:
            reserved=planned.get((route['rid'],route['mid']))
            if reserved:
                points=[(scale*(x-transform['x_min'])+transform['padding'],
                         scale*(transform['y_max']-y)+transform['padding']) for x,y in reserved['points']]
                if reserved['source']!=sources[route['marker']]:points.reverse()
                route.update(points=points,routing_locked=True)
    parallel_report = separate_parallel(
        routes, grid=meta['coordinate_transform']['scale'],
        obstacles=obstacles) if enabled else None
    node_report = avoid_nodes(routes, grid=meta['coordinate_transform']['scale'],
        obstacles=obstacles,
        route_offset=route_offset) if enabled else None
    requested_red_spacing=red_plan['requested_spacing_grid'] if red_plan and red_plan.get('enabled') else .3
    simplify_report=simplify_important(routes,grid=meta['coordinate_transform']['scale'],obstacles=obstacles,
        route_offset=route_offset,red_spacing=requested_red_spacing) if enabled else None
    gray_report=yield_gray(routes,grid=meta['coordinate_transform']['scale'],obstacles=obstacles,
        route_offset=route_offset) if enabled else None
    red_check=spacing_report(routes,red_plan['requested_spacing_grid'],meta['coordinate_transform']['scale']) if red_plan else None
    if red_plan and red_plan.get('enabled') and red_plan.get('verified') and not red_check['verified']:
        raise RuntimeError('Preallocated red spacing changed during Escher rendering')
    map_nodes={key:old_nodes[key] for key in sources}
    node_number=max(map(int,sources))+1
    segment_number=max(map(int,body['reactions']))+1
    source_edges={};priorities={};routing=[]
    for reaction in body['reactions'].values():reaction['segments']={}
    for route in sorted(routes,key=lambda r:(int(r['reaction']),r['mid'])):
        path=[route['marker']]
        for point in route['points'][1:-1]:
            key=str(node_number);node_number+=1
            map_nodes[key]=dict(node_type='multimarker',x=point[0],y=point[1]);routing.append(key);path.append(key)
        path.append(route['target'])
        for a,b in zip(path,path[1:]):
            sid=str(segment_number);segment_number+=1
            body['reactions'][route['reaction']]['segments'][sid]=dict(from_node_id=a,to_node_id=b,b1=None,b2=None)
            source_edges[sid]=route['pair'];priorities[sid]=route['priority']
    body['nodes']=map_nodes
    meta.update(routing_node_ids=routing,segment_source_edges=source_edges,segment_priorities=priorities)
    meta['overlap_routing']=dict(routing_preference=POLICY,important_simplification=simplify_report,gray_yield=gray_report,
                                red_layout_signature=(red_plan or {}).get('signature'),red_parallel_spacing=red_check,
                                trunk_grouping='positive_length_overlap',
                                node_obstacle_policy='red_orange_structural_nodes',obstacle_node_ids=obstacle_keys,
                                node_detours=node_report,parallel_detours=parallel_report,enabled=enabled,base_map_id=base_id,lane_spacing=step if enabled else 0,
                                shifted_legs=sum(bool(l['offset']) for l in legs),
                                max_offset=max((abs(l['offset']) for l in legs),default=0),
                                priority=['red','orange','gray'],shared_same_reaction_trunks=True,
                                max_offset_by_priority={name:max((abs(l['offset']) for l in legs if l['route']['priority']==priority),default=0) for name,priority in [('red',2),('orange',1),('gray',0)]},
                                primary_backbone_nodes_fixed=True,
                                primary_backbone_max_offset=max(max((abs(l['offset']) for l in legs if l['route']['backbone']),default=0),
                                    (node_report or {}).get('backbone_max_offset_grid',0)*meta['coordinate_transform']['scale']))
    if enabled:
        signature=[base_id,[(r['rid'],r['mid'],r['priority'],r['points']) for r in routes]]
        header['map_id']=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:24]
        canvas=body['canvas'];padding=meta['coordinate_transform']['padding']
        x1=min(canvas['x'],min(n['x'] for n in map_nodes.values())-padding)
        y1=min(canvas['y'],min(n['y'] for n in map_nodes.values())-padding)
        x2=max(canvas['x']+canvas['width'],max(n['x'] for n in map_nodes.values())+padding)
        y2=max(canvas['y']+canvas['height'],max(n['y'] for n in map_nodes.values())+padding)
        canvas.update(x=x1,y=y1,width=x2-x1,height=y2-y1)
    else:header['map_id']=base_id
    return data
