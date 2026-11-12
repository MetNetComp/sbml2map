"""Hard red crossing constraints and length-first orange routing.

Biological endpoints stay fixed. Same-reaction branches may share a connected
trunk; unrelated touching/overlap is forbidden whenever either branch is red.
"""
import bisect
import heapq
import math
import copy
import random
import hashlib
import json
from collections import defaultdict
import numpy as np
from tca_curves import arc_controls, curved_pairs, sample_curve

EPS=1e-7
POLICY='exchange_boundary_red_crossings_forbidden_orange_length_first'


def clean(points):
    result=[]
    for p in points:
        p=tuple(map(float,p))
        if result and math.dist(result[-1],p)<EPS:continue
        while len(result)>1:
            a,b=result[-2:]
            if not ((abs(a[0]-b[0])<EPS and abs(b[0]-p[0])<EPS) or
                    (abs(a[1]-b[1])<EPS and abs(b[1]-p[1])<EPS)):break
            result.pop()
        result.append(p)
    return result


def _curve_line(curve,a,b):
    axis=1 if abs(a[1]-b[1])<EPS else 0
    values=[curve['points'][0],*curve['curve_controls'],curve['points'][-1]]
    p0,p1,p2,p3=[p[axis] for p in values]
    polynomial=[-p0+3*p1-3*p2+p3,3*p0-6*p1+3*p2,-3*p0+3*p1,p0-a[axis]]
    nonzero=next((i for i,v in enumerate(polynomial) if abs(v)>1e-10),len(polynomial))
    if nonzero==len(polynomial):return []
    out=[]
    for root in np.roots(polynomial[nonzero:]):
        if abs(root.imag)>1e-7 or not -EPS<=root.real<=1+EPS:continue
        t=max(0.,min(1.,float(root.real)))
        point=tuple(sum(w*p[k] for w,p in zip(((1-t)**3,3*(1-t)**2*t,3*(1-t)*t*t,t**3),values)) for k in (0,1))
        side=1-axis
        if min(a[side],b[side])-EPS<=point[side]<=max(a[side],b[side])+EPS:
            if not any(math.dist(point,q)<1e-6 for q in out):out.append(point)
    return out


def intersections(first,second):
    """Return point contacts and positive overlaps of the actual saved geometry."""
    if first.get('curve_controls') and second.get('curve_controls'):
        # TCA arcs share one validated circle; source order selects disjoint arcs.
        return [(p,p) for p in (first['points'][0],first['points'][-1])
                for q in (second['points'][0],second['points'][-1]) if math.dist(p,q)<EPS]
    if first.get('curve_controls') or second.get('curve_controls'):
        curve,line=(first,second) if first.get('curve_controls') else (second,first)
        return [(p,p) for a,b in zip(line['points'],line['points'][1:]) for p in _curve_line(curve,a,b)]
    found=[]
    for a,b in zip(first['points'],first['points'][1:]):
        ha=abs(a[1]-b[1])<EPS
        for c,d in zip(second['points'],second['points'][1:]):
            hb=abs(c[1]-d[1])<EPS
            if ha==hb:
                axis=0 if ha else 1
                if abs(a[1-axis]-c[1-axis])>EPS:continue
                lo=max(min(a[axis],b[axis]),min(c[axis],d[axis]));hi=min(max(a[axis],b[axis]),max(c[axis],d[axis]))
                if lo>hi+EPS:continue
                p=list(a);q=list(a);p[axis]=lo;q[axis]=max(lo,hi);found.append((tuple(p),tuple(q)))
            else:
                h,v=(a,c) if ha else (c,a)
                x,y=v[0],h[1]
                if (min(a[0],b[0])-EPS<=x<=max(a[0],b[0])+EPS and min(a[1],b[1])-EPS<=y<=max(a[1],b[1])+EPS and
                    min(c[0],d[0])-EPS<=x<=max(c[0],d[0])+EPS and min(c[1],d[1])-EPS<=y<=max(c[1],d[1])+EPS):found.append(((x,y),(x,y)))
    return found


def _on(point,a,b):
    return abs(math.dist(a,point)+math.dist(point,b)-math.dist(a,b))<1e-6


def measure(paths):
    important=[p for p in paths if p['priority']>0]
    counts={'red_red':0,'red_orange':0,'orange_orange':0};details=[]
    for i,first in enumerate(important):
        for second in important[i+1:]:
            hits=intersections(first,second)
            if not hits:continue
            shared=set((first['source'],first['target']))&set((second['source'],second['target']))
            common=[first['points'][0] if k==first['source'] else first['points'][-1] for k in shared]
            # Only a connected shared prefix of the same reaction is a legitimate trunk.
            trunk=[]
            if first['rid']==second['rid']:
                remaining=[(a,b) for a,b in hits if math.dist(a,b)>EPS]
                reached=list(common)
                while remaining:
                    ready=[(a,b) for a,b in remaining if any(_on(p,a,b) for p in reached)]
                    if not ready:break
                    for segment in ready:remaining.remove(segment);trunk.append(segment);reached.extend(segment)
            bad=[]
            for a,b in hits:
                if any(_on(a,c,d) and _on(b,c,d) for c,d in trunk):continue
                if math.dist(a,b)<1e-6 and any(math.dist(a,p)<1e-6 for p in common):continue
                if not any(math.dist(a,c)<1e-6 and math.dist(b,d)<1e-6 for c,d in bad):bad.append((a,b))
            # An overlap and its two endpoints describe one conflict, not three.
            overlaps=[(a,b) for a,b in bad if math.dist(a,b)>1e-6]
            bad=overlaps+[(a,b) for a,b in bad if math.dist(a,b)<=1e-6 and not any(_on(a,c,d) for c,d in overlaps)]
            if bad:
                category='red_red' if first['priority']==second['priority']==2 else 'red_orange' if 2 in (first['priority'],second['priority']) else 'orange_orange'
                counts[category]+=len(bad)
                details.append(dict(category=category,reactions=[first['rid'],second['rid']],metabolites=[first['mid'],second['mid']],contacts=bad))
    return dict(**counts,hard_verified=counts['red_red']==counts['red_orange']==0,violations=details)


def orange_length(paths):
    """Physical orange length; connected same-reaction trunks count once."""
    groups=defaultdict(list)
    for path in paths:
        if path['priority']!=1:continue
        for a,b in zip(path['points'],path['points'][1:]):
            axis=0 if abs(a[1]-b[1])<EPS else 1
            groups[(path['rid'],axis,round(a[1-axis],8))].append(sorted((a[axis],b[axis])))
    total=0.
    for segments in groups.values():
        low=high=None
        for a,b in sorted(segments):
            if high is None:low,high=a,b
            elif a<=high+EPS:high=max(high,b)
            else:total+=high-low;low,high=a,b
        if high is not None:total+=high-low
    return total


def validate_boundary(paths,boundary):
    """Check all colors; cubic extrema, rather than only endpoints, are tested."""
    xmin,xmax,ymin,ymax=boundary
    bad=[];extents=[]
    for path in paths:
        points=list(path['points'])
        if path.get('curve_controls'):
            controls=[points[0],*path['curve_controls'],points[-1]]
            for axis in (0,1):
                p0,p1,p2,p3=[p[axis] for p in controls]
                polynomial=[3*(-p0+3*p1-3*p2+p3),2*(3*p0-6*p1+3*p2),-3*p0+3*p1]
                first=next((i for i,v in enumerate(polynomial) if abs(v)>1e-12),len(polynomial))
                for root in np.roots(polynomial[first:]) if first<len(polynomial) else ():
                    if abs(root.imag)<EPS and 0<root.real<1:
                        t=float(root.real);points.append(tuple(sum(w*p[a] for w,p in zip(((1-t)**3,3*(1-t)**2*t,3*(1-t)*t*t,t**3),controls)) for a in (0,1)))
        extents.extend(points)
        if any(not (xmin-EPS<=p[0]<=xmax+EPS and ymin-EPS<=p[1]<=ymax+EPS) for p in points):
            bad.append(dict(reaction_id=path['rid'],metabolite_id=path['mid'],priority=path['priority']))
    result=dict(boundary=list(boundary),outside_path_count=len(bad),outside_paths=bad,all_colors_checked=True,cubic_extrema_checked=True,verified=not bad)
    if extents:result['route_extents']=[min(p[0] for p in extents),max(p[0] for p in extents),min(p[1] for p in extents),max(p[1] for p in extents)]
    if bad:raise ValueError('Paths extend beyond the exchange boundary: '+str(bad))
    return result


def initial_paths(nodes,edges,xy,backbones,major):
    orange=set((backbones.get('exchange_connectors') or {}).get('edge_indices',()))
    curves=curved_pairs(backbones);region=backbones.get('tca_exclusion',{}).get('region')
    backbone={frozenset(('rxn:'+step['reaction_id'],'met:'+mid)) for r in backbones.get('routes',[]) for step in r.get('steps',[]) for mid in (step['source'],step['target'])}
    out=[]
    for i,e in enumerate(edges):
        marker='rxn:'+e['reaction_id'];target=e['target'] if e['source']==marker else e['source'];a,b=tuple(xy[marker]),tuple(xy[target])
        p=dict(rid=e['reaction_id'],mid=e['metabolite_id'],source=marker,target=target,
               priority=2 if e['reaction_id'] in major else 1 if i in orange else 0,
               backbone=frozenset((marker,target)) in backbone,points=clean([a,(b[0],a[1]),b]))
        if frozenset((marker,target)) in curves:
            p['curve_controls']=arc_controls(a,b,region['center']);p['points']=[a,b]
        out.append(p)
    return out


class RoutingGrid:
    def __init__(self,paths,xy,margin=3.,boundary=None):
        self.paths=paths;self.xy=xy
        self.keys={k for p in paths if p['priority'] or p.get('straight_locked') for k in (p['source'],p['target'])}
        positions=[xy[k] for k in self.keys]
        # Coordinate compression offers local lanes without a huge uniform raster.
        self.x=np.array(sorted({round(p[0]+d,8) for p in positions for d in (-.15,0,.15)}|{min(p[0] for p in positions)-margin,max(p[0] for p in positions)+margin}))
        self.y=np.array(sorted({round(p[1]+d,8) for p in positions for d in (-.15,0,.15)}|{min(p[1] for p in positions)-margin,max(p[1] for p in positions)+margin}))
        fixed_points=[point for p in paths if p.get('red_fixed') for point in p['points']]
        self.x=np.array(sorted(set(self.x)|{float(p[0]) for p in fixed_points}))
        self.y=np.array(sorted(set(self.y)|{float(p[1]) for p in fixed_points}))
        self.boundary=boundary
        if boundary is not None:
            xmin,xmax,ymin,ymax=boundary
            self.x=np.array(sorted({xmin,xmax,xmin+.15,xmax-.15}|{float(v) for v in self.x if xmin-EPS<=v<=xmax+EPS}))
            self.y=np.array(sorted({ymin,ymax,ymin+.15,ymax-.15}|{float(v) for v in self.y if ymin-EPS<=v<=ymax+EPS}))
        self.nx=len(self.x);self.ny=len(self.y);self.size=self.nx*self.ny
        self.node=np.zeros(self.size,dtype=bool)
        self.lookup={k:self.index(p) for k,p in xy.items() if k in self.keys}
        for v in self.lookup.values():self.node[v]=True
        self.traces={};self.curves=[p for p in paths if p.get('curve_controls')]
        self.arc_h=np.zeros((self.ny,self.nx-1),dtype=np.int16);self.arc_v=np.zeros((self.ny-1,self.nx),dtype=np.int16)
        self.events=defaultdict(list)
        for curve in self.curves:self.add_curve(curve)

    def index(self,p):
        ix=int(np.argmin(abs(self.x-p[0])));iy=int(np.argmin(abs(self.y-p[1])))
        return iy*self.nx+ix

    def point(self,index):return (float(self.x[index%self.nx]),float(self.y[index//self.nx]))

    def add_curve(self,curve):
        controls=[curve['points'][0],*curve['curve_controls'],curve['points'][-1]]
        for axis,levels in ((1,self.y),(0,self.x)):
            low=min(p[axis] for p in controls)-EPS;high=max(p[axis] for p in controls)+EPS
            for row in range(bisect.bisect_left(levels,low),bisect.bisect_right(levels,high)):
                if axis==1:a,b=(self.x[0],levels[row]),(self.x[-1],levels[row])
                else:a,b=(levels[row],self.y[0]),(levels[row],self.y[-1])
                for point in _curve_line(curve,a,b):
                    values=self.x if axis==1 else self.y;value=point[1-axis];slot=bisect.bisect_right(values,value)-1
                    slots={max(0,min(len(values)-2,slot))}
                    near=int(np.argmin(abs(values-value)))
                    if abs(values[near]-value)<1e-6:slots.update(k for k in (near-1,near) if 0<=k<len(values)-1)
                    for col in slots:
                        key=('h',row,col) if axis==1 else ('v',col,row)
                        array=self.arc_h if axis==1 else self.arc_v
                        array[key[1],key[2]]+=1;self.events[key].append((curve,point))

    def reserve(self,path,trace):
        self.traces[id(path)]=trace

    def masks(self,path):
        blocked=self.node.copy();start=self.lookup[path['source']];end=self.lookup[path['target']];blocked[[start,end]]=False
        if self.boundary is not None:
            # Boundary lanes cannot serve as through routes. EX endpoints only
            # connect to an inward neighbor; the outer side never exists.
            blocked[:self.nx]=True;blocked[-self.nx:]=True
            blocked[::self.nx]=True;blocked[self.nx-1::self.nx]=True
            blocked[[start,end]]=False

        hard_h=self.arc_h.copy();hard_v=self.arc_v.copy()
        for key,events in self.events.items():
            allowed=0
            for curve,point in events:
                common=set((path['source'],path['target']))&set((curve['source'],curve['target']))
                if any(math.dist(point,self.xy[k])<1e-6 for k in common):allowed+=1
            if allowed:(hard_h if key[0]=='h' else hard_v)[key[1],key[2]]-=allowed
        hard_h=hard_h>0;hard_v=hard_v>0
        soft=np.zeros(self.size,dtype=np.int16);soft_h=np.zeros_like(hard_h,dtype=np.int16);soft_v=np.zeros_like(hard_v,dtype=np.int16)
        own_h=np.zeros_like(hard_h,dtype=bool);own_v=np.zeros_like(hard_v,dtype=bool)
        for other in self.paths:
            trace=self.traces.get(id(other))
            if trace is None:continue
            if other['rid']==path['rid']:
                oh,ov=own_h,own_v
            elif other['priority']==2 or (path['priority']==2 and other.get('straight_locked')):
                blocked[trace]=True;oh,ov=hard_h,hard_v
            else:
                soft[trace]+=1;oh,ov=soft_h,soft_v
            for a,b in zip(trace,trace[1:]):
                if a//self.nx==b//self.nx:oh[a//self.nx,min(a%self.nx,b%self.nx)]=True if oh.dtype==bool else oh[a//self.nx,min(a%self.nx,b%self.nx)]+1
                else:ov[min(a//self.nx,b//self.nx),a%self.nx]=True if ov.dtype==bool else ov[min(a//self.nx,b//self.nx),a%self.nx]+1
        blocked[[start,end]]=False;soft[[start,end]]=0
        # Red must also avoid already routed orange when used for a later repair.
        if path['priority']==2:
            blocked|=soft>0;blocked[[start,end]]=False;hard_h|=soft_h>0;hard_v|=soft_v>0
        return blocked,hard_h,hard_v,soft,soft_h,soft_v,own_h,own_v

    def search(self,path):
        start=self.lookup[path['source']];target=self.lookup[path['target']]
        blocked,hh,hv,soft,sh,sv,oh,ov=self.masks(path)
        x,y=self.x,self.y;nx=self.nx;ny=self.ny
        # Exact lexicographic Dijkstra: length first, then bends. Orange
        # occupancy never contributes a crossing penalty or blocks orange.
        initial=start*3;distance={initial:(0.,0)};previous={};queue=[(0.,0,initial)]
        while queue:
            length_cost,bends,state=heapq.heappop(queue);cost=(length_cost,bends)
            if cost!=distance.get(state):continue
            v,axis=divmod(state,3);ix,iy=v%nx,v//nx
            if v==target:
                trace=[v]
                while state!=initial:state=previous[state];trace.append(state//3)
                trace.reverse();return trace
            options=[]
            if ix:options.append((v-1,1,x[ix]-x[ix-1],hh[iy,ix-1],oh[iy,ix-1]))
            if ix+1<nx:options.append((v+1,1,x[ix+1]-x[ix],hh[iy,ix],oh[iy,ix]))
            if iy:options.append((v-nx,2,y[iy]-y[iy-1],hv[iy-1,ix],ov[iy-1,ix]))
            if iy+1<ny:options.append((v+nx,2,y[iy+1]-y[iy],hv[iy,ix],ov[iy,ix]))
            for other,direction,length,wall,own in options:
                if wall or blocked[other]:continue
                # A shared orange trunk has already paid its physical length.
                step=0. if own and path['priority']==1 else float(length)*(.05 if own else 1.)
                value=(round(length_cost+step,8),bends+int(bool(axis and axis!=direction)))
                following=other*3+direction
                if value>=distance.get(following,(math.inf,math.inf)):continue
                distance[following]=value;previous[following]=state
                heapq.heappush(queue,(*value,following))
        return None


def route(paths,xy,margin=3.,progress=print,attempts=24,seed=17,boundary=None):
    important=[p for p in paths if (p['priority'] or p.get('straight_locked')) and not p.get('curve_controls')]
    if boundary is not None:validate_boundary(paths,boundary)
    if not important:return dict(policy=POLICY,verified=True,metrics=measure(paths),nodes_moved=False,attempts=0,
        boundary=boundary,boundary_validation=validate_boundary(paths,boundary) if boundary is not None else None,orange_length_grid=orange_length(paths),orange_crossings_optimized=False)
    grid=RoutingGrid(paths,xy,margin,boundary)
    locked=[p for p in important if p.get('straight_locked') or p.get('red_fixed')]
    for p in locked:
        trace=[]
        for first,last in zip(p['points'],p['points'][1:]):
            a,b=grid.index(first),grid.index(last)
            step=(1 if b>a else -1) if a//grid.nx==b//grid.nx else (grid.nx if b>a else -grid.nx)
            trace.extend(range(a,b,step))
        trace.append(grid.index(p['points'][-1]));grid.reserve(p,trace)
    ordered=sorted([p for p in important if not p.get('straight_locked') and not p.get('red_fixed')],key=lambda p:(-p['priority'],not p['backbone'],math.dist(p['points'][0],p['points'][-1]),p['rid'],p['mid']))
    for index,path in enumerate(ordered):
        trace=grid.search(path)
        if trace is None:
            return dict(verified=False,reason='No crossing-free path in the current fixed embedding',reaction_id=path['rid'],metabolite_id=path['mid'],routed=index,total=len(ordered))
        path['points']=clean([grid.point(v) for v in trace]);grid.reserve(path,trace)
        if index%20==0:progress(f'Crossing constraints: routed {index+1}/{len(ordered)}',flush=True)
    metrics=measure(paths)
    baseline=copy.deepcopy(metrics)
    first_length=orange_length(paths)
    def score():
        m=measure(paths)
        return (m['red_red']+m['red_orange'],round(orange_length(paths),8),
                sum(len(p['points'])-2 for p in paths if p['priority']==1))
    best_score=score();best_points={id(p):list(p['points']) for p in important}
    orange=[p for p in ordered if p['priority']==1]
    red_traces={id(p):grid.traces[id(p)] for p in important if p['priority']==2 or p.get('straight_locked')}
    generator=random.Random(seed)
    for attempt in range(attempts):
        grid.traces=dict(red_traces)
        ordering=list(orange)
        if attempt==0:ordering.sort(key=lambda p:-math.dist(p['points'][0],p['points'][-1]))
        else:generator.shuffle(ordering)
        for path in ordering:
            trace=grid.search(path)
            if trace is None:break
            path['points']=clean([grid.point(v) for v in trace]);grid.reserve(path,trace)
        else:
            candidate=score()
            if candidate<best_score:
                best_score=candidate;best_points={id(p):list(p['points']) for p in important}
                progress(f'Orange length search {attempt+1}/{attempts}: {best_score[1]:.4f} grid units',flush=True)
    for path in important:path['points']=best_points[id(path)]
    metrics=measure(paths)
    return dict(policy=POLICY,verified=metrics['hard_verified'],metrics=metrics,initial_feasible_metrics=baseline,nodes_moved=False,grid_shape=[grid.nx,grid.ny],
                boundary=boundary,boundary_validation=validate_boundary(paths,boundary) if boundary is not None else None,
                attempts=attempts,objective_order=['red crossings forbidden','orange total physical length','bends'],
                orange_crossings_optimized=False,orange_length_grid=orange_length(paths),first_feasible_orange_length_grid=first_length,
                length_definition='Union of orange segments within each reaction; different reactions count separately',global_optimality_proven=False)


def plan_routes(nodes,edges,xy,backbones,major,attempts=24,progress=print):
    paths=initial_paths(nodes,edges,xy,backbones,major)
    red_phase=backbones.get('red_first')
    if red_phase:
        saved={(p['rid'],p['mid']):p for p in red_phase['routes']}
        for path in paths:
            if path['priority']==2:
                record=saved[path['rid'],path['mid']]
                path['points']=list(record['points']);path['red_fixed']=True
    corridors=backbones.get('exchange_corridors')
    if corridors:
        from exchange_corridors import lock_corridor_paths
        lock_corridor_paths(paths,corridors)
    before=measure(paths)
    boundary=backbones.get('routing_boundary')
    report=route(paths,xy,progress=progress,attempts=attempts,boundary=boundary)
    report['before']=before
    report['red_first']=bool(red_phase)
    if red_phase:report['red_routes_preserved']=all(p['points']==saved[p['rid'],p['mid']]['points'] for p in paths if p['priority']==2)
    if report['verified'] and corridors:
        from exchange_corridors import validate_corridors
        report['straight_exchange_validation']=validate_corridors(paths,corridors)
    if not report['verified']:
        raise ValueError('Hard red crossing constraints are infeasible in the current embedding: '+str(report))
    return dict(routes=paths,report=report,signature=hashlib.sha256(json.dumps(paths,sort_keys=True).encode()).hexdigest()[:24])


def apply_to_escher(data,plan):
    from tca_curves import transform_point
    header,body=data;meta=header['metadata'];transform=meta['coordinate_transform']
    sources=meta['source_node_ids'];lookup={v:k for k,v in sources.items()};nodes={k:body['nodes'][k] for k in sources}
    records={(p['rid'],p['mid']):p for p in plan['routes']};next_node=max(map(int,nodes))+1;next_segment=1
    edges={};priorities={};routing=[];curves=[]
    for reaction in body['reactions'].values():
        reaction['segments']={}
        for item in reaction['metabolites']:
            p=records[reaction['bigg_id'],item['bigg_id']]
            points=[transform_point(point,transform) for point in p['points']]
            chain=[lookup[p['source']]]
            for point in points[1:-1]:
                key=str(next_node);next_node+=1;nodes[key]=dict(node_type='multimarker',x=point[0],y=point[1]);routing.append(key);chain.append(key)
            chain.append(lookup[p['target']])
            controls=[dict(zip(('x','y'),transform_point(c,transform))) for c in p.get('curve_controls',())]
            for a,b in zip(chain,chain[1:]):
                sid=str(next_segment);next_segment+=1
                reaction['segments'][sid]=dict(from_node_id=a,to_node_id=b,b1=controls[0] if controls else None,b2=controls[1] if controls else None)
                edges[sid]=[p['source'],p['target']];priorities[sid]=p['priority']
                if controls:curves.append(sid)
    body['nodes']=nodes
    meta.update(routing_node_ids=routing,segment_source_edges=edges,segment_priorities=priorities,
                tca_curves=dict(style='circular_bezier',segment_ids=curves,curved_incidences=len(curves)),
                crossing_routing=plan['report'],crossing_routing_signature=plan['signature'])
    # The same saved geometry is used for previews and static rendering.
    from routing_simplification import POLICY as LEGACY_POLICY
    from red_route_layout import plan_from_metadata
    meta['overlap_routing']=dict(enabled=True,routing_preference=LEGACY_POLICY,node_obstacle_policy='red_orange_structural_nodes',
        trunk_grouping='positive_length_overlap',red_layout_signature=(plan_from_metadata(meta) or {}).get('signature'),hard_crossing_policy=plan['report'].get('policy',POLICY),
        preallocated_parallel_channels_applied=False)
    header['map_id']=plan['signature']
    canvas=body['canvas'];padding=transform['padding']
    xs=[n['x'] for n in nodes.values()];ys=[n['y'] for n in nodes.values()]
    x1=min(canvas['x'],min(xs)-padding);x2=max(canvas['x']+canvas['width'],max(xs)+padding)
    y1=min(canvas['y'],min(ys)-padding);y2=max(canvas['y']+canvas['height'],max(ys)+padding)
    canvas.update(x=x1,y=y1,width=x2-x1,height=y2-y1)
    if plan['report'].get('uniform_node_sizes'):
        fit_escher_canvas(body)
    return data


def validate_saved_geometry(data,metadata=None):
    """Reconstruct actual saved paths; never trust a stored crossing count."""
    from tca_curves import backbone_from_metadata
    meta=metadata or data[0]['metadata'];body=data[1];sources=meta['source_node_ids'];t=meta['coordinate_transform']
    backbones=backbone_from_metadata(meta);major={rid for r in backbones.get('routes',[]) for rid in r['reaction_ids']}
    canonical={(p['rid'],p['mid']):p for p in backbones.get('crossing_routing',{}).get('routes',())}
    inverse=lambda p:((p['x']-t['padding'])/t['scale']+t['x_min'],t['y_max']-(p['y']-t['padding'])/t['scale'])
    paths=[]
    for reaction in body['reactions'].values():
        rid=reaction['bigg_id'];marker=meta['reactions'][rid]['marker_node_id'];adj=defaultdict(list)
        for seg in reaction['segments'].values():
            a,b=seg['from_node_id'],seg['to_node_id'];adj[a].append((b,seg));adj[b].append((a,seg))
        for target in (k for k in adj if k in sources and k!=marker):
            mid=body['nodes'][target]['bigg_id'];chain=[target];segments=[];previous=None
            while chain[-1]!=marker:
                options=[(v,s) for v,s in adj[chain[-1]] if v!=previous]
                if len(options)!=1:raise ValueError('Invalid saved incidence for crossing verification')
                following,segment=options[0];segments.append(segment);previous=chain[-1];chain.append(following)
            chain.reverse();segments.reverse()
            record=canonical[rid,mid]
            p=dict(rid=rid,mid=mid,source=sources[marker],target=sources[target],priority=record['priority'] if backbones.get('adaptive_layout') else 2 if rid in major else record['priority'],points=[inverse(body['nodes'][n]) for n in chain])
            p['backbone']=record.get('backbone',False)
            if record.get('auxiliary_lane'):p['auxiliary_lane']=True
            if record.get('local_cofactor'):p['local_cofactor']=True
            if len(segments)==1 and segments[0]['b1'] is not None:
                controls=[segments[0]['b1'],segments[0]['b2']]
                if segments[0]['from_node_id']!=marker:controls.reverse()
                p['curve_controls']=[inverse(c) for c in controls]
            paths.append(p)
    if backbones.get('adaptive_layout'):
        from adaptive_routing import validate_paths
        xy={source:inverse(body['nodes'][key]) for key,source in sources.items()}
        result=validate_paths(paths,xy,backbones['routing_boundary'],backbones['tca_exclusion'].get('region'),
                              backbones['crossing_routing']['report']['conditional_crossings'])
        optimization=backbones.get('red_optimization')
        if optimization:
            from red_objective import preservation_report,relaxation_level
            reference=dict(positions=optimization['reference_positions'],routes=optimization['reference_routes'],
                           length=optimization['optimized_length_grid'],signature=optimization['reference_signature'])
            check=preservation_report(reference,xy,paths)
            if abs(check['final_length_grid']-optimization['preservation']['final_length_grid'])>1e-6:
                raise ValueError('Saved red length differs from its preservation report')
            if relaxation_level(reference,dict(positions=xy,routes=paths))!=optimization['relaxation_level']:
                raise ValueError('Saved red geometry violates its preservation level')
            result['red_preservation']=check
        return result
    result=measure(paths)
    boundary=backbones.get('routing_boundary')
    if boundary is not None:result['boundary_validation']=validate_boundary(paths,boundary)
    if backbones.get('exchange_corridors'):
        from exchange_corridors import validate_corridors
        result['straight_exchange_validation']=validate_corridors(paths,backbones['exchange_corridors'])
    if not result['hard_verified']:raise ValueError('Saved geometry violates hard red crossing constraints: '+str(result['violations']))
    return {k:v for k,v in result.items() if k!='violations'}


def fit_escher_canvas(body):
    canvas=body['canvas']
    # Fit the saved geometry and full label widths rather than keeping
    # the initial canvas's fixed 350-unit right-hand allowance.
    bounds=[(n['x']-10,n['y']-10,n['x']+10,n['y']+10) for n in body['nodes'].values()]
    bounds += [(n['label_x'],n['label_y']-20,n['label_x']+len(n['bigg_id'])*12,n['label_y']+5)
               for n in body['nodes'].values() if n['node_type']=='metabolite']
    bounds += [(r['label_x'],r['label_y']-20,r['label_x']+len(r['bigg_id'])*12,r['label_y']+5) for r in body['reactions'].values()]
    bounds += [(t['x'],t['y']-20,t['x']+len(t['text'])*12,t['y']+5) for t in body['text_labels'].values()]
    for reaction in body['reactions'].values():
        for segment in reaction['segments'].values():
            for key in ('b1','b2'):
                if segment[key] is not None:
                    p=segment[key];bounds.append((p['x'],p['y'],p['x'],p['y']))
    margin=12.
    x1=min(b[0] for b in bounds)-margin;y1=min(b[1] for b in bounds)-margin
    x2=max(b[2] for b in bounds)+margin;y2=max(b[3] for b in bounds)+margin
    canvas.update(x=x1,y=y1,width=x2-x1,height=y2-y1)
