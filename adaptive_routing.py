"""Integer-grid routing with disjoint orange lanes and audited red crossings."""
import copy
import hashlib
import heapq
import json
import math
from collections import defaultdict
import numpy as np
from crossing_routing import EPS,clean,initial_paths,intersections,measure,orange_length,validate_boundary,_on
from exchange_corridors import place_corridors,lock_corridor_paths,validate_corridors
from tca_geometry import allowed_points,validate_exclusion

POLICY='adaptive_grid_disjoint_orange_conditional_red_crossings'
PORT_RADIUS=.5
PARALLEL_SPACING=1.
PARALLEL_LENGTH=2.


def make_paths(nodes,edges,xy,backbones,major):
    paths=initial_paths(nodes,edges,xy,backbones,major)
    for p in paths:
        if nodes[p['target']]['kind']=='currency':
            p['priority']=0
            p['local_cofactor']=True
            region=backbones.get('tca_exclusion',{}).get('region')
            if p['source'] in set((region or {}).get('node_ids',[])):
                a,b=xy[p['source']],xy[p['target']];options=[]
                for axis in (0,1):
                    port=list(a);port[axis]=b[axis]
                    if not integer(a[1-axis]) or not integer(port[axis]):continue
                    if not incident_exempt(a,port,p,xy,region):continue
                    if segment_clearance(port,b,region)<1-EPS:continue
                    options.append(clean([a,port,b]))
                if options:p['points']=min(options,key=lambda points:sum(math.dist(a,b) for a,b in zip(points,points[1:])))
    return paths


def integer(x):return abs(x-round(x))<EPS

def in_common_port(a,b,first,second,xy):
    common={first['source'],first['target']}&{second['source'],second['target']}
    return any(all(math.dist(p,xy[k])<=PORT_RADIUS+EPS for p in (a,b)) for k in common)


def parallel_conflicts(paths):
    """Check maximal straight runs of every color, including coincident lines."""
    legs=[];conflicts=[]
    for index,p in enumerate(paths):
        if p.get('curve_controls'):continue
        for a,b in zip(clean(p['points']),clean(p['points'])[1:]):
            axis=0 if abs(a[1]-b[1])<EPS else 1
            if abs(a[1-axis]-b[1-axis])>EPS:continue
            lo,hi=sorted((a[axis],b[axis]))
            if hi-lo>=PARALLEL_LENGTH-EPS:legs.append((index,axis,a[1-axis],lo,hi))
    for i,(pid,axis,level,lo,hi) in enumerate(legs):
        for qid,other_axis,other_level,left,right in legs[i+1:]:
            if pid==qid or axis!=other_axis or abs(level-other_level)>=PARALLEL_SPACING-EPS:continue
            overlap=min(hi,right)-max(lo,left)
            if overlap>=PARALLEL_LENGTH-EPS:
                conflicts.append(dict(first=[paths[pid]['rid'],paths[pid]['mid']],second=[paths[qid]['rid'],paths[qid]['mid']],spacing_grid=abs(level-other_level),length_grid=overlap))
    return conflicts


def validate_node_grid(nodes,xy,region):
    ring=set((region or {}).get('node_ids',[]))
    bad=[k for k,p in xy.items() if k not in ring and nodes[k]['kind']!='currency' and not all(integer(v) for v in p)]
    if bad:raise ValueError('Non-TCA ordinary nodes must be grid intersections: '+str(bad))
    return dict(verified=True,exceptions=['TCA perimeter','reaction-local cofactors'],checked_nodes=sum(k not in ring and nodes[k]['kind']!='currency' for k in xy))


def perimeter_distance(point, region):
    if not region:return math.inf
    if region['shape']=='circle':return abs(math.dist(point,region['center'])-region['radius'])
    x,y=point;l,t,s=region['left'],region['top'],region['side_length']
    if l<=x<=l+s and t-s<=y<=t:return min(x-l,l+s-x,y-(t-s),t-y)
    return math.hypot(max(l-x,0.,x-l-s),max(t-s-y,0.,y-t))


def segment_clearance(a,b,region):
    if not region:return math.inf
    if region['shape']=='circle':
        c=region['center'];v=np.asarray(b)-a;vv=float(v@v)
        u=max(0.,min(1.,float((np.asarray(c)-a)@v)/vv)) if vv else 0.
        low=math.dist(np.asarray(a)+u*v,c);high=max(math.dist(a,c),math.dist(b,c));r=region['radius']
        return 0. if low<=r<=high else min(abs(low-r),abs(high-r))
    l,t,s=region['left'],region['top'],region['side_length'];sides=[([l,t],[l+s,t]),([l+s,t],[l+s,t-s]),([l+s,t-s],[l,t-s]),([l,t-s],[l,t])]
    distance=math.inf
    for c,d in sides:
        first=dict(points=[a,b]);second=dict(points=[c,d])
        if intersections(first,second):return 0.
        for p,u,v in [(a,c,d),(b,c,d),(c,a,b),(d,a,b)]:
            vec=np.asarray(v)-u;den=float(vec@vec);k=max(0.,min(1.,float((np.asarray(p)-u)@vec)/den)) if den else 0.
            distance=min(distance,math.dist(p,np.asarray(u)+k*vec))
    return distance


def incident_exempt(a,b,path,xy,region):
    if not region:return False
    ring=set(region['node_ids'])
    if region['shape']=='square' and path.get('backbone') and path['source'] in ring and path['target'] in ring:return True
    for key in (path['source'],path['target']):
        if key not in ring:continue
        p=xy[key]
        if region['shape']=='circle':
            axis=0 if abs(a[1]-b[1])<EPS else 1
            if abs(a[1-axis]-b[1-axis])<EPS:
                center=region['center'];radius=region['radius'];perpendicular=a[1-axis]-center[1-axis]
                lo,hi=sorted((a[axis],b[axis]));closest=max(lo,min(hi,center[axis]))
                minimum=math.hypot(closest-center[axis],perpendicular)
                if minimum>=radius-EPS and abs(perpendicular)<=radius+1:
                    extent=math.sqrt(max(0.,(radius+1)**2-perpendicular**2))
                    left=max(lo,center[axis]-extent);right=min(hi,center[axis]+extent)
                    near=[]
                    for value in (left,right):
                        q=list(a);q[axis]=value;near.append(q)
                    tangent=list(a);tangent[axis]=closest
                    touches_elsewhere=minimum<=radius+EPS and math.dist(tangent,p)>EPS
                    if left<=right+EPS and not touches_elsewhere and max(math.dist(p,q) for q in near)<=1.5+EPS:return True
        if region['shape']=='square':
            l,t,side=region['left'],region['top'],region['side_length']
            # Several incidences need independent ports at a square-side
            # anchor. Permit short outward fan-out before reaching the full
            # one-grid clearance; never enter another side of the cycle.
            outward=((abs(p[0]-l)<EPS and max(a[0],b[0])<=l+EPS) or
                     (abs(p[0]-l-side)<EPS and min(a[0],b[0])>=l+side-EPS) or
                     (abs(p[1]-t)<EPS and min(a[1],b[1])>=t-EPS) or
                     (abs(p[1]-t+side)<EPS and max(a[1],b[1])<=t-side+EPS))
            if outward:
                # Validate the part inside the side's clearance band, even
                # after collinear graph edges have been merged for export.
                bands=[]
                if abs(p[0]-l)<EPS:bands.append((l-1,l,t-side-1,t+1))
                if abs(p[0]-l-side)<EPS:bands.append((l+side,l+side+1,t-side-1,t+1))
                if abs(p[1]-t)<EPS:bands.append((l-1,l+side+1,t,t+1))
                if abs(p[1]-t+side)<EPS:bands.append((l-1,l+side+1,t-side-1,t-side))
                for xmin,xmax,ymin,ymax in bands:
                    if abs(a[1]-b[1])<EPS and ymin-EPS<=a[1]<=ymax+EPS:
                        lo=max(min(a[0],b[0]),xmin);hi=min(max(a[0],b[0]),xmax)
                        near=[(lo,a[1]),(hi,a[1])] if lo<=hi+EPS else []
                    elif abs(a[0]-b[0])<EPS and xmin-EPS<=a[0]<=xmax+EPS:
                        lo=max(min(a[1],b[1]),ymin);hi=min(max(a[1],b[1]),ymax)
                        near=[(a[0],lo),(a[0],hi)] if lo<=hi+EPS else []
                    else:near=[]
                    if near and max(math.dist(p,q) for q in near)<=1.5+EPS:return True
        for axis in (0,1):
            if not(integer(p[1-axis]) and abs(a[1-axis]-p[1-axis])<EPS and abs(b[1-axis]-p[1-axis])<EPS):continue
            da,db=a[axis]-p[axis],b[axis]-p[axis]
            if da*db < -EPS:continue
            sign=1 if da+db>=0 else -1
            if region['shape']=='circle':
                projection=(p[axis]-region['center'][axis])*sign
                if projection<0:
                    discriminant=projection*projection-2*region['radius']+1
                    if discriminant<0:continue
                    # An inward lead may enter the clear interior, but cannot
                    # re-enter the opposite side's one-grid clearance band.
                    limit=-projection+math.sqrt(discriminant)
                    if max(abs(da),abs(db))>limit+EPS:continue
            else:
                l,t,side=region['left'],region['top'],region['side_length']
                lower=l if axis==0 else t-side;upper=l+side if axis==0 else t
                outward=(abs(p[axis]-lower)<EPS and sign<0) or (abs(p[axis]-upper)<EPS and sign>0)
                if outward:return True
                # Tangential runs along a square side are not attachment leads.
                if axis==0 and (abs(p[1]-t)<EPS or abs(p[1]-t+side)<EPS):continue
                if axis==1 and (abs(p[0]-l)<EPS or abs(p[0]-l-side)<EPS):continue
                inward=(abs(p[axis]-lower)<EPS and sign>0) or (abs(p[axis]-upper)<EPS and sign<0)
                if inward and max(abs(da),abs(db))>side-1+EPS:continue
            return True
    return False


def attach_grid_cofactors(nodes,xy,region,boundary,reserved_paths=()):
    groups=defaultdict(list)
    for k,n in nodes.items():
        if n['kind']=='currency':groups[n['parent_id']].append(k)
    for parent,keys in groups.items():
        p=xy[parent];candidates=[];candidate_routes=[];ports=[p]
        if region and parent in region['node_ids']:
            ports=[]
            for axis in (0,1):
                if not integer(p[1-axis]):continue
                for sign in (-1,1):
                    for distance in range(1,10):
                        port=list(p);port[axis]=float(math.ceil(p[axis])+distance-1 if sign>0 else math.floor(p[axis])-distance+1)
                        if perimeter_distance(port,region)>=1-EPS and math.dist(port,p)>.2:
                            dummy=dict(source=parent,target=parent)
                            if incident_exempt(p,port,dummy,xy,region):ports.append(port)
                            break
            ports.sort(key=lambda q:math.dist(q,p))
            extended=[]
            for port in ports:
                stem=0 if abs(port[1]-p[1])<EPS else 1
                sign=1 if port[stem]>p[stem] else -1
                for extra in range(5):
                    q=list(port);q[stem]+=sign*extra;extended.append(q)
            ports=extended
        for port in ports:
            lengths=(.45,.75) if region and parent in region['node_ids'] else (.45,.75,1.25,1.75,2.25,2.75)
            for length in lengths:
                for axis in (0,1):
                    if region and parent in region['node_ids']:
                        stem=0 if abs(port[1]-p[1])<EPS else 1
                        if axis==stem:continue
                    if not integer(port[1-axis]):continue
                    for sign in (-1,1):
                        q=list(port);q[axis]+=sign*length
                        if not (region and parent in region['node_ids']):q[1-axis]+=.2
                        if not(boundary[0]<q[0]<boundary[1] and boundary[2]<q[1]<boundary[3]):continue
                        if not allowed_points([q],region)[0]:continue
                        if any(math.dist(q,v)<.16-EPS for k,v in xy.items() if k!=parent and k not in keys):continue
                        elbow=[q[0],port[1]]
                        if min(segment_clearance(port,elbow,region),segment_clearance(elbow,q,region))<1-EPS:continue
                        candidates.append(q)
                        candidate_routes.append(clean([p,port,q]))
        if len(candidates)<len(keys):raise ValueError('No grid-line cofactor positions at '+parent)
        if region and parent in region['node_ids']:
            chosen=[];paths=[]
            for key in sorted(keys):
                for q,points in zip(candidates,candidate_routes):
                    if any(math.dist(q,v)<.16-EPS for v in chosen):continue
                    path=dict(rid=parent,mid=key,source=parent,target=key,priority=0,points=points)
                    if parallel_conflicts([*reserved_paths,*paths,path]):continue
                    chosen.append(q);paths.append(path);break
                else:raise ValueError('No local cofactor branch without long parallel overlap at '+parent)
            candidates=chosen
        for key,q in zip(sorted(keys),candidates):xy[key]=q


class Grid:
    def __init__(self,xy,boundary,region,auxiliary=False):
        self.xy=xy;self.region=region;self.boundary=boundary;self.auxiliary=auxiliary;self.reserved=[];self.parallel_intervals=defaultdict(list)
        xmin,xmax,ymin,ymax=boundary
        xs=set(float(x) for x in range(math.ceil(xmin),math.floor(xmax)+1));ys=set(float(y) for y in range(math.ceil(ymin),math.floor(ymax)+1))
        if auxiliary:
            xs.update(x+d for x in range(math.ceil(xmin),math.floor(xmax)) for d in (.25,.75))
            ys.update(y+d for y in range(math.ceil(ymin),math.floor(ymax)) for d in (.25,.75))
        xs.update(float(p[0]) for p in xy.values());ys.update(float(p[1]) for p in xy.values())
        self.x=sorted(xs);self.y=sorted(ys);self.nx=len(self.x);self.ny=len(self.y)
        self.lookup={k:self.index(p) for k,p in xy.items()};self.owners=defaultdict(set)
        for k,i in self.lookup.items():self.owners[i].add(k)
        self.red_edges=defaultdict(list);self.red_vertices=defaultdict(list);self.orange_edges=defaultdict(list)
        self.clearance={}
    def index(self,p):return min(range(len(self.y)),key=lambda i:abs(self.y[i]-p[1]))*self.nx+min(range(len(self.x)),key=lambda i:abs(self.x[i]-p[0]))
    def point(self,v):return (self.x[v%self.nx],self.y[v//self.nx])
    def neighbors(self,v):
        ix,iy=v%self.nx,v//self.nx;p=self.point(v)
        if integer(p[1]) or (self.auxiliary and any(abs(p[1]-math.floor(p[1])-f)<EPS for f in (.25,.75))):
            if ix:yield v-1,1
            if ix+1<self.nx:yield v+1,1
        if integer(p[0]) or (self.auxiliary and any(abs(p[0]-math.floor(p[0])-f)<EPS for f in (.25,.75))):
            if iy:yield v-self.nx,2
            if iy+1<self.ny:yield v+self.nx,2
    def reserve(self,p):
        if p.get('curve_controls'):return
        self.reserved.append(p)
        points=clean(p['points'])
        for a,b in zip(points,points[1:]):
            axis=0 if abs(a[1]-b[1])<EPS else 1
            lo,hi=sorted((a[axis],b[axis]));level=a[1-axis]
            if hi-lo>=PARALLEL_LENGTH-EPS:
                for bucket in range(math.floor(level)-1,math.floor(level)+2):
                    self.parallel_intervals[axis,bucket].append((level,lo,hi,p))
        edges=self.red_edges if p['priority']==2 else self.orange_edges if p['priority']==1 or p.get('straight_locked') else None
        if edges is None:return
        for a,b in zip(p['points'],p['points'][1:]):
            start,end=self.index(a),self.index(b);step=(1 if end>start else -1) if start//self.nx==end//self.nx else (self.nx if end>start else -self.nx)
            for v in range(start,end,step):
                w=v+step;key=tuple(sorted((v,w)));edges[key].append(p)
                if p['priority']==2:
                    self.red_vertices[v].append((p,1 if abs(step)==1 else 2));self.red_vertices[w].append((p,1 if abs(step)==1 else 2))
    def search(self,p,allow_red=False):
        start,end=self.lookup[p['source']],self.lookup[p['target']];endpoints={p['source'],p['target']}
        target=self.point(end);origin=(start,0,start);dist={origin:(0.,0,0)};prev={};queue=[(math.dist(self.point(start),target),0.,0,0,start,0,start)]
        while queue:
            _,length,bends,crossings,v,axis,run_start=heapq.heappop(queue);state=(v,axis,run_start)
            if dist.get(state)!=(length,bends,crossings):continue
            if v==end:
                trace=[v]
                while state!=origin:state=prev[state];trace.append(state[0])
                return clean([self.point(v) for v in reversed(trace)]),length
            a=self.point(v)
            for w,direction in self.neighbors(v):
                b=self.point(w);key=tuple(sorted((v,w)))
                signed_direction=direction if w>v else -direction
                if axis==-signed_direction:continue
                along=0 if direction==1 else 1
                next_start=run_start if axis==signed_direction else v
                run_origin=self.point(next_start)
                # A run may turn away before two grids. Keep its origin in
                # the search state so collinear graph edges count together.
                if axis==signed_direction and (a[along]-run_origin[along])*(b[along]-a[along]) < -EPS:continue
                low,high=sorted((run_origin[along],b[along]));level=a[1-along]
                nearby=[(lo,hi) for other,lo,hi,q in self.parallel_intervals[along,math.floor(level)]
                        if abs(level-other)<PARALLEL_SPACING-EPS]
                if any(abs(level-other)<PARALLEL_SPACING-EPS and min(high,hi)-max(low,lo)>=PARALLEL_LENGTH-EPS
                       for other,lo,hi,q in self.parallel_intervals[along,math.floor(level)]):continue
                if not any(min(high,hi)-max(low,lo)>EPS for lo,hi in nearby):next_start=w
                if axis and abs(axis)!=direction and not self.auxiliary and not (integer(a[0]) and integer(a[1])):continue
                if w not in (start,end) and {k for k in self.owners[w]-endpoints if not (p['priority'] and k.startswith('copy:'))}:continue
                if w not in (start,end) and (b[0] in (self.boundary[0],self.boundary[1]) or b[1] in (self.boundary[2],self.boundary[3])):continue
                # Positive-length orange overlap is forbidden even within a reaction.
                if p['priority'] in (1,2) and any(not in_common_port(a,b,p,q,self.xy) for q in self.orange_edges[key]):continue
                red=self.red_edges[key]
                if p['priority'] in (1,2) and any(not in_common_port(a,b,p,q,self.xy) for q in red):continue
                if key not in self.clearance:self.clearance[key]=segment_clearance(a,b,self.region)
                if self.clearance[key]<1-EPS and not incident_exempt(a,b,p,self.xy,self.region):continue
                hits=0;bad=False
                if p['priority'] in (1,2):
                    for vertex in (v,w):
                        entries=[(q,d) for q,d in self.red_vertices[vertex] if q['rid']!=p['rid']]
                        if not entries:continue
                        common=any(k in endpoints for q,_ in entries for k in (q['source'],q['target']) if self.lookup[k]==vertex)
                        if common or all(in_common_port(self.point(vertex),self.point(vertex),p,q,self.xy) for q,d in entries):continue
                        orientations={d for q,d in entries}
                        if p['priority']==2 or not allow_red or len(orientations)!=1 or direction in orientations:bad=True;break
                        if vertex==v and axis and abs(axis)!=direction:bad=True;break
                        hits+=int(vertex==w)
                if bad:continue
                step=abs(a[0]-b[0])+abs(a[1]-b[1]);value=(round(length+step,9),bends+int(bool(axis and abs(axis)!=direction)),crossings+hits)
                following=(w,signed_direction,next_start)
                if value>=dist.get(following,(math.inf,math.inf,math.inf)):continue
                dist[following]=value;prev[following]=state
                h=abs(b[0]-target[0])+abs(b[1]-target[1])
                heapq.heappush(queue,(value[0]+h,*value,w,signed_direction,next_start))
        return None


def route_paths(paths,xy,boundary,region,progress=lambda *a,**k:None,*,red_order=0,reservations=()):
    grid=Grid(xy,boundary,region);exceptions=[]
    for p in reservations:grid.reserve(p)
    ring=set((region or {}).get('node_ids',[]))
    for p in paths:
        if p.get('red_fixed') or p.get('local_cofactor'):grid.reserve(p)
    ordered=sorted([p for p in paths if not p.get('curve_controls') and not p.get('red_fixed') and not p.get('local_cofactor')],
                   key=lambda p:(-int(p['priority']==2),not p.get('straight_locked',False),-p['priority'],
                                 not (p.get('backbone',False) or p['source'] in ring or p['target'] in ring),
                                 (-1 if red_order and p['priority']==2 else 1)*math.dist(p['points'][0],p['points'][-1]),p['rid'],p['mid']))
    for i,p in enumerate(ordered):
        result=grid.search(p,allow_red=False)
        direct=math.dist(xy[p['source']],xy[p['target']])
        if p.get('straight_locked') and result and abs(result[1]-direct)>EPS:result=None
        if not result and not grid.auxiliary and not p.get('straight_locked'):
            aux=Grid(xy,boundary,region,auxiliary=True)
            for saved in grid.reserved:aux.reserve(saved)
            candidate=aux.search(p,allow_red=False)
            grid=aux
            if candidate:result=candidate
        if p['priority']==1:
            shorter=grid.search(p,allow_red=True)
            if p.get('straight_locked') and shorter and abs(shorter[1]-direct)>EPS:shorter=None
            strict=result[1] if result else math.inf
            if shorter and shorter[1]<strict-EPS and (strict>=2*shorter[1]-EPS or strict-shorter[1]>=4-EPS):
                result=shorter;exceptions.append(dict(reaction_id=p['rid'],metabolite_id=p['mid'],uncrossed_length_grid=None if math.isinf(strict) else strict,crossed_length_grid=shorter[1],no_uncrossed_path=math.isinf(strict)))
        if result is None:raise ValueError('No disjoint grid route: '+p['rid']+'/'+p['mid'])
        if p.get('straight_locked'):
            if abs(result[1]-math.dist(xy[p['source']],xy[p['target']]))>EPS:raise ValueError('Blocked straight exchange corridor: '+p['rid'])
        p['points']=result[0]
        if any(not integer(a[1] if abs(a[1]-b[1])<EPS else a[0]) for a,b in zip(p['points'],p['points'][1:])):p['auxiliary_lane']=True
        grid.reserve(p)
        if i%30==0:progress(f'Adaptive grid: routed {i+1}/{len(ordered)}',flush=True)
    return exceptions


def auxiliary_channel_load(paths,xy=None):
    channels=defaultdict(list)
    for p in paths:
        if p['priority']!=1:continue
        for a,b in zip(p['points'],p['points'][1:]):
            axis=0 if abs(a[1]-b[1])<EPS else 1
            level=a[1-axis];low,high=sorted((a[axis],b[axis]))
            channels[axis,math.floor(level+EPS)].append((low,high,(p['rid'],p['mid']),not integer(level),p,level,axis))
    peak=0
    for segments in channels.values():
        if not any(entry[3] for entry in segments):continue
        events={v for low,high,*_ in segments for v in (low,high)}
        if xy:
            for low,high,pid,aux,p,level,axis in segments:
                for key in (p['source'],p['target']):
                    center=xy[key];distance=abs(level-center[1-axis])
                    if distance>PORT_RADIUS:continue
                    extent=math.sqrt(max(0.,PORT_RADIUS**2-distance**2))
                    events.update(max(low,min(high,center[axis]+sign*extent)) for sign in (-1,1))
        levels=sorted(events)
        for low,high in zip(levels,levels[1:]):
            if high-low<EPS:continue
            mid=(low+high)/2
            active=[entry for entry in segments if entry[0]+EPS<mid<entry[1]-EPS]
            users=set()
            for a,b,pid,aux,p,level,axis in active:
                point=[mid,level] if axis==0 else [level,mid]
                local=False
                if xy:
                    for _,_,qid,_,q,qlevel,_ in active:
                        if qid==pid:continue
                        other=[mid,qlevel] if axis==0 else [qlevel,mid]
                        if in_common_port(point,other,p,q,xy):local=True;break
                if not local:users.add(pid)
            peak=max(peak,len(users))
    return peak


def outside_port_contacts(first,second,xy):
    result=[]
    for a,b in intersections(first,second):
        if math.dist(a,b)<EPS:
            if not in_common_port(a,b,first,second,xy):result.append((a,b))
            continue
        axis=0 if abs(a[1]-b[1])<EPS else 1;segments=[sorted((a[axis],b[axis]))]
        common={first['source'],first['target']}&{second['source'],second['target']}
        for key in common:
            c=xy[key]
            perpendicular=abs(a[1-axis]-c[1-axis])
            if perpendicular>PORT_RADIUS+EPS:continue
            extent=math.sqrt(max(0.,PORT_RADIUS**2-perpendicular**2))
            lo,hi=c[axis]-extent,c[axis]+extent;remaining=[]
            for left,right in segments:
                if lo>left+EPS:remaining.append([left,min(lo,right)])
                if hi<right-EPS:remaining.append([max(hi,left),right])
            segments=[(l,r) for l,r in remaining if r-l>EPS]
        for lo,hi in segments:
            u,v=list(a),list(a);u[axis]=lo;v[axis]=hi;result.append((tuple(u),tuple(v)))
    return result


def validate_paths(paths,xy,boundary,region,exceptions):
    validate_boundary(paths,boundary)
    from orange_reposition import orange_overlaps
    raw_overlap=orange_overlaps(paths)['length_grid'];visible=0.;port_nodes=set();port_contacts=0
    metrics=measure(paths)
    if any(outside_port_contacts(p,q,xy) for i,p in enumerate(paths) if p['priority']==2 for q in paths[i+1:] if q['priority']==2 and p['rid']!=q['rid']):raise ValueError('Red/red crossing remains')
    if any(math.dist(a,b)>EPS for i,p in enumerate(paths) if p['priority']==2 for q in paths[i+1:] if q['priority']==2 for a,b in outside_port_contacts(p,q,xy)):
        raise ValueError('Red/red overlap outside node')
    for i,p in enumerate(paths):
        for q in paths[i+1:]:
            if not p['priority'] or not q['priority']:continue
            raw=intersections(p,q);outside=outside_port_contacts(p,q,xy)
            if any(in_common_port(a,b,p,q,xy) and (math.dist(a,b)>EPS or not any(math.dist(a,xy[k])<EPS for k in {p['source'],p['target']}&{q['source'],q['target']})) for a,b in raw):
                port_contacts+=1;port_nodes.update({p['source'],p['target']}&{q['source'],q['target']})
            if p['priority']==q['priority']==1:
                visible+=sum(math.dist(a,b) for a,b in outside)
    if visible>EPS:raise ValueError('Visible orange overlap remains: '+str(visible))
    lane_load=auxiliary_channel_load(paths,xy)
    conflicts=parallel_conflicts(paths)
    if conflicts:raise ValueError('Long close parallel lines remain: '+str(conflicts[:3]))
    for p in paths:
        if p.get('curve_controls'):continue
        for a,b in zip(p['points'],p['points'][1:]):
            if abs(a[0]-b[0])>EPS and abs(a[1]-b[1])>EPS:raise ValueError('Non-orthogonal line: '+p['rid'])
            if abs(a[1]-b[1])<EPS:valid=integer(a[1])
            elif abs(a[0]-b[0])<EPS:valid=integer(a[0])
            else:valid=False
            if not valid:
                level=a[1] if abs(a[1]-b[1])<EPS else a[0]
                valid=p.get('auxiliary_lane') and any(abs(level-math.floor(level)-f)<EPS for f in (.25,.75))
            if not valid and not p.get('local_cofactor'):raise ValueError('Non-grid line: '+p['rid'])
            if segment_clearance(a,b,region)<1-EPS and not incident_exempt(a,b,p,xy,region):raise ValueError('TCA clearance violated: '+p['rid'])
    authorized={(e['reaction_id'],e['metabolite_id']) for e in exceptions};crossings=0
    for p in paths:
        if p['priority']!=1:continue
        for q in paths:
            if q['priority']!=2:continue
            hits=outside_port_contacts(p,q,xy)
            if not hits:continue
            if (p['rid'],p['mid']) not in authorized:raise ValueError('Unaudited red/orange crossing: '+p['rid']+'/'+p['mid'])
            if any(math.dist(a,b)>EPS for a,b in hits):raise ValueError('Red/orange overlap outside node')
            crossings+=len(hits)
    metrics={k:v for k,v in metrics.items() if k!='violations'}
    metrics['raw_red_red']=metrics['red_red'];metrics['red_red']=0
    metrics['strict_red_crossings_absent']=metrics['hard_verified'];metrics['hard_verified']=True
    metrics['raw_red_orange']=metrics['red_orange'];metrics['red_orange']=crossings
    return dict(policy=POLICY,verified=True,metrics=metrics,orange_overlap_grid=0.,raw_orange_overlap_grid=raw_overlap,
                node_port_nodes=sorted(port_nodes),node_port_half_width_grid=PORT_RADIUS,node_port_contacts=port_contacts,
                uniform_node_sizes=True,parallel_spacing_grid=PARALLEL_SPACING,parallel_length_grid=PARALLEL_LENGTH,parallel_conflict_count=0,
                grid_lines_verified=not any(p.get('auxiliary_lane') for p in paths),primary_or_auxiliary_grid_verified=True,
                thin_lanes_used=sum(bool(p.get('auxiliary_lane')) for p in paths),tca_clearance_grid=1.,tca_clearance_verified=True,
                boundary_validation=validate_boundary(paths,boundary),conditional_crossings=exceptions,
                orange_length_grid=orange_length(paths),orange_crossings_optimized=False,red_first=True,
                objective_order=['no overlaps outside common node neighborhoods','no long close parallel runs','TCA clearance','short routes with audited red crossings','bends'],
                thin_lane_limit=None,thin_lane_policy='Short fan-out allowed; close parallel runs of two grids or more forbidden',
                maximum_thin_corridor_load=lane_load,global_optimality_proven=False)


def adaptive_layout(model,nodes,edges,xy,backbones,exchange_nodes,groups,cols,rows,major,progress=print,*,red_time_limit=10.,red_candidates=10):
    original_xy=copy.deepcopy(xy);original=copy.deepcopy(backbones);failures=[]
    # Leave primary grid channels to the left of the vertical red backbone.
    for key in original['fixed_node_offsets']:
        if key not in set((original.get('tca_exclusion',{}).get('region') or {}).get('node_ids',[])):
            original_xy[key][0]+=2
            offset=[original_xy[key][0],original_xy[key][1]-(rows-1)]
            original['fixed_node_offsets'][key]=offset
            for route in original['routes']:
                if key in route.get('node_offsets',{}):route['node_offsets'][key]=offset
    from red_objective import optimize_red_core,preservation_report,relaxation_level
    red_keys={k for p in make_paths(nodes,edges,original_xy,original,major) if p['priority']==2 for k in (p['source'],p['target'])}
    access_groups=[g for g in groups if g.get('internal') in red_keys or g['exchange'] in red_keys]
    access_exchanges={g['exchange'] for g in access_groups}
    access_keys=set(red_keys)
    from exchange_corridors import discover_corridors
    access_records=discover_corridors(model,nodes,access_exchanges,access_groups) if access_exchanges else []
    for g in access_records:
        access_keys.update(g[k] for k in ('exchange','external','internal') if k in g)
        for path in g.get('transport_paths',[]):access_keys.update(path)
    access_keys &= set(nodes)
    def check_exchange_access(positions,trial,red):
        if not access_exchanges:return
        preview=copy.deepcopy(trial);points={k:list(positions[k]) for k in access_keys}
        preview_nodes={k:nodes[k] for k in access_keys}
        preview['adaptive_grid']=True
        preview['routing_boundary']=[0,2*(cols-1),-(rows-1),rows-1]
        preview['red_first']=dict(routes=copy.deepcopy(red or []))
        if red:
            for p in red:
                for k in (p['source'],p['target']):preview['fixed_node_offsets'][k]=[points[k][0],points[k][1]-(rows-1)]
        preview['fixed_node_offsets']={k:v for k,v in preview['fixed_node_offsets'].items() if k in access_keys}
        specification=place_corridors(model,preview_nodes,points,preview,access_exchanges,access_groups,cols,rows)
        preview_edges=[e for e in edges if e['source'] in access_keys and e['target'] in access_keys]
        # make_paths uses edge indices for orange membership; locked access
        # incidences are explicitly marked below instead of reusing those IDs.
        paths=initial_paths(preview_nodes,preview_edges,points,preview,major)
        locked=lock_corridor_paths(paths,specification)
        for p in locked:
            if p['priority']!=2:p['priority']=1
        if red is None:
            return dict(paths=[p for p in locked if p['priority']!=2],positions=points,boundary=preview['routing_boundary'])
        red_pairs={(p['rid'],p['mid']) for p in red}
        paths=copy.deepcopy(red)+[p for p in locked if (p['rid'],p['mid']) not in red_pairs]
        for p in paths:
            if p['priority']==2:p['red_fixed']=True
        keys={k for p in paths for k in (p['source'],p['target'])};points={k:points[k] for k in keys}
        region=preview.get('tca_exclusion',{}).get('region')
        exceptions=route_paths(paths,points,preview['routing_boundary'],region)
        validate_paths(paths,points,preview['routing_boundary'],region,exceptions)
    cores,red_optimization=optimize_red_core(nodes,edges,original_xy,original,major,cols,rows,
                                           red_time_limit,red_candidates,progress,validate_candidate=check_exchange_access)
    reference=cores[0]
    # Exhaust non-red expansion first. Only then try equal-length rerouting,
    # equal-length relocation and the shortest longer red alternatives.
    def relaxation(core):
        return relaxation_level(reference,core)
    alternatives=sorted(cores[1:],key=lambda c:(round(c['length'],7),relaxation(c),
                                              preservation_report(reference,c['positions'],c['routes'])['moved_node_count']))[:7]
    schedule=[(reference,i) for i in range(4)]
    schedule.extend((core,expansion) for core in alternatives for expansion in (2,4))
    for attempt,(core,expansion) in enumerate(schedule):
        trial=copy.deepcopy(core['backbones']);positions=copy.deepcopy(original_xy)
        positions.update(copy.deepcopy(core['positions']))
        for key,p in core['positions'].items():trial['fixed_node_offsets'][key]=[p[0],p[1]-(rows-1)]
        right=math.ceil((cols-1)*expansion/2);down=math.ceil((rows-1)*expansion/2);boundary=[0,cols-1+right,-down,rows-1]
        trial['routing_boundary']=boundary;trial['render_frame']=boundary;trial['adaptive_grid']=True
        region=trial.get('tca_exclusion',{}).get('region')
        if region:
            region['clearance_grid']=1.
            trial['tca_geometry']['clearance_grid']=1.
        ring=set((region or {}).get('node_ids',[]));shift=core['shift']
        fixed=set(trial['fixed_node_offsets']);occupied={tuple(positions[k]) for k in fixed}
        # Spread ordinary nodes rather than merely translating the whole map.
        for k in sorted(positions):
            if k in fixed or nodes[k]['kind']=='currency' or k in exchange_nodes:continue
            p=positions[k];preferred=[round(p[0]+max(0,p[0]-3)*expansion/2),round((rows-1)-(rows-1-p[1])*(1+expansion/2))]
            available=[]
            for distance in range(0,20):
                for dx in range(-distance,distance+1):
                    for dy in set((distance-abs(dx),-distance+abs(dx))):
                        q=[float(preferred[0]+dx),float(preferred[1]+dy)]
                        if not(boundary[0]<q[0]<boundary[1] and boundary[2]<q[1]<boundary[3]):continue
                        if tuple(q) in occupied or perimeter_distance(q,region)<1-EPS:continue
                        if not allowed_points([q],region)[0]:continue
                        available.append(q)
                if available:break
            if not available:break
            positions[k]=min(available,key=lambda q:(math.dist(q,p),q));occupied.add(tuple(positions[k]))
        try:
            # Keep an upper external entry immediately above its red endpoint.
            for ex,chain in trial.get('exchange_connectors',{}).get('exchange_paths',{}).items():
                if len(chain)!=4 or chain[0] not in fixed or chain[2] in fixed or chain[1] in fixed:continue
                anchor=positions[chain[0]]
                if anchor[1]<rows-10 or not all(integer(v) for v in anchor):continue
                for distance in (2,3,4):
                    exterior=[anchor[0]+distance,anchor[1]+1];transport=[anchor[0]+distance,anchor[1]]
                    own={chain[1],chain[2]}
                    if any(math.dist(q,p)<.16-EPS for q in (exterior,transport) for k,p in positions.items() if k not in own and nodes[k]['kind']!='currency'):continue
                    positions[chain[2]]=exterior;positions[chain[1]]=transport;break
            attach_grid_cofactors(nodes,positions,region,boundary,core['routes'])
            red=copy.deepcopy(core['routes'])
            trial['red_first']=dict(routes=copy.deepcopy(red),report=dict(verified=True,policy=POLICY))
            if groups or exchange_nodes:
                trial['exchange_corridors']=place_corridors(model,nodes,positions,trial,exchange_nodes,groups,cols,rows)
            important=make_paths(nodes,edges,positions,trial,major)
            protected={k for p in important if p['priority'] for k in (p['source'],p['target'])}
            protected.update(trial.get('exchange_corridors',{}).get('straight_node_ids',[]));protected.update(exchange_nodes)
            occupied={tuple(positions[k]) for k in protected}
            for key in sorted(positions):
                if key in protected or nodes[key]['kind']=='currency' or key in exchange_nodes:continue
                base=positions[key];options=[]
                for distance in range(0,40):
                    for dx in range(-distance,distance+1):
                        for dy in set((distance-abs(dx),-distance+abs(dx))):
                            q=[float(round(base[0])+dx),float(round(base[1])+dy)]
                            if not(boundary[0]<q[0]<boundary[1] and boundary[2]<q[1]<boundary[3]):continue
                            if tuple(q) in occupied or perimeter_distance(q,region)<1-EPS:continue
                            if any(_on(q,*g['corridor']) for g in trial.get('exchange_corridors',{}).get('groups',[])):continue
                            if any(_on(q,a,b) for p in red if not p.get('curve_controls') for a,b in zip(p['points'],p['points'][1:])):continue
                            if any(math.dist(q,positions[k])<1.01 for k in protected):continue
                            if region and not allowed_points([q],dict(region,interior_policy='forbid'))[0]:continue
                            options.append(q)
                    if options:break
                if not options:raise ValueError('Congested ordinary nodes need more grid space')
                positions[key]=min(options,key=lambda q:(math.dist(q,base),q));occupied.add(tuple(positions[key]))
            # Reattach movable reaction-local currency copies on grid lines.
            attach_grid_cofactors(nodes,positions,region,boundary,red)
            paths=make_paths(nodes,edges,positions,trial,major)
            saved={(p['rid'],p['mid']):p for p in red}
            for p in paths:
                if p['priority']==2:
                    p['points']=copy.deepcopy(saved[p['rid'],p['mid']]['points']);p['red_fixed']=True
                    if saved[p['rid'],p['mid']].get('auxiliary_lane'):p['auxiliary_lane']=True
            if trial.get('exchange_corridors'):lock_corridor_paths(paths,trial['exchange_corridors'])
            exceptions=route_paths(paths,positions,boundary,region,progress)
            report=validate_paths(paths,positions,boundary,region,exceptions)
            report['node_grid']=validate_node_grid(nodes,positions,region)
            report['red_routes_preserved']=all(p['points']==saved[p['rid'],p['mid']]['points'] for p in paths if p['priority']==2)
            if not preservation_report(core,positions,paths)['exactly_preserved']:raise ValueError('Frozen red core changed during non-red placement')
            if trial.get('exchange_corridors'):report['straight_exchange_validation']=validate_corridors(paths,trial['exchange_corridors'])
            validate_exclusion(nodes,positions,region)
            signature=hashlib.sha256(json.dumps(paths,sort_keys=True).encode()).hexdigest()[:24]
            trial['crossing_routing']=dict(routes=paths,report=report,signature=signature)
            trial['adaptive_layout']=dict(policy=POLICY,expansion_steps=expansion,before_boundary=[0,cols-1,0,rows-1],after_boundary=boundary,
                tca_translation=shift,rejected_layouts=failures,congestion_rule='Short fan-out allowed; close parallel runs of two grids or more forbidden',thin_lanes_used=report['thin_lanes_used'])
        except (ValueError,RuntimeError) as error:
            failures.append(dict(attempt=attempt,expansion_steps=expansion,red_length_grid=core['length'],relaxation_level=relaxation(core),boundary=boundary,reason=str(error)[:400]))
            progress(f'Adaptive grid expansion {attempt}: {error}',flush=True);continue
        xy.clear();xy.update(positions);backbones.clear();backbones.update(trial)
        compact_unused_space(model,nodes,edges,xy,backbones,exchange_nodes,groups,cols,rows,major,progress)
        final_paths=backbones['crossing_routing']['routes']
        if not preservation_report(core,xy,final_paths)['exactly_preserved']:raise ValueError('Compaction changed the frozen red core')
        red_optimization.update(preservation=preservation_report(reference,xy,final_paths),relaxation_level=relaxation(core),
                                relaxation_reason=None if core is reference else 'Previous red cores failed full placement; see rejected_layouts',
                                rejected_full_layouts=failures,accepted_reference_signature=core['signature'])
        backbones['red_optimization']=red_optimization
        backbones['crossing_routing']['report']['red_preservation']=red_optimization['preservation']
        return backbones['adaptive_layout']
    raise ValueError('Could not satisfy adaptive grid layout after expansion: '+str(failures))


def compact_unused_space(model,nodes,edges,xy,backbones,exchange_nodes,groups,cols,rows,major,progress=print):
    """Retain expanded node spacing, but bring remote EX edges back toward it."""
    if compact_terminal_corridors(nodes,xy,backbones,progress):return
    original=backbones['routing_boundary'];baseline=backbones['crossing_routing']['report']['orange_length_grid']
    corridor_nodes=set();internal=set()
    for g in backbones.get('exchange_corridors',{}).get('groups',[]):
        corridor_nodes.update(g['node_ids'])
        if 'internal' in g:internal.add(g['internal'])
    excluded=(corridor_nodes-internal)|set(exchange_nodes)
    excluded.update(k for k,n in nodes.items() if n['kind']=='currency' and n['parent_id'] in excluded)
    content=[p for k,p in xy.items() if k not in excluded]
    if not content:return
    rejections=[]
    for margin in (1,2,3,5,8):
        boundary=[original[0],min(original[1],math.ceil(max(p[0] for p in content))+margin),max(original[2],math.floor(min(p[1] for p in content))-margin),original[3]]
        if boundary==original:continue
        trial=copy.deepcopy(backbones);positions=copy.deepcopy(xy);trial['routing_boundary']=boundary
        region=trial['tca_exclusion']['region'];red=trial['red_first']['routes']
        try:
            trial['exchange_corridors']=place_corridors(model,nodes,positions,trial,exchange_nodes,groups,cols,rows)
            candidates=make_paths(nodes,edges,positions,trial,major)
            protected={k for p in candidates if p['priority'] for k in (p['source'],p['target'])}
            protected.update(trial['exchange_corridors']['straight_node_ids']);protected.update(exchange_nodes)
            occupied={tuple(positions[k]) for k in protected}
            for key in sorted(positions):
                if key in protected or nodes[key]['kind']=='currency':continue
                base=positions[key];options=[]
                for radius in range(40):
                    for dx in range(-radius,radius+1):
                        for dy in set((radius-abs(dx),-radius+abs(dx))):
                            q=[float(round(base[0])+dx),float(round(base[1])+dy)]
                            if not(boundary[0]<q[0]<boundary[1] and boundary[2]<q[1]<boundary[3]):continue
                            if tuple(q) in occupied or perimeter_distance(q,region)<1-EPS:continue
                            if any(_on(q,*g['corridor']) for g in trial['exchange_corridors']['groups']):continue
                            if any(_on(q,a,b) for p in red if not p.get('curve_controls') for a,b in zip(p['points'],p['points'][1:])):continue
                            if any(math.dist(q,positions[k])<1.01 for k in protected):continue
                            if region and not allowed_points([q],dict(region,interior_policy='forbid'))[0]:continue
                            options.append(q)
                    if options:break
                if not options:raise ValueError('No compact position for gray node')
                positions[key]=min(options,key=lambda q:(math.dist(q,base),q));occupied.add(tuple(positions[key]))
            attach_grid_cofactors(nodes,positions,region,boundary,red)
            paths=make_paths(nodes,edges,positions,trial,major);saved={(p['rid'],p['mid']):p for p in red}
            for p in paths:
                if p['priority']==2:
                    p['points']=copy.deepcopy(saved[p['rid'],p['mid']]['points']);p['red_fixed']=True
                    if saved[p['rid'],p['mid']].get('auxiliary_lane'):p['auxiliary_lane']=True
            lock_corridor_paths(paths,trial['exchange_corridors'])
            exceptions=route_paths(paths,positions,boundary,region)
            report=validate_paths(paths,positions,boundary,region,exceptions)
            report['node_grid']=validate_node_grid(nodes,positions,region)
            report['straight_exchange_validation']=validate_corridors(paths,trial['exchange_corridors']);report['red_routes_preserved']=True
            validate_exclusion(nodes,positions,region)
            if report['orange_length_grid']>baseline+EPS:raise ValueError('Compaction would lengthen orange routes')
            trial['crossing_routing']=dict(routes=paths,report=report,signature=hashlib.sha256(json.dumps(paths,sort_keys=True).encode()).hexdigest()[:24])
            trial['render_frame']=boundary
            trial['adaptive_layout'].update(after_boundary=boundary,expanded_boundary=original,unused_margin_compacted=True,
                compaction_rejections=rejections,thin_lanes_used=report['thin_lanes_used'])
        except (ValueError,RuntimeError) as error:
            rejections.append(dict(boundary=boundary,reason=str(error)[:200]));continue
        xy.clear();xy.update(positions);backbones.clear();backbones.update(trial)
        progress(f'Adaptive unused margins: {original} -> {boundary}; orange length {baseline:.2f} -> {report["orange_length_grid"]:.2f}',flush=True)
        return
    backbones['adaptive_layout'].update(unused_margin_compacted=False,compaction_rejections=rejections)


def compact_terminal_corridors(nodes,xy,backbones,progress=print):
    """Shorten exterior exchange tails without disturbing validated interior routes."""
    specification=backbones.get('exchange_corridors')
    if not specification:return False
    excluded=set()
    for g in specification['groups']:
        if g['face'] not in ('right','bottom'):continue
        excluded.update(k for k in g['node_ids'] if k!=g.get('internal'))
    excluded.update(k for k,n in nodes.items() if n['kind']=='currency' and n['parent_id'] in excluded)
    original=backbones['routing_boundary'];plan=backbones['crossing_routing']
    content=[p for k,p in xy.items() if k not in excluded]
    content.extend(q for p in plan['routes'] if p['source'] not in excluded and p['target'] not in excluded for q in p['points'])
    for margin in (1,2,3):
        boundary=[original[0],min(original[1],math.ceil(max(p[0] for p in content))+margin),
                  max(original[2],math.floor(min(p[1] for p in content))-margin),original[3]]
        if boundary==original:return False
        trial=copy.deepcopy(backbones);positions=copy.deepcopy(xy)
        try:
            for g in trial['exchange_corridors']['groups']:
                if g['face'] not in ('right','bottom'):continue
                axis,edge,sign=(0,boundary[1],-1) if g['face']=='right' else (1,boundary[2],1)
                chain=g.get('chain_node_ids',[])[:]
                if chain and g.get('internal'):chain=chain[:-1]
                if not chain:
                    chain=[g['exchange'],g['external']]
                    if g['primary_transport']:chain.append('rxn:'+g['primary_transport'])
                for offset,key in enumerate(chain):positions[key][axis]=edge+sign*offset
                target=g.get('internal',g['external'])
                g['corridor']=[list(positions[g['exchange']]),list(positions[target])]
                if g.get('internal') and sign*(positions[target][axis]-positions[chain[-1]][axis])<1-EPS:raise ValueError('Insufficient exchange tail space')
            for key,n in nodes.items():
                if n['kind']=='currency':
                    parent=n['parent_id']
                    positions[key]=[xy[key][a]+positions[parent][a]-xy[parent][a] for a in (0,1)]
            fixed=trial['fixed_node_offsets']
            if any(math.dist(positions[k],xy[k])>EPS for k in fixed):raise ValueError('Fixed exchange cannot move')
            paths=trial['crossing_routing']['routes']
            for p in paths:
                if not any(math.dist(positions[k],xy[k])>EPS for k in (p['source'],p['target'])):continue
                a,b=positions[p['source']],positions[p['target']]
                if p.get('straight_locked'):p['points']=clean([a,b])
                elif p.get('local_cofactor'):p['points']=clean([a,[b[0],a[1]],b])
                else:
                    for key,reverse in ((p['source'],False),(p['target'],True)):
                        if math.dist(positions[key],xy[key])<EPS:continue
                        points=list(reversed(p['points'])) if reverse else p['points'];q=positions[key]
                        index=next((i for i,(u,v) in enumerate(zip(points,points[1:])) if _on(q,u,v)),None)
                        if index is None:
                            if len(points)<3:raise ValueError('Exchange tail has another incidence: '+p['rid']+'/'+p['mid'])
                            elbow=([points[1][0],q[1]] if abs(points[0][1]-points[1][1])<EPS
                                   else [q[0],points[1][1]])
                            points=clean([q,elbow,*points[2:]])
                        else:points=clean([q,*points[index+1:]])
                        p['points']=list(reversed(points)) if reverse else points
            region=trial.get('tca_exclusion',{}).get('region')
            exceptions=trial['crossing_routing']['report']['conditional_crossings']
            report=validate_paths(paths,positions,boundary,region,exceptions)
            report['node_grid']=validate_node_grid(nodes,positions,region)
            report['straight_exchange_validation']=validate_corridors(paths,trial['exchange_corridors'])
            report['red_routes_preserved']=True
            validate_exclusion(nodes,positions,region)
            if any(not(boundary[0]-EPS<=p[0]<=boundary[1]+EPS and boundary[2]-EPS<=p[1]<=boundary[3]+EPS) for p in positions.values()):raise ValueError('Node outside compact frame')
            trial['crossing_routing'].update(report=report,signature=hashlib.sha256(json.dumps(paths,sort_keys=True).encode()).hexdigest()[:24])
            trial['routing_boundary']=boundary;trial['render_frame']=boundary
            trial['adaptive_layout'].update(after_boundary=boundary,expanded_boundary=original,unused_margin_compacted=True,compaction_method='shorten exterior exchange tails')
        except (ValueError,RuntimeError) as error:
            progress(f'Adaptive terminal compaction rejected: {error}',flush=True)
            continue
        xy.clear();xy.update(positions);backbones.clear();backbones.update(trial)
        progress(f'Adaptive terminal compaction: {original} -> {boundary}',flush=True)
        return True
    return False
