"""Keep structural red/orange paths simple; lower-priority gray paths yield."""
import math
from collections import defaultdict
from parallel_routing import EPS,_legs,_overlap,_distance,_node_detour

POLICY='important_bends_then_length_gray_yields'


def clean(points):
    out=[]
    for point in points:
        point=tuple(point)
        if out and math.dist(out[-1],point)<EPS:continue
        while len(out)>1:
            a,b=out[-2:]
            aligned=(abs(a[0]-b[0])<EPS and abs(b[0]-point[0])<EPS or
                     abs(a[1]-b[1])<EPS and abs(b[1]-point[1])<EPS)
            if not aligned or sum((b[k]-a[k])*(point[k]-b[k]) for k in (0,1))<0:break
            out.pop()
        out.append(point)
    return out


def candidates(path,offset,grid):
    a,b=tuple(path['points'][0]),tuple(path['points'][-1]);seen=set()
    def remember(points):
        points=clean(points);key=tuple(points)
        if key in seen:return None
        seen.add(key);return points
    raw=[path['points'],[a,(a[0],b[1]),b],[a,(b[0],a[1]),b]]
    xs={a[0],b[0],(a[0]+b[0])/2,a[0]+offset,a[0]-offset,b[0]+offset,b[0]-offset}
    ys={a[1],b[1],(a[1]+b[1])/2,a[1]+offset,a[1]-offset,b[1]+offset,b[1]-offset}
    for x in sorted(xs):raw.append([a,(x,a[1]),(x,b[1]),b])
    for y in sorted(ys):raw.append([a,(a[0],y),(b[0],y),b])
    if abs(a[0]-b[0])<EPS or abs(a[1]-b[1])<EPS:
        leg=next(_legs([a,b]),None)
        if leg:
            for inset in [(True,False),(False,True),(True,True)]:
                for sign in (1,-1):raw.append(_node_detour(leg,sign*offset,grid,*inset))
    for points in raw:
        result=remember(points)
        if result is not None:yield result


def length(points):
    return sum(math.dist(a,b) for a,b in zip(points,points[1:]))


def node_hits(path,points,obstacles,grid):
    endpoints={tuple(path['points'][0]),tuple(path['points'][-1])}
    return [p for p in obstacles if p not in endpoints and any(_distance(p,s)<.12*grid-EPS for s in _legs(points))]



def line_index(paths):
    index=defaultdict(list)
    for path in paths:
        for leg in _legs(path['points']):
            index[(leg['axis'],math.floor(leg['level']/EPS))].append(dict(leg,rid=path['rid']))
    return index


def coincident(leg,index,rid):
    bucket=math.floor(leg['level']/EPS)
    for k in (bucket-1,bucket,bucket+1):
        for other in index.get((leg['axis'],k),()):
            if other['rid']!=rid and abs(leg['level']-other['level'])<EPS and _overlap(leg,other)>EPS:yield other


def interactions(path,points,others,grid,red_spacing,gray_index=None):
    """Important coincidences are hard; gray interference is a soft cost."""
    crossings=0;gray_cost=0
    original=list(_legs(path['points']));endpoints={tuple(path['points'][0]),tuple(path['points'][-1])}
    for segment in _legs(points):
        if gray_index is not None:
            gray_cost+=sum(_overlap(segment,current)/grid for current in coincident(segment,gray_index,path['rid']))
        for other in others:
            if path is other or path['rid']==other['rid']:continue
            for current in _legs(other['points']):
                if segment['axis']==current['axis']:
                    overlap=_overlap(segment,current);gap=abs(segment['level']-current['level'])
                    if overlap<=EPS:continue
                    if other['priority']==0:
                        if gap<EPS:gray_cost+=overlap/grid
                        continue
                    required=red_spacing if path['priority']==other['priority']==2 else .3
                    if overlap>=grid-EPS and gap<required*grid-EPS:return None
                    if gap<EPS:
                        retained=max((_overlap(old,current) for old in original if old['axis']==current['axis'] and abs(old['level']-current['level'])<EPS),default=0)
                        at_endpoint=any(_distance(p,segment)<EPS and _distance(p,current)<EPS for p in endpoints)
                        if overlap>.15*grid+EPS or overlap>retained+EPS or not at_endpoint:return None
                else:
                    h,v=(segment,current) if segment['axis']==0 else (current,segment)
                    if h['lo']+EPS<v['level']<h['hi']-EPS and v['lo']+EPS<h['level']<v['hi']-EPS:
                        if other['priority']>0:crossings+=1
    return crossings,gray_cost


def shared_trunk(path,others):
    legs=list(_legs(path['points']))
    for other in others:
        if other is path or other['rid']!=path['rid'] or other['priority']!=path['priority']:continue
        for a in legs:
            for b in _legs(other['points']):
                if a['axis']==b['axis'] and abs(a['level']-b['level'])<EPS and _overlap(a,b)>EPS:return True
    return False


def simplify_important(paths,obstacles=(),grid=1.0,route_offset=.3,red_spacing=.3):
    obstacles=[tuple(p) for p in obstacles]
    report=dict(policy=POLICY,changed_paths=0,removed_bends=0,details=[],unresolved=[])
    if route_offset==0:return report
    important=[p for p in paths if p['priority']>0]
    gray_index=line_index([p for p in paths if p['priority']==0])
    for path in sorted(paths,key=lambda p:(-p['priority'],p['rid'],p.get('mid',''))):
        if not path['priority'] or path.get('routing_locked'):continue
        original=clean(path['points']);path['points']=original
        # Keep positively shared reaction trunks bundled rather than rerouting one copy.
        if shared_trunk(path,paths):continue
        viable=[]
        for index,points in enumerate(candidates(path,route_offset*grid,grid)):
            if node_hits(path,points,obstacles,grid):continue
            cost=interactions(path,points,important,grid,red_spacing,gray_index)
            if cost is not None:
                viable.append(((len(points)-2,round(length(points)/grid,8),*cost,index),points))
        if not viable:
            report['unresolved'].append(dict(reaction_id=path['rid'],reason='No simpler feasible important route'))
            continue
        _,chosen=min(viable,key=lambda c:c[0])
        if chosen!=original:
            removed=len(original)-len(chosen)
            path['points']=chosen;report['changed_paths']+=1;report['removed_bends']+=removed
            report['details'].append(dict(reaction_id=path['rid'],bends_before=len(original)-2,bends_after=len(chosen)-2))
    return report


def gray_conflicts(path,points,important,index=None):
    index=line_index(important) if index is None else index
    return sum(1 for leg in _legs(points) for _ in coincident(leg,index,path['rid']))


def yield_gray(paths,obstacles=(),grid=1.0,route_offset=.3):
    """Move only conflicting gray routes; red/orange geometry stays untouched."""
    important=[p for p in paths if p['priority']>0];obstacles=[tuple(p) for p in obstacles]
    important_index=line_index(important)
    report=dict(changed_paths=0,details=[],unresolved=[])
    if route_offset==0:return report
    for path in sorted(paths,key=lambda p:(p['rid'],p.get('mid',''))):
        if path['priority'] or not gray_conflicts(path,path['points'],important,important_index):continue
        original=clean(path['points']);path['points']=original
        if shared_trunk(path,paths):
            report['unresolved'].append(dict(reaction_id=path['rid'],reason='Shared gray trunk retained'));continue
        viable=[]
        for index,points in enumerate(candidates(path,max(.15,route_offset/2)*grid,grid)):
            if node_hits(path,points,obstacles,grid) or gray_conflicts(path,points,important,important_index):continue
            viable.append(((len(points)-2,round(length(points)/grid,8),index),points))
        if viable:
            _,chosen=min(viable,key=lambda c:c[0]);path['points']=chosen
            if chosen!=original:
                report['changed_paths']+=1;report['details'].append(dict(reaction_id=path['rid'],bends_before=len(original)-2,bends_after=len(chosen)-2))
        else:report['unresolved'].append(dict(reaction_id=path['rid'],reason='No free gray route found'))
    return report
