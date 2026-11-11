"""Allocate red routing channels during layout, before structural rendering."""
import hashlib
import json
import math
import time
from collections import defaultdict

from parallel_routing import _legs, _overlap, _distance, _node_detour


def validate_red_spacing(value):
    if not math.isfinite(value) or value < 0:
        raise ValueError('red-parallel-spacing must be finite and nonnegative')
    return value


def spacing_report(paths, spacing, grid=1.0):
    """Check all actual same-axis red segments, excluding same-reaction trunks."""
    segments=[dict(leg,rid=p['rid']) for p in paths if p.get('priority',2)==2 for leg in _legs(p['points'])]
    failures=[];minimum=None;count=0
    for i,a in enumerate(segments):
        for b in segments[i+1:]:
            if a['rid']==b['rid'] or a['axis']!=b['axis'] or _overlap(a,b)<grid-1e-7:continue
            gap=abs(a['level']-b['level'])/grid
            minimum=gap if minimum is None else min(minimum,gap);count+=1
            if gap<spacing-1e-7:
                failures.append(dict(reactions=[a['rid'],b['rid']],axis=a['axis'],gap_grid=gap,
                                     overlap_grid=_overlap(a,b)/grid))
    return dict(requested_spacing_grid=spacing,min_parallel_length_grid=1,
                minimum_parallel_gap_grid=minimum,parallel_pairs=count,violations=failures,
                verified=not failures)


def plan_from_metadata(meta):
    for key in ('sbml2map','sbml2escher','escher_grid'):
        plan=meta.get(key,{}).get('major_backbones',{}).get('red_parallel_layout')
        if plan:return plan
    return None


def plan_red_routes(xy, incidences, spacing=1.0, obstacle_node_ids=(), time_limit=3.0, use_milp=True):
    """Choose compatible channels with a bounded greedy / finite-candidate MILP.

    Fixed biological positions are retained. Long red shafts reserve their full
    original spans, so shortening terminal stubs cannot evade the spacing rule.
    Candidate geometry also checks perpendicular shoulders and protected nodes.
    """
    started=time.monotonic();validate_red_spacing(spacing)
    protected=sorted(set(obstacle_node_ids)&set(xy))
    report=dict(enabled=spacing>0,requested_spacing_grid=spacing,min_parallel_length_grid=1,
                stage='grid layout routing-channel allocation',obstacle_node_ids=protected,
                routes=[],allocations=[],unresolved=[],nodes_moved=False,global_optimality_proven=False)
    if not spacing or not incidences:
        report.update(verified=True,elapsed_seconds=time.monotonic()-started,method='Disabled or no red grid incidences')
        return report
    if not math.isfinite(spacing*(len(incidences)+3)):
        raise ValueError('red-parallel-spacing is too large for routing channels')
    deadline=started+max(.05,time_limit)
    paths=[]
    for edge in sorted(incidences,key=lambda e:(e['rid'],e['mid'])):
        a,b=tuple(xy[edge['source']]),tuple(xy[edge['target']])
        if edge.get('reverse'):a,b=b,a
        elbow=(b[0],a[1]);points=[a]+([elbow] if elbow!=a and elbow!=b else [])+[b]
        if edge.get('reverse'):points.reverse()
        paths.append(dict(edge,points=points,priority=2))
    buckets=defaultdict(list);legs=[]
    for index,path in enumerate(paths):
        for part,leg in enumerate(_legs(path['points'])):
            leg.update(path=index,part=part,rid=path['rid'],backbone=path.get('backbone',False),
                       endpoints={tuple(path['points'][0]),tuple(path['points'][-1])})
            buckets[(leg['rid'],leg['axis'],round(leg['level'],8))].append(leg);legs.append(leg)
    groups=[]
    for key,bucket in sorted(buckets.items()):
        group=[];end=-math.inf
        for leg in sorted(bucket,key=lambda a:(a['lo'],a['hi'])):
            if group and leg['lo']>=end-1e-7:groups.append(group);group=[]
            group.append(leg);end=max(end,leg['hi'])
        if group:groups.append(group)
    groups.sort(key=lambda g:(not any(l['backbone'] for l in g),-max(l['hi']-l['lo'] for l in g),g[0]['rid'],g[0]['axis'],g[0]['level']))
    obstacles=[tuple(xy[k]) for k in protected]
    def make_candidate(group,offset,mode):
        shapes=[_node_detour(leg,offset,1,*mode) if offset else [leg['a'],leg['b']] for leg in group]
        for leg,points in zip(group,shapes):
            if any(_distance(node,s)<.12-1e-7 for node in obstacles if node not in leg['endpoints'] for s in _legs(points)):return None
        actual=[dict(s,rid=group[0]['rid']) for points in shapes for s in _legs(points)]
        reserved=[dict(leg,level=leg['level']+offset) for leg in group]
        # Keep backbone centerlines when possible; routed side branches take channels first.
        weight=100 if any(l['backbone'] for l in group) else 1
        return dict(offset=offset,mode=mode,shapes=shapes,actual=actual,reserved=reserved,
                    cost=weight*abs(offset)+.001*(sum(len(p)-2 for p in shapes)))
    def compatible(a,b):
        if time.monotonic()>deadline:return False
        if a['reserved'][0]['rid']==b['reserved'][0]['rid']:return True
        for first,second in ((a['reserved'],b['reserved']),(a['actual'],b['actual'])):
            for x in first:
                for y in second:
                    if x['axis']==y['axis'] and _overlap(x,y)>=1-1e-7 and abs(x['level']-y['level'])<spacing-1e-7:return False
        return True
    step=max(1,math.ceil(spacing));domains=[];timed_out=False
    for group in groups:
        candidates=[]
        # Outside channels expand the routing envelope, without moving EX nodes.
        for offset in [0]+[sign*step*k for k in range(1,len(groups)+3) for sign in (1,-1)]:
            if time.monotonic()>deadline:timed_out=True;break
            for mode in ([(False,False)] if not offset else [(False,False),(True,False),(False,True),(True,True)]):
                c=make_candidate(group,offset,mode)
                if c is not None:candidates.append(c)
            if len(candidates)>=20:break
        candidates.sort(key=lambda c:(c['cost'],c['offset']<0,c['mode']))
        domains.append(candidates)
        if not candidates:break
    chosen=None;method='Constraint-checked greedy routing channels'
    if len(domains)==len(groups) and all(domains):
        for ordering in (list(range(len(groups))),sorted(range(len(groups)),key=lambda i:len(domains[i]))):
            trial={}
            for i in ordering:
                if time.monotonic()>deadline:break
                found=next((c for c in domains[i] if all(compatible(c,b) for b in trial.values())),None)
                if found is None:break
                trial[i]=found
            if len(trial)==len(groups):chosen=trial;break
    if chosen is None and use_milp and not timed_out and len(domains)==len(groups) and all(domains) and time.monotonic()<deadline:
        import numpy as np
        from scipy.optimize import milp,Bounds,LinearConstraint
        from scipy.sparse import coo_matrix
        variables=[(i,j) for i,d in enumerate(domains) for j in range(len(d))]
        lookup={pair:k for k,pair in enumerate(variables)};cuts=set()
        for _ in range(20):
            remaining=deadline-time.monotonic()
            if remaining<=.01:break
            rr=[];cc=[];vv=[];low=[];high=[]
            for i,domain in enumerate(domains):
                row=len(low);low.append(1);high.append(1)
                for j in range(len(domain)):rr.append(row);cc.append(lookup[i,j]);vv.append(1)
            for a,b in sorted(cuts):
                row=len(low);low.append(-np.inf);high.append(1)
                rr.extend([row,row]);cc.extend([a,b]);vv.extend([1,1])
            matrix=coo_matrix((vv,(rr,cc)),shape=(len(low),len(variables))).tocsc()
            result=milp(np.array([domains[i][j]['cost'] for i,j in variables]),integrality=np.ones(len(variables)),
                        bounds=Bounds(0,1),constraints=LinearConstraint(matrix,low,high),options={'time_limit':remaining})
            if result.x is None:break
            selected={i:(j,k) for k,(i,j) in enumerate(variables) if result.x[k]>.5}
            if len(selected)!=len(groups):break
            conflicts=[]
            for i in range(len(groups)):
                for j in range(i+1,len(groups)):
                    ai,ak=selected[i];bi,bk=selected[j]
                    if not compatible(domains[i][ai],domains[j][bi]):conflicts.append((min(ak,bk),max(ak,bk)))
            if not conflicts:
                chosen={i:domains[i][selected[i][0]] for i in range(len(groups))};method='Finite-candidate routing MILP / HiGHS';break
            cuts.update(conflicts)
    if chosen is None:
        original=spacing_report(paths,spacing)
        report.update(original)
        report.update(verified=False,method=method,elapsed_seconds=time.monotonic()-started,
                      unresolved=[dict(reason='No compatible red-channel assignment found within the layout time budget')])
        return report
    by_path=defaultdict(list)
    for index,group in enumerate(groups):
        candidate=chosen[index]
        report['allocations'].append(dict(rid=group[0]['rid'],axis=group[0]['axis'],original_level=group[0]['level'],
            lane_level=group[0]['level']+candidate['offset'],offset=candidate['offset'],inset_start=candidate['mode'][0],inset_end=candidate['mode'][1]))
        for leg,points in zip(group,candidate['shapes']):by_path[leg['path']].append((leg['part'],points))
    for index,path in enumerate(paths):
        points=[]
        for _,part in sorted(by_path[index]):points.extend(part if not points else part[1:])
        report['routes'].append({k:path[k] for k in ('rid','mid','source','target')}|dict(points=points,priority=2))
    check=spacing_report(report['routes'],spacing)
    report.update(check);report['verified']=check['verified'];report['method']=method
    report['elapsed_seconds']=time.monotonic()-started
    report['signature']=hashlib.sha256(json.dumps(report['routes'],sort_keys=True).encode()).hexdigest()[:24]
    points=[p for route in report['routes'] for p in route['points']]
    report['routing_envelope']=dict(x_min=min(p[0] for p in points),x_max=max(p[0] for p in points),
                                    y_min=min(p[1] for p in points),y_max=max(p[1] for p in points))
    return report
