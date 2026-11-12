"""Length-first optimization and immutable references for the red core."""
import copy
import hashlib
import json
import math
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from crossing_routing import EPS
from tca_geometry import allowed_points


def route_length(path):
    """Drawn length, including the actual cubic rather than its chord."""
    if not path.get('curve_controls'):
        return sum(math.dist(a,b) for a,b in zip(path['points'],path['points'][1:]))
    p=np.asarray([path['points'][0],*path['curve_controls'],path['points'][-1]],dtype=float)
    points,weights=np.polynomial.legendre.leggauss(16)
    t=(points+1)/2
    derivative=(3*(1-t[:,None])**2*(p[1]-p[0])+
                6*(1-t[:,None])*t[:,None]*(p[2]-p[1])+
                3*t[:,None]**2*(p[3]-p[2]))
    return float(np.dot(weights,np.linalg.norm(derivative,axis=1))/2)


def red_length(paths):
    # Each biological incidence counts once, including its short shared port.
    return sum(route_length(p) for p in paths if p['priority']==2)


def translation_seed(paths,xy,ring,bounds,time_limit):
    """MILP lower-bound seed for translating a rigid TCA against the backbone."""
    terms=[]
    for p in paths:
        if p.get('curve_controls'):continue
        a,b=p['source'],p['target'];coefficient=int(a in ring)-int(b in ring)
        if not coefficient:continue
        for axis in (0,1):terms.append((axis,coefficient,xy[a][axis]-xy[b][axis]))
    if not terms:return (0,0),dict(method='No movable red attachment',optimality_proven=True)
    size=2+len(terms);cost=np.r_[np.zeros(2),np.ones(len(terms))]
    lower=np.r_[[b[0] for b in bounds],np.zeros(len(terms))]
    upper=np.r_[[b[1] for b in bounds],np.full(len(terms),np.inf)]
    matrix=[];limits=[]
    for i,(axis,coefficient,constant) in enumerate(terms):
        for sign in (-1,1):
            row=np.zeros(size);row[axis]=sign*coefficient;row[2+i]=-1
            matrix.append(row);limits.append(-sign*constant)
    result=milp(cost,integrality=np.r_[np.ones(2),np.zeros(len(terms))],
                bounds=Bounds(lower,upper),
                constraints=LinearConstraint(matrix,-np.inf,limits),
                options=dict(time_limit=max(.05,min(3.,time_limit)),mip_rel_gap=0.))
    seed=(tuple(int(round(v)) for v in result.x[:2]) if result.x is not None else (0,0))
    return seed,dict(method='MILP rigid-TCA translation, nonarc L1 attachment objective',
                     optimality_proven=result.status==0,status=int(result.status),
                     objective=None if result.fun is None else float(result.fun),
                     limitation='Lower-bound seed only; routed red optimality is not implied.')


def translated_core(xy,backbones,ring,shift,rows,*,global_dx=0,red_keys=()):
    positions=copy.deepcopy(xy);trial=copy.deepcopy(backbones)
    for k in red_keys:
        positions[k][0]+=global_dx
        trial['fixed_node_offsets'][k]=[positions[k][0],positions[k][1]-(rows-1)]
    for k in ring:
        positions[k]=[xy[k][a]+shift[a]+(global_dx if a==0 else 0) for a in (0,1)]
        trial['fixed_node_offsets'][k]=[positions[k][0],positions[k][1]-(rows-1)]
    for route in trial['routes']:
        for k in ring|set(red_keys):
            if k in route.get('node_offsets',{}):route['node_offsets'][k]=trial['fixed_node_offsets'][k][:]
    region=trial.get('tca_exclusion',{}).get('region')
    if region:
        geometry=trial['tca_geometry']
        if region['shape']=='circle':
            region['center']=[region['center'][a]+shift[a]+(global_dx if a==0 else 0) for a in (0,1)]
            offset=[region['center'][0],region['center'][1]-(rows-1)]
            region['center_offset']=offset[:];geometry['center_offset']=offset[:]
            trial['circle']['center_offset']=offset[:]
        else:
            region['left']+=shift[0]+global_dx;region['top']+=shift[1]
            geometry.update(left=region['left'],top_offset=region['top']-(rows-1))
            trial['square'].update(left=region['left'],top_offset=region['top']-(rows-1))
    return positions,trial


def reference_signature(positions,paths):
    return hashlib.sha256(json.dumps(dict(positions=positions,routes=paths),sort_keys=True).encode()).hexdigest()[:24]


def preservation_report(reference,positions,paths):
    actual=[p for p in paths if p['priority']==2]
    saved={(p['rid'],p['mid']):p for p in reference['routes']}
    if {(p['rid'],p['mid']) for p in actual}!=set(saved):raise ValueError('Red incidence set changed')
    moved=[k for k,p in reference['positions'].items() if math.dist(p,positions[k])>EPS]
    def same_points(a,b):
        return len(a)==len(b) and all(math.dist(p,q)<=EPS for p,q in zip(a,b))
    changed=[(p['rid'],p['mid']) for p in actual
             if not same_points(p['points'],saved[p['rid'],p['mid']]['points']) or
             not same_points(p.get('curve_controls',[]),saved[p['rid'],p['mid']].get('curve_controls',[]))]
    length=red_length(actual);baseline=reference['length']
    return dict(initial_length_grid=baseline,final_length_grid=length,length_increase_grid=length-baseline,
                moved_node_count=len(moved),moved_node_ids=moved,changed_route_count=len(changed),
                changed_routes=changed,exactly_preserved=not moved and not changed,
                reference_signature=reference['signature'])


def relaxation_level(reference,candidate):
    check=preservation_report(reference,candidate['positions'],candidate['routes'])
    return (3 if check['length_increase_grid']>EPS else 2 if check['moved_node_count']
            else 1 if check['changed_route_count'] else 0)


def optimize_red_core(nodes,edges,xy,backbones,major,cols,rows,time_limit=10.,candidate_limit=10,progress=print,*,validate_candidate=None):
    """Evaluate routed candidates before orange/gray nodes become obstacles."""
    from adaptive_routing import make_paths,route_paths,validate_paths,validate_node_grid
    started=time.monotonic();deadline=started+time_limit
    base_paths=[p for p in make_paths(nodes,edges,xy,backbones,major) if p['priority']==2]
    keys={k for p in base_paths for k in (p['source'],p['target'])}
    ring=set((backbones.get('tca_exclusion',{}).get('region') or {}).get('node_ids',[]))
    ring &= keys
    bounds=[(0,0),(0,0)]
    if ring:
        bounds=[(math.ceil(1-min(xy[k][a] for k in ring)),
                 math.floor((cols-2 if a==0 else rows-2)-max(xy[k][a] for k in ring))) for a in (0,1)]
    if any(lo>hi for lo,hi in bounds):bounds=[(0,0),(0,0)]
    seed,solver=translation_seed(base_paths,xy,ring,bounds,time_limit)
    # The MILP supplies a lower-bound seed; the actual route lengths decide.
    translations={(0,0),seed}
    area=math.prod(hi-lo+1 for lo,hi in bounds)
    if area<=20000:
        translations.update((x,y) for x in range(bounds[0][0],bounds[0][1]+1)
                            for y in range(bounds[1][0],bounds[1][1]+1))
    else:
        translations.update((seed[0]+x,seed[1]+y) for x in range(-16,17) for y in range(-16,17)
                            if all(bounds[a][0]<=v<=bounds[a][1] for a,v in enumerate((seed[0]+x,seed[1]+y))))
    def proxy(shift):
        return sum(sum(abs(xy[p['source']][a]-xy[p['target']][a]+
                           shift[a]*(int(p['source'] in ring)-int(p['target'] in ring))) for a in (0,1))
                   for p in base_paths if not p.get('curve_controls'))
    ranked=sorted(((s,g) for s in translations for g in (0,2,4,6)),
                  key=lambda item:(proxy(item[0]),abs(item[0][0]-seed[0])+abs(item[0][1]-seed[1]),item[1],item[0]))
    candidates=[];failures=[];evaluated=0;baseline=None
    # Keep the original geometry as an actual routed baseline when feasible.
    ranked=[((0,0),0),*[item for item in ranked if item!=((0,0),0)]]
    for shift,global_dx in ranked:
        if evaluated>=candidate_limit or (time.monotonic()>deadline and candidates):break
        positions,trial=translated_core(xy,backbones,ring,shift,rows,global_dx=global_dx,red_keys=keys)
        region=trial.get('tca_exclusion',{}).get('region')
        red_xy={k:positions[k] for k in keys}
        if len({tuple(p) for p in red_xy.values()})!=len(red_xy):continue
        if any(not (0<=p[0]<=cols-1 and 0<=p[1]<=rows-1) for p in red_xy.values()):continue
        if any(not allowed_points([p],region)[0] for k,p in red_xy.items() if k not in ring):continue
        access=None
        if validate_candidate:
            try:access=validate_candidate(positions,trial,None)
            except (ValueError,RuntimeError):continue
        evaluated+=1
        progress(f'Red core candidate {evaluated}: TCA translation {shift}, core x translation {global_dx}',flush=True)
        for order in (0,1):
            if order and time.monotonic()>deadline and candidates:break
            try:
                red=[p for p in make_paths(nodes,edges,positions,trial,major) if p['priority']==2]
                routing_xy=dict(red_xy,**(access['positions'] if access else {}))
                boundary=access['boundary'] if access else [0,cols-1,0,rows-1]
                route_paths(red,routing_xy,boundary,region,red_order=order,reservations=access['paths'] if access else ())
                validate_paths(red,red_xy,boundary,region,[])
                validate_node_grid({k:nodes[k] for k in keys},red_xy,region)
                length=red_length(red)
                if shift==(0,0) and global_dx==0 and order==0:baseline=length
                if validate_candidate:validate_candidate(positions,trial,red)
                signature=reference_signature(red_xy,red)
                if any(c['signature']==signature for c in candidates):continue
                candidates.append(dict(positions=red_xy,routes=copy.deepcopy(red),backbones=trial,
                                       length=length,shift=[shift[0]+global_dx,shift[1]],relative_tca_shift=list(shift),
                                       global_translation=[global_dx,0],order=order,signature=signature))
            except (ValueError,RuntimeError) as error:
                failures.append(dict(shift=list(shift),order=order,reason=str(error)[:240]))
    if not candidates:raise ValueError('No valid length-first red core: '+str(failures[-4:]))
    candidates.sort(key=lambda c:(round(c['length'],7),sum(abs(v) for v in c['shift']),c['order']))
    reference=candidates[0]
    report=dict(objective='Sum of drawn red incidence lengths, including cubic TCA arcs',
                objective_order=['hard constraints','red total length','red reference changes','orange length','margins'],
                solver=solver,initial_routed_length_grid=baseline,optimized_length_grid=reference['length'],
                evaluated_placements=evaluated,feasible_candidates=len(candidates),rejected_candidates=failures,
                time_limit_seconds=time_limit,elapsed_seconds=time.monotonic()-started,
                global_optimality_proven=False,reference_positions=copy.deepcopy(reference['positions']),
                reference_routes=copy.deepcopy(reference['routes']),reference_signature=reference['signature'])
    progress(f'Red core: {baseline} -> {reference["length"]:.3f} grid units; {len(candidates)} valid candidates',flush=True)
    return candidates,report
