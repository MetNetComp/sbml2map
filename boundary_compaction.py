"""Shrink lower/right EX boundaries, preserving red geometry and orange length."""
import copy
import math
from crossing_routing import plan_routes,validate_boundary,orange_length,EPS
from exchange_corridors import place_corridors
from exchange_alignment import validate_alignment
from tca_geometry import validate_exclusion


def compact_boundary(model,nodes,edges,xy,backbones,exchange_nodes,groups,cols,rows,
                     major,attach,attempts=12,margin=1.,max_candidates=24,progress=print):
    original_boundary=list(backbones.get('routing_boundary',[0,cols-1,0,rows-1]))
    backbones['routing_boundary']=original_boundary
    baseline=plan_routes(nodes,edges,xy,backbones,major,attempts=attempts,progress=progress)
    backbones['crossing_routing']=baseline
    before_length=orange_length(baseline['routes'])
    red=backbones.get('red_first',{}).get('routes',[])
    report=dict(enabled=True,before_boundary=original_boundary,after_boundary=original_boundary[:],
                before_orange_length_grid=before_length,after_orange_length_grid=before_length,
                margin_grid=margin,accepted=False,red_coordinates_preserved=True,red_routes_preserved=True,
                attempted_candidates=0,rejected_candidates=[],global_optimality_proven=False)
    if not red or not groups:
        report['reason']='No frozen red geometry or matched exchange triplets'
        backbones['boundary_compaction']=report;return report
    red_keys={k for p in red for k in (p['source'],p['target'])}
    protected=validate_boundary(red,original_boundary)['route_extents']
    target_right=min(original_boundary[1],max(original_boundary[0]+3,math.ceil(protected[1]+margin)))
    target_bottom=max(original_boundary[2],min(original_boundary[3]-3,math.floor(protected[2]-margin)))
    report['protected_red_extents']=protected
    candidates=[]
    for right in range(int(target_right),int(original_boundary[1])+1):
        for bottom in range(int(target_bottom),int(original_boundary[2])-1,-1):
            b=[original_boundary[0],right,bottom,original_boundary[3]]
            if b==original_boundary:continue
            candidates.append(b)
    candidates.sort(key=lambda b:((b[1]-b[0])*(b[3]-b[2]),b[1],-b[2]))
    region=backbones.get('tca_exclusion',{}).get('region')
    original_xy=copy.deepcopy(xy);original_backbones=copy.deepcopy(backbones)
    for boundary in candidates[:max_candidates]:
        report['attempted_candidates']+=1
        trial_xy=copy.deepcopy(original_xy);trial=copy.deepcopy(original_backbones)
        trial.pop('crossing_routing',None);trial['routing_boundary']=boundary
        try:
            trial['exchange_corridors']=place_corridors(model,nodes,trial_xy,trial,exchange_nodes,groups,cols,rows)
            attach(nodes,trial_xy,grid_shape=(cols,rows),forbidden_region=region,boundary=boundary)
            if any(math.dist(trial_xy[k],original_xy[k])>1e-9 for k in red_keys):raise ValueError('Frozen red node moved')
            if any(not (boundary[0]-EPS<=p[0]<=boundary[1]+EPS and boundary[2]-EPS<=p[1]<=boundary[3]+EPS) for p in trial_xy.values()):raise ValueError('Node outside candidate boundary')
            if any(not (abs(trial_xy[k][0]-boundary[0])<EPS or abs(trial_xy[k][0]-boundary[1])<EPS or abs(trial_xy[k][1]-boundary[2])<EPS or abs(trial_xy[k][1]-boundary[3])<EPS) for k in exchange_nodes):raise ValueError('Exchange node left the perimeter')
            ordinary=[k for k in trial_xy if nodes[k]['kind']!='currency']
            if len({tuple(trial_xy[k]) for k in ordinary})!=len(ordinary):raise ValueError('Ordinary node positions overlap')
            validate_exclusion(nodes,trial_xy,region)
            validate_alignment(groups,trial_xy,cols,rows,boundary)
            plan=plan_routes(nodes,edges,trial_xy,trial,major,attempts=attempts,progress=lambda *a,**k:None)
            length=orange_length(plan['routes'])
            if length>before_length+1e-7:raise ValueError('Orange length would increase: '+str(length))
            if not plan['report'].get('red_routes_preserved'):raise ValueError('Frozen red routing changed')
            trial['crossing_routing']=plan;trial['render_frame']=boundary
        except (ValueError,RuntimeError) as error:
            report['rejected_candidates'].append(dict(boundary=boundary,reason=str(error)[:400]));continue
        report.update(accepted=True,after_boundary=boundary[:],after_orange_length_grid=length,
                      lower_margin_removed_grid=boundary[2]-original_boundary[2],right_margin_removed_grid=original_boundary[1]-boundary[1],
                      before_area_grid=(original_boundary[1]-original_boundary[0])*(original_boundary[3]-original_boundary[2]),
                      after_area_grid=(boundary[1]-boundary[0])*(boundary[3]-boundary[2]))
        xy.clear();xy.update(trial_xy);backbones.clear();backbones.update(trial)
        progress(f'Boundary compacted: {original_boundary} -> {boundary}; orange length {before_length:.3f} -> {length:.3f}',flush=True)
        break
    backbones['boundary_compaction']=report
    return report
