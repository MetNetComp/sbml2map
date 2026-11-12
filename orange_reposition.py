"""Bounded local relocation to remove positive-length orange overlaps.

Point crossings are not penalized. Same-reaction shared trunks are evaluated.
Every accepted candidate is fully routed and checked against the frozen red
geometry, exchange corridors, TCA exclusion and active exchange boundary.
"""
import copy
import math
from collections import defaultdict
from crossing_routing import EPS, orange_length, plan_routes
from exchange_alignment import validate_alignment
from exchange_corridors import validate_corridors
from tca_geometry import validate_exclusion


def orange_overlaps(paths):
    """Positive-length overlaps, including different incidences of one reaction.

    Distinct reaction pairs retain union counting. Within a reaction each pair
    of biological incidences is counted separately. Point crossings cost zero.
    """
    lines = defaultdict(lambda: defaultdict(list))
    for p in paths:
        if p['priority'] != 1:
            continue
        incidence = (p['rid'], p['mid'])
        for a, b in zip(p['points'], p['points'][1:]):
            if math.dist(a, b) < EPS:
                continue
            axis = 0 if abs(a[1]-b[1]) < EPS else 1
            lines[axis, round(a[1-axis], 8)][incidence].append(sorted((a[axis], b[axis])))
    overlaps = []
    for (axis, level), incidences in sorted(lines.items()):
        ids = sorted(incidences)
        intervals = defaultdict(list)
        for i, left in enumerate(ids):
            for right in ids[i+1:]:
                group = (left, right) if left[0] == right[0] else ((left[0], ''), (right[0], ''))
                for a, b in incidences[left]:
                    for c, d in incidences[right]:
                        low, high = max(a, c), min(b, d)
                        if high-low > EPS:
                            intervals[group].append((low, high))
        for pair, values in sorted(intervals.items()):
            merged = []
            for low, high in sorted(values):
                if merged and low <= merged[-1][1]+EPS:
                    merged[-1][1] = max(merged[-1][1], high)
                else:
                    merged.append([low, high])
            for low, high in merged:
                a, b = [0., 0.], [0., 0.]
                a[axis], b[axis] = low, high
                a[1-axis] = b[1-axis] = level
                overlaps.append(dict(reactions=[p[0] for p in pair], incidences=[list(p) for p in pair],
                                     same_reaction=pair[0][0] == pair[1][0], points=[a, b], length_grid=high-low))
    return dict(length_grid=sum(p['length_grid'] for p in overlaps), segments=overlaps)


def relocation_candidates(nodes, xy, backbones, paths, radius=3, overlap=None):
    """Move affected ordinary endpoints; pair corridors retain their boundary axis."""
    overlaps = orange_overlaps(paths)
    affected = set(overlap['reactions']) if overlap else {r for s in overlaps['segments'] for r in s['reactions']}
    self_reactions={s['reactions'][0] for s in ([overlap] if overlap else overlaps['segments']) if s.get('same_reaction')}
    fixed = set(backbones.get('fixed_node_offsets', {}))
    fixed.update(k for p in backbones.get('red_first', {}).get('routes', []) for k in (p['source'], p['target']))
    groups = backbones.get('exchange_corridors', {}).get('groups', [])
    pairs = {g['external']: g for g in groups if g['kind'] == 'pair'}
    locked = {k for g in groups for k in g['node_ids']} - set(pairs)
    keys = {k for p in paths if p['priority'] == 1 and p['rid'] in affected for k in (p['source'], p['target'])}
    boundary = backbones['routing_boundary']
    region = backbones.get('tca_exclusion', {}).get('region')
    out = []
    for key in sorted(keys, key=lambda k: (not (k.startswith('rxn:') and k[4:] in self_reactions), k not in pairs, k)):
        if key in fixed | locked or nodes[key]['kind'] == 'currency':
            continue
        base = xy[key]
        group = pairs.get(key)
        axis = None
        if group:
            a, b = group['corridor']
            axis = 0 if abs(a[1]-b[1]) < EPS else 1
        offsets = [(dx, dy) for dx in range(-radius, radius+1) for dy in range(-radius, radius+1)
                   if 0 < abs(dx)+abs(dy) <= radius and (axis is None or (dy == 0 if axis == 0 else dx == 0))]
        # Upward candidates precede downward ones in ties (Cartesian y grows up).
        offsets.sort(key=lambda d: (abs(d[0])+abs(d[1]), d[0], -d[1]))
        for dx, dy in offsets:
            point = [base[0]+dx, base[1]+dy]
            if not (boundary[0]+.16 < point[0] < boundary[1]-.16 and boundary[2]+.16 < point[1] < boundary[3]-.16):
                continue
            if any(math.dist(point, p) < .16-EPS for k, p in xy.items() if k != key):
                continue
            if group:
                a = xy[group['exchange']]
                # Keep the same face and an inward-only straight pair.
                direction = base[axis]-a[axis]
                if (point[axis]-a[axis])*direction <= 0:
                    continue
            test = dict(xy);test[key] = point
            try:
                validate_exclusion(nodes, test, region)
            except (ValueError, RuntimeError):
                continue
            out.append((key, point))
    return out


def _pair_length(paths, overlap):
    pair=overlap['incidences']
    return sum(s['length_grid'] for s in orange_overlaps(paths)['segments'] if s['incidences']==pair)


def _trial_relocation(nodes, edges, xy, backbones, key, point, major, groups, cols, rows, attempts):
    trial_xy=copy.deepcopy(xy);trial=copy.deepcopy(backbones)
    old=trial_xy[key];trial_xy[key]=list(point)
    for child,node in nodes.items():
        if node['kind']=='currency' and node['parent_id']==key:
            trial_xy[child]=[xy[child][i]+point[i]-old[i] for i in (0,1)]
    for g in trial.get('exchange_corridors',{}).get('groups',[]):
        if key==g['external'] and g['kind']=='pair':
            g['corridor']=[list(trial_xy[g['exchange']]),list(point)]
    boundary=trial['routing_boundary']
    if any(not (boundary[0]<=p[0]<=boundary[1] and boundary[2]<=p[1]<=boundary[3]) for p in trial_xy.values()):
        raise ValueError('Node outside exchange boundary')
    if any(math.dist(trial_xy[a],trial_xy[b])<.16-EPS for a in trial_xy for b in trial_xy if a<b and (a==key or b==key)):
        raise ValueError('Node positions overlap')
    validate_exclusion(nodes,trial_xy,trial.get('tca_exclusion',{}).get('region'))
    validate_alignment(groups,trial_xy,cols,rows,boundary)
    candidate=plan_routes(nodes,edges,trial_xy,trial,major,attempts=attempts,progress=lambda *a,**k:None)
    validate_corridors(candidate['routes'],trial.get('exchange_corridors',{'groups':[]}))
    return trial_xy,trial,candidate


def reposition_orange(nodes, edges, xy, backbones, major, groups, cols, rows,
                      attempts=4, max_candidates=4, progress=print, max_passes=2, radii=(1,3,6)):
    """Give every overlap its own budget; retry with wider neighborhoods.

    max_candidates is a budget PER overlap and radius, never a global cutoff.
    Accepted candidates must reduce both their targeted pair and total overlap.
    """
    if max_candidates<1 or max_passes<1 or not radii or any(r<1 for r in radii):
        raise ValueError('Positive per-overlap budget, passes and radii are required')
    backbones.setdefault('routing_boundary',[0,cols-1,0,rows-1])
    plan=backbones.get('crossing_routing')
    if plan is None:plan=plan_routes(nodes,edges,xy,backbones,major,attempts=attempts,progress=progress)
    initial=orange_overlaps(plan['routes'])
    report=dict(enabled=True,before_overlap_grid=initial['length_grid'],after_overlap_grid=initial['length_grid'],
                before_orange_length_grid=orange_length(plan['routes']),after_orange_length_grid=orange_length(plan['routes']),
                attempted_candidates=0,accepted_moves=[],rejected_candidates=[],overlap_searches=[],
                objective_order=['all_orange_overlap_length','orange_length','displacement'],
                point_crossings_penalized=False,shared_same_reaction_trunks_penalized=True,
                global_optimality_proven=False,candidate_limit_per_overlap_radius=max_candidates,
                search_radii=list(radii),max_passes=max_passes,routing_attempts_per_candidate=min(attempts,1))
    tried=set();completed_passes=0
    for sweep in range(max_passes):
        improved=False;completed_passes=sweep+1
        # A reaction-pair group is revisited after accepted moves in later sweeps.
        segments=orange_overlaps(plan['routes'])['segments']
        unique={tuple(tuple(p) for p in s['incidences']):s for s in segments}
        ordered=sorted(unique.values(),key=lambda s:(not s['same_reaction'],s['length_grid'],s['incidences']))
        for overlap in ordered:
            current=orange_overlaps(plan['routes'])['length_grid'];pair_before=_pair_length(plan['routes'],overlap)
            if pair_before<=EPS:continue
            search=dict(incidences=overlap['incidences'],pass_index=sweep+1,attempted_candidates=0,
                        before_pair_overlap_grid=pair_before,accepted=False,radii_tried=[],feasible_candidates=0,infeasible_candidates=0)
            best=None
            for radius in radii:
                search['radii_tried'].append(radius);tested=0
                candidates=relocation_candidates(nodes,xy,backbones,plan['routes'],radius,overlap)
                # Round-robin the affected nodes, so one immovable region does
                # not consume another endpoint's entire candidate budget.
                by_node=defaultdict(list)
                for key,point in candidates:by_node[key].append(point)
                flattened=[]
                for index in range(max(map(len,by_node.values()),default=0)):
                    for key,points in by_node.items():
                        if index<len(points):flattened.append((key,points[index]))
                for key,point in flattened:
                    signature=(plan['signature'],key,tuple(point))
                    if signature in tried:continue
                    tried.add(signature);tested+=1;search['attempted_candidates']+=1;report['attempted_candidates']+=1
                    try:
                        trial_xy,trial,candidate=_trial_relocation(nodes,edges,xy,backbones,key,point,major,groups,cols,rows,min(attempts,1))
                        search['feasible_candidates']+=1
                        total=orange_overlaps(candidate['routes'])['length_grid'];pair_after=_pair_length(candidate['routes'],overlap)
                        score=(round(total,7),round(orange_length(candidate['routes']),7),math.dist(xy[key],point))
                        if total<current-EPS and pair_after<pair_before-EPS and (best is None or score<best[0]):
                            best=(score,key,list(xy[key]),list(point),trial_xy,trial,candidate,pair_after)
                        if best and pair_after<=EPS and orange_length(candidate['routes'])<=orange_length(plan['routes'])+EPS:
                            break
                    except (ValueError,RuntimeError) as error:
                        search['infeasible_candidates']+=1
                        report['rejected_candidates'].append(dict(node=key,position=point,reason=str(error)[:250]))
                    if tested>=max_candidates:break
                if best is not None:break
            if best:
                score,key,before,after,trial_xy,trial,plan,pair_after=best
                xy.clear();xy.update(trial_xy);backbones.clear();backbones.update(trial)
                report['accepted_moves'].append(dict(node=key,before=before,after=after,overlap_grid=score[0],orange_length_grid=score[1],
                                                     target_incidences=overlap['incidences']))
                search.update(accepted=True,after_pair_overlap_grid=pair_after);improved=True
                progress(f'Orange overlap relocation: {key} {before} -> {after}; overlap {current:.3f} -> {score[0]:.3f}',flush=True)
            else:
                search['reason']=('No movable endpoint or allowed candidate' if not search['attempted_candidates'] else
                                  'All tested candidates violate hard constraints; wider bounded search may still find a solution' if not search['feasible_candidates'] else
                                  'No improving candidate in bounded search; untested positions may improve the layout')
            report['overlap_searches'].append(search)
        if not improved:break
    backbones['crossing_routing']=plan
    remaining=orange_overlaps(plan['routes'])
    report.update(after_overlap_grid=remaining['length_grid'],after_orange_length_grid=orange_length(plan['routes']),
                  remaining_overlaps=remaining['segments'],red_routes_preserved=plan['report'].get('red_routes_preserved',False),
                  completed_passes=completed_passes,termination='pass_limit' if improved and completed_passes==max_passes else 'no_tested_improvement')
    unresolved=[]
    for segment in remaining['segments']:
        searches=[s for s in report['overlap_searches'] if s['incidences']==segment['incidences']]
        unresolved.append(dict(incidences=segment['incidences'],length_grid=segment['length_grid'],
                               reason=searches[-1].get('reason','Overlap remains after bounded passes') if searches else 'Created after the last search; pass limit reached'))
    report['unresolved_overlaps']=unresolved
    backbones['orange_reposition']=report
    return report


def place_gray_inside_tca(nodes, edges, xy, backbones):
    """Use newly allowed interior cells for unconstrained gray nodes when shorter.

    Orange and red endpoints keep a feasible initial embedding; orange moves
    into the interior are subsequently checked by the full route planner.
    """
    region=backbones.get('tca_exclusion',{}).get('region')
    if not region or region.get('interior_policy')!='allow':return dict(enabled=False,moves=[])
    from tca_geometry import allowed_points
    from crossing_routing import initial_paths
    major={rid for route in backbones.get('routes',[]) for rid in route['reaction_ids']}
    paths=initial_paths(nodes,edges,xy,backbones,major)
    fixed={k for p in paths if p['priority'] for k in (p['source'],p['target'])}
    fixed.update(backbones.get('fixed_node_offsets',{}))
    fixed.update(backbones.get('exchange_corridors',{}).get('straight_node_ids',[]))
    fixed.update(k for k in xy if k.startswith('rxn:EX_'))
    boundary=backbones['routing_boundary'];outside=dict(region,interior_policy='forbid')
    cells=[[float(x),float(y)] for x in range(math.floor(boundary[0])+1,math.ceil(boundary[1]))
           for y in range(math.floor(boundary[2])+1,math.ceil(boundary[3]))]
    cells=[q for q in cells if allowed_points([q],region)[0] and not allowed_points([q],outside,0)[0]]
    neighbors=defaultdict(list)
    for e in edges:neighbors[e['source']].append(e['target']);neighbors[e['target']].append(e['source'])
    moves=[]
    for key in sorted(xy):
        if key in fixed or nodes[key]['kind']=='currency':continue
        def cost(point):return sum(abs(point[i]-xy[k][i]) for k in neighbors[key] if nodes[k]['kind']!='currency' for i in (0,1))
        before=list(xy[key]);baseline=cost(before)
        for point in sorted(cells,key=lambda q:(cost(q),math.dist(q,before),q)):
            if cost(point)>=baseline-EPS:break
            trial=dict(xy);trial[key]=point
            for child,node in nodes.items():
                if node['kind']=='currency' and node['parent_id']==key:
                    trial[child]=[xy[child][i]+point[i]-before[i] for i in (0,1)]
            changed={k for k in trial if trial[k]!=xy[k]}
            if any(math.dist(trial[a],trial[b])<.16-EPS for a in changed for b in trial if a!=b):continue
            if any(not (boundary[0]<=trial[k][0]<=boundary[1] and boundary[2]<=trial[k][1]<=boundary[3]) for k in changed):continue
            try:validate_exclusion(nodes,trial,region)
            except (ValueError,RuntimeError):continue
            xy.update(trial);moves.append(dict(node=key,before=before,after=list(point),before_length_grid=baseline,after_length_grid=cost(point)));break
    # Gray coordinates do not contribute to compressed orange routing axes.
    # Update only gray paths; preserve every saved red/orange route exactly.
    if backbones.get('crossing_routing'):
        plan=backbones['crossing_routing']
        fresh={(p['rid'],p['mid']):p for p in initial_paths(nodes,edges,xy,backbones,major)}
        for index,p in enumerate(plan['routes']):
            if not p['priority']:plan['routes'][index]=fresh[p['rid'],p['mid']]
        import hashlib,json
        plan['signature']=hashlib.sha256(json.dumps(plan['routes'],sort_keys=True).encode()).hexdigest()[:24]
    return dict(enabled=True,moves=moves,orange_coordinates_preserved=True,red_coordinates_preserved=True)
