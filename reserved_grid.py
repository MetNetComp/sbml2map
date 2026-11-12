"""Lattice placement with fixed circular anchors, forbidden interiors and EX triplets.

Both the incumbent heuristic and finite-cell MILP use exactly the same allowed
positions. Fixed TCA anchors and aligned partners may use fractional coordinates.
"""
from collections import defaultdict
import math
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree
from tca_geometry import absolute_region, allowed_points, EPS
from exchange_alignment import triplet_face, validate_alignment


def optimize_reserved_grid(node_ids, edge_pairs, cols=None, rows=None, time_limit=60, seed=0,
        progress=print, top_left_node=None, boundary_nodes=(), fixed_offsets=None,
        top_left_margin=0, node_top_margins=None, max_milp_nodes=None,
        exclusive_boundary=False, initial_positions=None, edge_weights=None,
        reserved_region=None, alignment_groups=()):
    from grid_layout2 import validated_edge_weights
    ids = list(node_ids)
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate grid node IDs')
    fixed_offsets = dict(fixed_offsets or {})
    margins = dict(node_top_margins or {})
    boundary = set(boundary_nodes)
    if (set(fixed_offsets)|set(margins)|boundary)-set(ids):
        raise ValueError('Unknown constrained grid node')
    for point in fixed_offsets.values():
        if len(point) != 2 or not all(math.isfinite(v) for v in point) or point[0] < 0 or point[1] > 0:
            raise ValueError('Backbone offsets must be finite [x>=0,y<=0] values')
    if top_left_margin < 0 or int(top_left_margin) != top_left_margin or any(v < 0 or int(v) != v for v in margins.values()):
        raise ValueError('Top margins must be nonnegative integers')
    for group in alignment_groups:
        if any(group[k] not in ids for k in ('exchange','external','internal')):
            raise ValueError('Unknown aligned node')
    if top_left_node is not None and top_left_node not in ids:
        raise ValueError('Unknown top-left node')
    auto_cols, auto_rows = cols is None, rows is None
    if cols is not None and (cols < 1 or int(cols) != cols) or rows is not None and (rows < 1 or int(rows) != rows):
        raise ValueError('Grid dimensions must be positive integers')
    n = len(ids)
    cols = cols if cols is not None else max(3, math.ceil(math.sqrt(max(1,n)*1.3)))
    rows = rows if rows is not None else max(3, math.ceil(max(1,n)*1.3/cols))
    if fixed_offsets:
        need_cols = math.ceil(max(x for x,y in fixed_offsets.values()))+2
        need_rows = math.ceil(-min(y for x,y in fixed_offsets.values()))+2
        if auto_cols:cols=max(cols,need_cols)
        elif cols < need_cols-1:raise ValueError('Grid has too few columns for TCA/backbone')
        if auto_rows:rows=max(rows,need_rows)
        elif rows < need_rows-1:raise ValueError('Grid has too few rows for TCA/backbone')
    index = {key:i for i,key in enumerate(ids)}
    edge_indices = np.asarray([(index[a],index[b]) for a,b in edge_pairs], dtype=int).reshape((-1,2))
    weights = validated_edge_weights(edge_weights,len(edge_pairs))
    neighbors = defaultdict(list)
    for (a,b),w in zip(edge_pairs,weights):neighbors[a].append((b,w));neighbors[b].append((a,w))
    def value(xy):
        return float(sum(w*(abs(xy[a][0]-xy[b][0])+abs(xy[a][1]-xy[b][1])) for (a,b),w in zip(edge_pairs,weights)))
    aligned = {g[k] for g in alignment_groups for k in ('exchange','external','internal')}

    def prepare():
        fixed = {key:(float(x),float(rows-1+y)) for key,(x,y) in fixed_offsets.items()}
        if top_left_node is not None:
            anchor = (float(int(exclusive_boundary)),float(rows-1-top_left_margin))
            if top_left_node in fixed and math.dist(anchor,fixed[top_left_node]) > EPS:
                raise ValueError('Backbone position conflicts with the top-left anchor')
            fixed[top_left_node] = anchor
        region = absolute_region(reserved_region,rows)
        ring = set((region or {}).get('node_ids',()))
        base = [(float(x),float(y)) for y in range(rows) for x in range(cols)]
        # Matching a fractional TCA anchor requires fractional rows/columns for
        # that exchange and extracellular node. Ordinary nodes stay on lattice.
        extra = set()
        for group in alignment_groups:
            for key in (group['internal'],group['external'],group['exchange']):
                if key in fixed:
                    x,y=fixed[key]
                    extra.update((float(col),y) for col in range(cols))
                    extra.update((x,float(row)) for row in range(rows))
        choices={}
        for key in ids:
            candidates=[fixed[key]] if key in fixed else base+(sorted(extra) if key in aligned else [])
            unique=sorted(set(candidates))
            positions=np.asarray(unique,dtype=float).reshape((-1,2))
            valid=np.ones(len(positions),dtype=bool)
            if len(positions):
                valid &= (positions[:,0]>=-EPS)&(positions[:,0]<=cols-1+EPS)&(positions[:,1]>=-EPS)&(positions[:,1]<=rows-1+EPS)
                on_boundary=(np.abs(positions[:,0])<EPS)|(np.abs(positions[:,0]-cols+1)<EPS)|(np.abs(positions[:,1])<EPS)|(np.abs(positions[:,1]-rows+1)<EPS)
                if key in boundary:valid &= on_boundary
                elif exclusive_boundary:valid &= ~on_boundary
                valid &= positions[:,1] <= rows-1-margins.get(key,0)+EPS
                if key not in ring:valid &= allowed_points(positions,region)
                for other,p in fixed.items():
                    if other!=key:valid &= np.linalg.norm(positions-np.asarray(p),axis=1)>=.16-EPS
            choices[key]=[tuple(v) for v in positions[valid]]
        return fixed,region,choices

    def group_candidates(group,xy,choices):
        ex,external,internal=(group[k] for k in ('exchange','external','internal'))
        valid_sets={key:set(choices[key]) for key in (ex,external,internal)}
        anchors=[xy[internal]] if internal in xy else choices[internal]
        for c in anchors:
            for face in ('left','right','bottom','top'):
                if face in ('left','right'):
                    a=(0.,c[1]) if face=='left' else (float(cols-1),c[1])
                    options=[(float(x),c[1]) for x in range(1,cols-1) if min(a[0],c[0])+.16-EPS<=x<=max(a[0],c[0])-.16+EPS]
                else:
                    a=(c[0],0.) if face=='bottom' else (c[0],float(rows-1))
                    options=[(c[0],float(y)) for y in range(1,rows-1) if min(a[1],c[1])+.16-EPS<=y<=max(a[1],c[1])-.16+EPS]
                if a not in valid_sets[ex] or ex in xy and math.dist(xy[ex],a)>EPS:continue
                for b in options:
                    if b not in valid_sets[external] or external in xy and math.dist(xy[external],b)>EPS:continue
                    result={ex:a,external:b,internal:c}
                    if not triplet_face(group,result,cols,rows):continue
                    yield result

    def initial_placement(fixed,choices):
        if any(not choices[k] for k in ids):raise ValueError('No allowed cell for a node')
        xy=dict(fixed)
        if initial_positions is not None:
            proposed={k:tuple(initial_positions[k]) for k in ids}
            if all(proposed[k] in choices[k] for k in ids) and all(triplet_face(g,proposed,cols,rows) for g in alignment_groups):
                if min((math.dist(proposed[a],proposed[b]) for i,a in enumerate(ids) for b in ids[i+1:]),default=1)>=.16-EPS:
                    return proposed
        for group in sorted(alignment_groups,key=lambda g:(-sum(g[k] in xy for k in ('exchange','external','internal')),g['reaction_id'])):
            best=None
            for candidate in group_candidates(group,xy,choices):
                if any(math.dist(p,q)<.16-EPS for k,p in candidate.items() for other,q in xy.items() if k!=other and other not in candidate):continue
                trial=dict(xy,**candidate)
                score=sum(w*math.dist(trial[a],trial[b]) for a in candidate for b,w in neighbors[a] if b in trial)
                score+=.001*sum(abs(p[0]-(cols-1)/2)+abs(p[1]-(rows-1)/2) for p in candidate.values())
                marker=(score,tuple(sorted(candidate.items())))
                if best is None or marker<best[0]:best=(marker,candidate)
            if best is None:raise ValueError('No feasible EX/extracellular/cytosolic straight placement: '+group['reaction_id'])
            xy.update(best[1])
        pending=set(ids)-set(xy)
        while pending:
            key=min(pending,key=lambda k:(-sum(b in xy for b,w in neighbors[k]),k))
            occupied=np.asarray(list(xy.values()),dtype=float).reshape((-1,2))
            candidates=np.asarray(choices[key])
            valid=np.ones(len(candidates),dtype=bool)
            for p in occupied:valid &= np.linalg.norm(candidates-p,axis=1)>=.16-EPS
            candidates=candidates[valid]
            if not len(candidates):raise ValueError('Not enough permitted cells outside TCA')
            known=[(xy[b],w) for b,w in neighbors[key] if b in xy]
            if known:
                costs=sum(w*np.abs(candidates-np.asarray(p)).sum(axis=1) for p,w in known)
            else:costs=np.abs(candidates-np.asarray([(cols-1)/2,(rows-1)/2])).sum(axis=1)
            xy[key]=tuple(candidates[int(np.argmin(costs))]);pending.remove(key)
        return xy

    for attempt in range(80):
        fixed,region,choices=prepare()
        try:xy=initial_placement(fixed,choices);break
        except ValueError as exc:
            if not auto_cols and not auto_rows:raise ValueError('Explicit grid cannot satisfy TCA/exchange constraints: '+str(exc)) from exc
            if auto_cols and (not auto_rows or cols<=rows):cols+=1
            else:rows+=1
    else:raise ValueError('Could not fit TCA/exchange constraints after automatic grid expansion')
    progress(f'Grid: {n} nodes on {cols} x {rows}; TCA interior {(region or {}).get('interior_policy','forbid')}; {len(alignment_groups)} straight exchange triplets',flush=True)
    initial=value(xy)
    pinned=set(fixed)|aligned
    movable=[k for k in ids if k not in pinned]
    # Local relocations/swaps optimize the original weighted L1 objective.
    for _ in range(3):
        improved=False
        owner={tuple(p):key for key,p in xy.items()}
        for key in movable:
            current=xy[key]
            targets=choices[key]
            if len(targets)>120:
                center=np.median([xy[b] for b,w in neighbors[key]],axis=0) if neighbors[key] else current
                targets=sorted(targets,key=lambda p:(abs(p[0]-center[0])+abs(p[1]-center[1]),p))[:120]
            best=None
            current_keys=list(xy)
            near=cKDTree(np.asarray(list(xy.values()))).query_ball_point(np.asarray(targets),.16-EPS) if targets else []
            for target,close in zip(targets,near):
                partner=owner.get(target)
                if partner==key or partner is not None and partner not in movable:continue
                if partner and current not in choices[partner]:continue
                changed={key:target}
                if partner:changed[partner]=current
                if any(current_keys[i] not in changed for i in close):continue
                affected=[(a,b,w) for (a,b),w in zip(edge_pairs,weights) if a in changed or b in changed]
                delta=sum(w*((abs(changed.get(a,xy[a])[0]-changed.get(b,xy[b])[0])+abs(changed.get(a,xy[a])[1]-changed.get(b,xy[b])[1]))-(abs(xy[a][0]-xy[b][0])+abs(xy[a][1]-xy[b][1]))) for a,b,w in affected)
                if delta < -EPS and (best is None or delta<best[0]):best=(delta,changed)
            if best:
                for k in best[1]:owner.pop(xy[k],None)
                xy.update(best[1]);owner.update({p:k for k,p in best[1].items()});improved=True
        if not improved:break
    incumbent=value(xy)
    report=dict(grid_cols=cols,grid_rows=rows,objective=incumbent,objective_metric='Weighted Manhattan incidence cost in grid units',
                optimality_proven=False,global_lower_bound=0.,relative_gap=1. if incumbent else 0.,
                method='Constraint-preserving lattice heuristic with fixed TCA anchors',selected_solution='validated heuristic incumbent',
                heuristic=dict(initial_length=initial,final_length=incumbent),
                position_constraints=dict(top_left_node=top_left_node,top_left_position=list(fixed[top_left_node]) if top_left_node else None,
                    top_left_verified=True,boundary_nodes_verified=True,boundary_node_ids=sorted(boundary),exclusive_boundary=exclusive_boundary,
                    fixed_backbone_positions={k:list(v) for k,v in fixed.items()},fixed_backbone_positions_verified=True,
                    tca_interior_excluded=bool(region) and region.get('interior_policy')!='allow',exchange_alignment=validate_alignment(alignment_groups,xy,cols,rows)))
    # Each allowed position is an explicit binary choice. The same masks enforce
    # the exclusion region in both the heuristic and MILP, including on timeout.
    variables=[(k,point) for k in ids if k not in fixed for point in choices[k]]
    if (max_milp_nodes is None or n<=max_milp_nodes) and len(variables)<=100000 and len(edge_pairs) and time_limit>0:
        vcount=len(variables);e=len(edge_pairs);groups=list(alignment_groups)
        size=2*n+vcount+2*e+4*len(groups)
        cost=np.zeros(size);cost[2*n+vcount:2*n+vcount+e]=weights;cost[2*n+vcount+e:2*n+vcount+2*e]=weights
        low=np.zeros(size);high=np.ones(size)
        high[:n]=cols-1;high[n:2*n]=rows-1;high[2*n+vcount:2*n+vcount+e]=cols-1;high[2*n+vcount+e:2*n+vcount+2*e]=rows-1
        integrality=np.zeros(size);integrality[2*n:2*n+vcount]=1;integrality[2*n+vcount+2*e:]=1
        for key,pos in fixed.items():
            i=index[key];low[i]=high[i]=pos[0];low[n+i]=high[n+i]=pos[1]
        ri,cj,data,lower,upper=[],[],[],[],[]
        def add(terms,lo=-np.inf,hi=np.inf):
            row=len(lower)
            for j,v in terms:ri.append(row);cj.append(j);data.append(v)
            lower.append(lo);upper.append(hi)
        bynode=defaultdict(list);bypoint=defaultdict(list)
        for j,(key,point) in enumerate(variables):bynode[key].append((2*n+j,point));bypoint[point].append(2*n+j)
        for key,options in bynode.items():
            i=index[key];add([(j,1) for j,p in options],1,1)
            for axis,offset in ((0,0),(1,n)):add([(i+offset,1)]+[(j,-p[axis]) for j,p in options],0,0)
        for point,options in bypoint.items():add([(j,1) for j in options],hi=1)
        point_keys=list(bypoint)
        if point_keys:
            for a,b in cKDTree(np.asarray(point_keys)).query_pairs(.16-EPS):add([(j,1) for j in bypoint[point_keys[a]]+bypoint[point_keys[b]]],hi=1)
        for k,(a,b) in enumerate(edge_indices):
            for axis,offset in ((0,0),(1,n)):
                d=2*n+vcount+axis*e+k
                add([(a+offset,1),(b+offset,-1),(d,-1)],hi=0);add([(a+offset,-1),(b+offset,1),(d,-1)],hi=0)
        for k,g in enumerate(groups):
            ex,ext,inside=(index[g[name]] for name in ('exchange','external','internal'))
            z=2*n+vcount+2*e+4*k
            add([(z+j,1) for j in range(4)],1,1)
            # z selects left, right, bottom or top. Equal transverse coordinates
            # and strict inward order hold for the selected face.
            for face in range(4):
                horizontal=face<2;axis=0 if horizontal else 1;trans=n if horizontal else 0
                along=0 if horizontal else n;extent=cols-1 if horizontal else rows-1;big=max(cols,rows)+1
                for a,b in ((ex,ext),(ext,inside)):
                    add([(a+trans,1),(b+trans,-1),(z+face,big)],hi=big)
                    add([(b+trans,1),(a+trans,-1),(z+face,big)],hi=big)
                    sign=1 if face in (0,2) else -1
                    add([(a+along,sign),(b+along,-sign),(z+face,big)],hi=big-.16)
                target=0 if face in (0,2) else extent
                add([(ex+along,1),(z+face,big)],hi=big+target)
                add([(ex+along,-1),(z+face,big)],hi=big-target)
        add([(j,float(v)) for j,v in enumerate(cost) if v],hi=incumbent+EPS)
        matrix=coo_matrix((data,(ri,cj)),shape=(len(lower),size)).tocsc()
        result=milp(cost,integrality=integrality,bounds=Bounds(low,high),constraints=LinearConstraint(matrix,lower,upper),options={'time_limit':time_limit,'mip_rel_gap':0.})
        if result.x is not None:
            proposal={k:tuple(fixed[k]) if k in fixed else variables[max((j for j,(name,point) in enumerate(variables) if name==k),key=lambda j:result.x[2*n+j])][1] for k in ids}
            if all(triplet_face(g,proposal,cols,rows) for g in groups) and value(proposal)<=incumbent+EPS:
                xy=proposal;report['selected_solution']='HiGHS MILP incumbent'
        bound=getattr(result,'mip_dual_bound',None)
        bound=max(0.,float(bound)) if bound is not None and np.isfinite(bound) else 0.
        report.update(method='Finite-cell MILP / SciPy HiGHS with fixed TCA anchors',solver_status=int(result.status),solver_message=result.message,variables=size,
                      objective=value(xy),global_lower_bound=bound,relative_gap=max(0.,value(xy)-bound)/max(1.,value(xy)),optimality_proven=bool((result.status==0 and report['selected_solution']=='HiGHS MILP incumbent') or value(xy)<=bound+EPS))
    # A final check never relaxes a failed hard constraint.
    if any(tuple(xy[k]) not in choices[k] for k in ids):raise RuntimeError('Forbidden node position')
    if len({tuple(p) for p in xy.values()})!=n:raise RuntimeError('Node positions overlap')
    if n>1 and cKDTree(np.asarray(list(xy.values()))).query_pairs(.16-EPS):raise RuntimeError('Node clearance violated')
    report['position_constraints']['exchange_alignment']=validate_alignment(alignment_groups,xy,cols,rows)
    report['distinct_grid_points_verified']=True
    report['fractional_fixed_node_ids']=sorted(k for k,p in xy.items() if not all(abs(v-round(v))<EPS for v in p))
    return {k:list(v) for k,v in xy.items()},report
