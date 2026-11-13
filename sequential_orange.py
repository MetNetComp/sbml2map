"""Place exchange branches against a frozen, growing colored network.

Each branch uses a layered grid-flow MILP. Layer changes place real biological
nodes; grid arcs draw real incidences. Existing colored geometry is immutable.
"""
import copy
import hashlib
import json
import tempfile
from pathlib import Path
import math
import time
from collections import defaultdict, deque

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

from crossing_routing import EPS, clean, _on
from red_objective import route_length


def organic_graph(nodes, edges):
    blocked = {k for k, n in nodes.items() if n['kind'] == 'currency' or n.get('shared_cofactor')}
    adjacency = defaultdict(set)
    lookup = {}
    for i, edge in enumerate(edges):
        a, b = edge['source'], edge['target']
        if a in blocked or b in blocked:
            continue
        adjacency[a].add(b); adjacency[b].add(a)
        lookup[frozenset((a, b))] = i
    return adjacency, lookup


def hop_distances(adjacency, roots):
    distance = {k: 0 for k in roots}
    queue = deque(sorted(roots))
    while queue:
        key = queue.popleft()
        for other in sorted(adjacency[key]):
            if other not in distance:
                distance[other] = distance[key] + 1
                queue.append(other)
    return distance


def candidate_paths(adjacency, terminal, connected, exchanges, limit=24):
    """Enumerate shortest-hop organic branches to any existing colored node.

    The geometry optimum is certified within these reported topology candidates,
    not over all longer biological paths or all sequential branch orders.
    """
    distance = hop_distances(adjacency, connected)
    if terminal not in distance:
        return []
    out = []; queue=deque([[terminal]])
    max_hops=distance[terminal]+4
    while queue and len(out)<limit:
        path=queue.popleft();key=path[-1]
        if key in connected:
            out.append(path);continue
        for other in sorted(adjacency[key]):
            if other in path or (other in exchanges and other!=terminal):continue
            if len(path)+distance.get(other,math.inf)>max_hops:continue
            queue.append(path+[other])
    return out



def branch_milp(chain, nodes, edges, lookup, fixed_xy, frozen, boundary, region, time_limit=20., auxiliary=False, fanout_extent=1.5):
    """Joint node placement and orthogonal routing on a finite primary grid.

    Positive grid arc costs and unit flow give a single branch. Zero-cost
    layer changes assign intermediate biological nodes to integer intersections.
    One change per intersection prevents coincident biological nodes.
    """
    from adaptive_routing import integer, in_common_port
    from sequential_grid import Grid
    endpoint = chain[-1]
    grid = Grid(fixed_xy, boundary, region, auxiliary=auxiliary, fanout_nodes=[endpoint] if auxiliary else (), fanout_extents={endpoint:fanout_extent}, extra_points=[q for p in frozen for q in p['points']])
    for path in frozen:
        grid.reserve(path)
    size = grid.nx * grid.ny
    layers = len(chain) - 1
    source = layers * size
    sink = source + 1
    tails = []; heads = []; costs = []; labels = []
    def arc(a, b, cost, label):
        tails.append(a); heads.append(b); costs.append(cost); labels.append(label)
    boundary_points = set()
    for v in range(size):
        x, y = grid.point(v)
        if integer(x) and integer(y) and (x in boundary[:2] or y in boundary[2:]) and v not in grid.owners:
            boundary_points.add(v)
    for v in sorted(boundary_points):
        arc(source, v, 0., ('start', v))
    occupied = set(grid.owners)
    target = grid.lookup[endpoint]
    switches = defaultdict(list)
    move_groups = defaultdict(list)
    for layer in range(layers):
        edge = edges[lookup[frozenset(chain[layer:layer+2])]]
        p = dict(rid=edge['reaction_id'], mid=edge['metabolite_id'], source=chain[layer], target=chain[layer+1], priority=1)
        h, vmask = grid.route_edges(p)
        for a in range(size):
            if a in occupied and not (layer == layers-1 and a == target):
                continue
            if a == target:
                continue
            iy, ix = divmod(a, grid.nx)
            options = []
            if ix and h[iy, ix-1]: options.append(a-1)
            if ix+1 < grid.nx and h[iy, ix]: options.append(a+1)
            if iy and vmask[iy-1, ix]: options.append(a-grid.nx)
            if iy+1 < grid.ny and vmask[iy, ix]: options.append(a+grid.nx)
            for b in options:
                if b in occupied and not (layer == layers-1 and b == target):
                    continue
                pa, pb = grid.point(a), grid.point(b)
                # Only the first layer may leave the EX perimeter.
                if a in boundary_points and layer != 0:
                    continue
                if b in boundary_points:
                    continue
                key = tuple(sorted((a, b)))
                blockers = [*grid.orange_edges.get(key, ()), *grid.red_edges.get(key, ())]
                if any(not (layer == layers-1 and in_common_port(pa, pb, p, q, fixed_xy)) for q in blockers):
                    continue
                # Red point crossings are disallowed in this sequential phase.
                if any(z != target and grid.red_vertices.get(z) and not (layer == layers-1 and all(in_common_port(grid.point(z),grid.point(z),p,q,fixed_xy) for q,_ in grid.red_vertices[z])) for z in (a,b)):
                    continue
                # Compare full reserved runs, including fractional TCA attachments.
                axis = 0 if abs(pa[1]-pb[1]) < EPS else 1
                level = pa[1-axis]
                lo, hi = sorted((pa[axis], pb[axis]))
                if any(abs(level-other) < 1-EPS and min(hi,right)-max(lo,left) >= 2-EPS
                       for other,left,right,q in grid.parallel_intervals.get((axis,math.floor(level)), ())):
                    continue
                move_groups[key].append(len(costs))
                arc(layer*size+a, layer*size+b, math.dist(pa,pb), ('move',layer,a,b))
        if layer+1 < layers:
            for a in range(size):
                point = grid.point(a)
                if a in occupied or a in boundary_points or not all(integer(x) for x in point):
                    continue
                if any(_on(point, x, y) for path in frozen for x,y in zip(path['points'],path['points'][1:])):
                    continue
                switches[a].append(len(costs))
                arc(layer*size+a, (layer+1)*size+a, 0., ('switch',layer,a))
    arc((layers-1)*size+target, sink, 0., ('end',target))
    n = len(costs)
    rows = []; columns = []; values = []
    lower = np.zeros(sink+1).tolist(); upper = lower.copy()
    lower[source] = upper[source] = 1.
    lower[sink] = upper[sink] = -1.
    for i,(a,b) in enumerate(zip(tails,heads)):
        rows.extend((a,b)); columns.extend((i,i)); values.extend((1.,-1.))
    for indices in [*switches.values(), *move_groups.values()]:
        row = len(lower); lower.append(-np.inf); upper.append(1.)
        rows.extend([row]*len(indices)); columns.extend(indices); values.extend([1.]*len(indices))
    matrix = coo_matrix((values,(rows,columns)), shape=(len(lower),n)).tocsc()
    deadline=time.monotonic()+time_limit
    for iteration in range(32):
        result = milp(np.array(costs),integrality=np.ones(n),bounds=Bounds(0.,1.),
                      constraints=LinearConstraint(matrix,lower,upper),
                      options=dict(time_limit=max(.1,deadline-time.monotonic()),mip_rel_gap=0.))
        if result.x is None:
            if not auxiliary:
                return branch_milp(chain,nodes,edges,lookup,fixed_xy,frozen,boundary,region,time_limit,True)
            if auxiliary and fanout_extent<3:
                return branch_milp(chain,nodes,edges,lookup,fixed_xy,frozen,boundary,region,time_limit,True,3.)
            return None, dict(status=int(result.status),message=result.message)
        outgoing=defaultdict(list)
        for i,a in enumerate(tails):
            if result.x[i]>.5:outgoing[a].append(i)
        queue=deque([source]);parent={source:None}
        while queue and sink not in parent:
            a=queue.popleft()
            for i in outgoing[a]:
                b=heads[i]
                if b not in parent:parent[b]=(a,i);queue.append(b)
        if sink not in parent:raise ValueError('MILP incumbent has no integral source-to-sink path')
        following={};cursor=sink
        while cursor!=source:
            a,i=parent[cursor];following[a]=i;cursor=a
        current=source; positions={}; pieces=[[] for _ in range(layers)]; visited=set()
        while current != sink:
            if current in visited or current not in following:
                raise ValueError('Invalid sequential MILP flow')
            visited.add(current); i=following[current]; label=labels[i]
            if label[0]=='start':
                positions[chain[0]]=list(grid.point(label[1]))
            elif label[0]=='move':
                _,layer,a,b=label
                if not pieces[layer]:pieces[layer].append(grid.point(a))
                pieces[layer].append(grid.point(b))
            elif label[0]=='switch':
                _,layer,a=label; positions[chain[layer+1]]=list(grid.point(a))
            current=heads[i]
        positions[endpoint]=list(fixed_xy[endpoint]); routes=[]
        for layer,points in enumerate(pieces):
            if len(points)<2:raise ValueError('Zero-length biological incidence in MILP branch')
            edge=edges[lookup[frozenset(chain[layer:layer+2])]]
            points=clean(points)
            if edge['source'] != chain[layer]:points=list(reversed(points))
            routes.append(dict(rid=edge['reaction_id'],mid=edge['metabolite_id'],source=edge['source'],target=edge['target'],
                               priority=1,points=points,sequential_fixed=True,auxiliary_lane=auxiliary,local_fanout=auxiliary,fanout_node_extents={endpoint:fanout_extent} if auxiliary else {}))
        branch=dict(positions=positions,routes=routes,length=sum(route_length(p) for p in routes))
        solver=dict(
            status=int(result.status),message=result.message,optimality_proven=result.status==0,
            objective=float(sum(route_length(p) for p in routes)),lower_bound=float(getattr(result,'mip_dual_bound',result.fun)),
            relative_gap=float(getattr(result,'mip_gap',0.)),variables=n)
        from adaptive_routing import validate_paths,parallel_conflicts
        try:
            validate_paths(frozen+routes,dict(fixed_xy,**positions),boundary,region,[])
        except ValueError as error:
            if time.monotonic()>=deadline:
                if auxiliary and fanout_extent<3:
                    return branch_milp(chain,nodes,edges,lookup,fixed_xy,frozen,boundary,region,time_limit,True,3.)
                return None,dict(solver,geometry_rejection=str(error))
            # Cut only simultaneous geometric runs that violate the spacing rule.
            # The same invalid run remains forbidden if nodes move along that run.
            by_key={(p['rid'],p['mid']):p for p in frozen+routes}
            new_keys={(p['rid'],p['mid']) for p in routes}
            cuts=[]
            for conflict in parallel_conflicts(frozen+routes):
                p,q=(by_key[tuple(conflict[k])] for k in ('first','second'))
                for a,b in zip(clean(p['points']),clean(p['points'])[1:]):
                    axis=0 if abs(a[1]-b[1])<EPS else 1
                    for c,d in zip(clean(q['points']),clean(q['points'])[1:]):
                        other=0 if abs(c[1]-d[1])<EPS else 1
                        if axis!=other or abs(a[1-axis]-c[1-axis])>=1-EPS:continue
                        if min(max(a[axis],b[axis]),max(c[axis],d[axis]))-max(min(a[axis],b[axis]),min(c[axis],d[axis]))<2-EPS:continue
                        segments=[]
                        if (p['rid'],p['mid']) in new_keys:segments.append((a,b))
                        if (q['rid'],q['mid']) in new_keys:segments.append((c,d))
                        physical=[key for key in move_groups if any(_on(grid.point(key[0]),u,v) and _on(grid.point(key[1]),u,v) for u,v in segments)]
                        if physical:cuts.append(([i for key in physical for i in move_groups[key]],len(physical)-1))
            if not cuts:
                active=[i for i in following.values()]
                cuts=[(active,len(active)-1)]
            from scipy.sparse import vstack
            extra_rows=[];extra_cols=[];extra_values=[];bounds=[]
            for row,(indices,bound) in enumerate(cuts):
                extra_rows.extend([row]*len(indices));extra_cols.extend(indices);extra_values.extend([1.]*len(indices));bounds.append(bound)
            extra=coo_matrix((extra_values,(extra_rows,extra_cols)),shape=(len(cuts),n)).tocsc()
            matrix=vstack([matrix,extra],format='csc')
            lower=list(lower)+[-np.inf]*len(cuts);upper=list(upper)+bounds
            continue
        solver['geometry_cut_rounds']=iteration
        solver['grid_domain']='local fan-out' if auxiliary else 'primary grid'
        endpoint_xy=fixed_xy[endpoint]
        hop_bound=len(chain)-1 if all(abs(x-round(x))<EPS for x in endpoint_xy) else len(chain)-2
        global_bound=max(hop_bound,min(abs(endpoint_xy[0]-boundary[0]),abs(endpoint_xy[0]-boundary[1]),abs(endpoint_xy[1]-boundary[2]),abs(endpoint_xy[1]-boundary[3])))
        if not auxiliary and branch['length']>global_bound+EPS:
            other,other_solver=branch_milp(chain,nodes,edges,lookup,fixed_xy,frozen,boundary,region,time_limit,True)
            domains=[dict(solver),dict(other_solver)]
            if other is not None and other['length']<branch['length']-EPS:
                branch,solver=other,other_solver
            solver['evaluated_grid_domains']=domains
            solver['optimality_proven']=other is not None and other_solver.get('optimality_proven',False)
        return branch,solver
    if auxiliary and fanout_extent<3:
        return branch_milp(chain,nodes,edges,lookup,fixed_xy,frozen,boundary,region,time_limit,True,3.)
    return None,dict(status=1,message='Geometry cut iteration limit')



def place_exchange_branches(nodes, edges, xy, red, exchanges, boundary, region, progress=print, time_limit=20.):
    from adaptive_routing import validate_paths
    adjacency,lookup=organic_graph(nodes,edges)
    roots={k for p in red for k in (p['source'],p['target']) if not nodes[k].get('shared_cofactor') and nodes[k]['kind']!='currency'}
    distance=hop_distances(adjacency,roots)
    order=sorted(exchanges,key=lambda k:(distance.get(k,math.inf),k))
    connected=set(roots); frozen=copy.deepcopy(red)
    fixed_xy={k:list(xy[k]) for p in red for k in (p['source'],p['target'])}
    chosen=set(); records=[]; disconnected=[]
    for terminal in order:
        if terminal in connected:
            records.append(dict(exchange=terminal,red_hops=distance[terminal],new_length_grid=0.,already_connected=True));continue
        chains=candidate_paths(adjacency,terminal,connected,exchanges)
        if not chains:
            disconnected.append(terminal);continue
        progress(f'Sequential orange: {terminal}, red hops {distance[terminal]}, {len(chains)} topology candidates',flush=True)
        best=None; attempts=[]
        for chain in chains:
            endpoint=fixed_xy[chain[-1]]
            hop_bound=len(chain)-1 if all(abs(x-round(x))<EPS for x in endpoint) else len(chain)-2
            bound=max(hop_bound,min(abs(endpoint[0]-boundary[0]),abs(endpoint[0]-boundary[1]),abs(endpoint[1]-boundary[2]),abs(endpoint[1]-boundary[3])))
            if best is not None and bound >= best[0]['length']-EPS:
                attempts.append(dict(path=chain,skipped='Proven lower bound cannot improve incumbent',lower_bound=bound));continue
            cache=Path(tempfile.gettempdir())/'sbml2map-sequential-orange-cache'
            cache.mkdir(exist_ok=True)
            digest=hashlib.sha256(json.dumps([chain,fixed_xy,frozen,boundary,region,'flow-v7-local-fan' ],sort_keys=True).encode()).hexdigest()
            cache_file=cache/(digest+'.json')
            if cache_file.exists():
                cached=json.loads(cache_file.read_text());branch,solver=cached['branch'],cached['solver']
            else:
                branch,solver=branch_milp(chain,nodes,edges,lookup,fixed_xy,frozen,boundary,region,time_limit)
                if branch is not None:cache_file.write_text(json.dumps(dict(branch=branch,solver=solver)))
            if branch is not None:
                trial_xy=dict(fixed_xy,**branch['positions'])
                try:
                    validate_paths(frozen+branch['routes'],trial_xy,boundary,region,[])
                except ValueError as error:
                    solver['geometry_rejection']=str(error);branch=None
            attempts.append(dict(path=chain,solver=solver))
            if branch is not None and (best is None or branch['length']<best[0]['length']-EPS):
                best=(branch,chain,solver)
        if best is None:raise ValueError('No sequential exchange branch: '+terminal+' '+str(attempts))
        branch,chain,solver=best
        target=chain[-1]; target_color='red' if target in roots else 'orange'
        before={k:list(v) for k,v in fixed_xy.items()}
        fixed_xy.update(branch['positions']);xy.update(branch['positions'])
        assert all(fixed_xy[k]==v for k,v in before.items())
        frozen.extend(branch['routes']);connected.update(chain)
        chosen.update(lookup[frozenset((a,b))] for a,b in zip(chain,chain[1:]))
        records.append(dict(exchange=terminal,red_hops=distance[terminal],red_intermediate_nodes=max(0,distance[terminal]-1),
                            attachment_node=target,attachment_color=target_color,path=chain,new_length_grid=branch['length'],
                            solver=solver,candidates=attempts,existing_colored_geometry_preserved=True))
        progress(f'Sequential orange: attached {terminal} to {target_color} {target}, new length {branch["length"]:g}',flush=True)
    return dict(routes=frozen,positions=fixed_xy,edge_indices=sorted(chosen),records=records,
                disconnected_exchange_node_ids=disconnected,root_node_ids=sorted(roots),
                policy='red_hop_order_then_shortest_new_branch_to_red_or_orange',
                ordering='Initial organic edge count to frozen red, then exchange ID; never recomputed against orange',
                objective='New drawn branch length only; existing colored length is excluded',
                topology_search='Up to 24 simple paths, at most four hops beyond the current minimum to the colored network',
                geometry_search='Layered grid-flow MILP; primary grid plus short endpoint fan-out, current fixed frame; no red point crossings',
                global_optimality_proven=False,coenzyme_transit_allowed=False)


def sequential_layout(nodes,edges,xy,backbones,exchanges,major,cols,rows,core,red_report,progress=print):
    from adaptive_routing import (attach_grid_cofactors,make_paths,route_paths,validate_paths,
                                  validate_node_grid,perimeter_distance,POLICY)
    from red_objective import preservation_report
    from tca_geometry import allowed_points,validate_exclusion
    original=copy.deepcopy(xy);failures=[]
    for expansion in (0,1,2,3,4,6):
        trial=copy.deepcopy(core['backbones']);positions=copy.deepcopy(original);positions.update(copy.deepcopy(core['positions']))
        boundary=[0,cols-1+math.ceil((cols-1)*expansion/2),-math.ceil((rows-1)*expansion/2),rows-1]
        region=trial.get('tca_exclusion',{}).get('region')
        try:
            selected=place_exchange_branches(nodes,edges,positions,core['routes'],exchanges,boundary,region,progress)
            fixed=selected['positions'];protected=set(fixed);occupied={tuple(p) for p in fixed.values()}
            paths=selected['routes']
            # Unconnected/inorganic exchanges remain gray, on unused perimeter cells.
            pending=sorted(k for k in nodes if k not in protected and nodes[k]['kind']!='currency')
            pending.sort(key=lambda k:(k not in exchanges,k))
            for key in pending:
                old=positions[key];choices=[]
                for y in range(math.ceil(boundary[2]),math.floor(boundary[3])+1):
                    for x in range(math.ceil(boundary[0]),math.floor(boundary[1])+1):
                        point=(float(x),float(y));outer=x in boundary[:2] or y in boundary[2:]
                        if outer!=(key in exchanges) or point in occupied:continue
                        if perimeter_distance(point,region)<1-EPS or not allowed_points([point],region)[0]:continue
                        if any(_on(point,a,b) for p in paths for a,b in zip(p['points'],p['points'][1:])):continue
                        if any(math.dist(point,p)<1.-EPS for p in fixed.values()):continue
                        choices.append(point)
                if not choices:raise ValueError('No remaining grid cell for '+key)
                point=min(choices,key=lambda p:(math.dist(p,old),p));positions[key]=list(point);occupied.add(point)
            trial['exchange_connectors']=dict(color='#f28e2b',edge_indices=selected['edge_indices'],
                root_node_ids=selected['root_node_ids'],exchange_paths={r['exchange']:r.get('path',[]) for r in selected['records']},
                all_exchanges_connected=not selected['disconnected_exchange_node_ids'],
                disconnected_exchange_node_ids=selected['disconnected_exchange_node_ids'],
                method=selected['policy'],objective=selected['objective'],global_optimality_proven=False)
            trial.pop('exchange_corridors',None)
            attach_grid_cofactors(nodes,positions,region,boundary,paths)
            full=make_paths(nodes,edges,positions,trial,major)
            saved={(p['rid'],p['mid']):p for p in paths}
            for i,p in enumerate(full):
                key=(p['rid'],p['mid'])
                if key in saved:full[i]=copy.deepcopy(saved[key])
            remaining=[p for p in full if (p['rid'],p['mid']) not in saved]
            exceptions=route_paths(remaining,positions,boundary,region,progress,reservations=paths)
            report=validate_paths(full,positions,boundary,region,exceptions)
            assert all(positions[k]==v for k,v in fixed.items())
            assert all(p['points']==saved[p['rid'],p['mid']]['points'] for p in full if p.get('sequential_fixed'))
            preservation=preservation_report(core,positions,full)
            if not preservation['exactly_preserved']:raise ValueError('Frozen red geometry changed')
            report.update(node_grid=validate_node_grid(nodes,positions,region),red_routes_preserved=True,
                          red_preservation=preservation,sequential_orange_preserved=True)
            validate_exclusion(nodes,positions,region)
        except (ValueError,RuntimeError) as error:
            failures.append(dict(expansion=expansion,reason=str(error)))
            progress(f'Sequential orange expansion {expansion}: {error}',flush=True);continue
        selected.pop('routes');selected.pop('positions')
        trial['sequential_orange']=selected
        trial['red_first']=dict(routes=copy.deepcopy(core['routes']),report=dict(verified=True))
        trial['routing_boundary']=boundary;trial['render_frame']=boundary;trial['adaptive_grid']=True
        trial['crossing_routing']=dict(routes=full,report=report,signature=hashlib.sha256(json.dumps(full,sort_keys=True).encode()).hexdigest()[:24])
        trial['adaptive_layout']=dict(policy=POLICY,exchange_routing='sequential',expansion_steps=expansion,
            before_boundary=[0,cols-1,0,rows-1],after_boundary=boundary,rejected_layouts=failures,
            thin_lanes_used=report['thin_lanes_used'],compaction_skipped='Frozen orange positions and perimeter are preserved')
        red_report=copy.deepcopy(red_report)
        red_report.update(preservation=preservation,relaxation_level=0,relaxation_reason=None,
                          rejected_full_layouts=failures,accepted_reference_signature=core['signature'])
        trial['red_optimization']=red_report
        xy.clear();xy.update(positions);backbones.clear();backbones.update(trial)
        return backbones['adaptive_layout']
    raise ValueError('Sequential orange layout failed: '+str(failures))
