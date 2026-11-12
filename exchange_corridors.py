"""Reserve true exchange/transport incidences as straight, inward corridors."""
import math
from collections import defaultdict,deque
from crossing_routing import EPS,_on,measure,validate_boundary,intersections
from tca_curves import arc_controls
from tca_geometry import allowed_points
from exchange_alignment import _compartment,_base,_identities


def transport_paths(model,nodes,external,internal):
    """Shortest real same-compound transport chains, including periplasm."""
    metabolites={m['id']:m for m in model['metabolites']}
    compartments={mid:(_compartment(m) or str(m.get('compartment','')).lower()) for mid,m in metabolites.items()}
    reference=metabolites[external];identity=_identities(reference);base=_base(reference,_compartment(reference))
    chemical={external,internal}
    chemical.update(mid for mid,m in metabolites.items() if 'met:'+mid in nodes and
                    ((identity and identity&_identities(m)) or
                     (base is not None and _base(m,compartments[mid])==base)))
    graph=defaultdict(list)
    for reaction in model['reactions']:
        if 'rxn:'+reaction['id'] not in nodes:continue
        mids=sorted(chemical&set(reaction['metabolites']))
        for a in mids:
            for b in mids:
                if compartments[a]!=compartments[b] and reaction['metabolites'][a]*reaction['metabolites'][b]<0:
                    graph[a].append((b,reaction['id']))
    queue=deque([(external,['met:'+external],{external})]);found=[];shortest=None
    while queue:
        mid,path,visited=queue.popleft()
        if shortest is not None and len(path)>shortest:break
        if mid==internal:
            shortest=len(path);found.append(path)
            if len(found)>=32:break
            continue
        for target,rid in sorted(graph[mid]):
            if target not in visited and 'rxn:'+rid not in path:queue.append((target,[*path,'rxn:'+rid,'met:'+target],visited|{target}))
    return found


def corridor_incidences(group):
    keys=[(group['reaction_id'],group['external'][4:])]
    path=group.get('transport_path')
    if path:
        for i in range(1,len(path),2):keys.extend((path[i][4:],path[j][4:]) for j in (i-1,i+1))
    elif group.get('primary_transport'):
        keys.extend((group['primary_transport'],group[k][4:]) for k in ('external','internal'))
    return keys


def discover_corridors(model,nodes,exchange_nodes,groups):
    metabolites={m['id']:m for m in model['metabolites']};reactions={r['id']:r for r in model['reactions']}
    matched={g['exchange']:g for g in groups};out=[]
    for key in sorted(exchange_nodes):
        reaction=reactions[nodes[key]['original_id']]
        external=[mid for mid in reaction['metabolites'] if 'met:'+mid in nodes and _compartment(metabolites[mid])=='e']
        if not external:continue  # Cytosolic sinks are not EX_a_e / a_e pairs.
        if len(external)!=1:raise ValueError('Cannot uniquely identify exchange extracellular partner: '+reaction['id'])
        mid=external[0];group=dict(matched.get(key,dict(reaction_id=reaction['id'],exchange=key,external='met:'+mid)))
        group['transports']=[]
        if 'internal' not in group:
            e=metabolites[mid];identity=_identities(e);base=_base(e,'e')
            candidates=[m['id'] for m in metabolites.values() if _compartment(m)=='c' and ((identity and identity&_identities(m)) or (base is not None and _base(m,'c')==base))]
            if candidates:raise ValueError('Exchange cytosolic identity is ambiguous or absent from display: '+reaction['id']+' '+str(candidates))
            group['kind']='pair';group['reason']='No corresponding cytosolic species exists in the source model'
        else:
            group['kind']='triplet';inside=group['internal'][4:]
            group['transport_paths']=transport_paths(model,nodes,mid,inside)
            group['transports']=sorted({path[i][4:] for path in group['transport_paths'] for i in range(1,len(path),2)})
            if not group['transport_paths']:raise ValueError('No displayed source transport path joins exchange triplet: '+reaction['id'])
        out.append(group)
    return out


def place_corridors(model,nodes,xy,backbones,exchange_nodes,groups,cols,rows):
    source_reactions={r['id']:r for r in model['reactions']}
    red_phase=backbones.get('red_first',{})
    fixed=set(backbones['fixed_node_offsets'])|{k for p in red_phase.get('routes',[]) for k in (p['source'],p['target'])}
    region=backbones['tca_exclusion']['region'] if backbones.get('tca_exclusion') else backbones.get('_placement_region')
    records=discover_corridors(model,nodes,exchange_nodes,groups);original={k:list(v) for k,v in xy.items()}
    curves=[]
    for r in backbones.get('routes',[]):
        if r['id']!='tca' or region is None or region['shape']!='circle':continue
        for step in r['steps']:
            for mid in (step['source'],step['target']):
                a,b='rxn:'+step['reaction_id'],'met:'+mid
                if a in xy and b in xy:curves.append(dict(rid=step['reaction_id'],mid=mid,source=a,target=b,priority=2,points=[xy[a],xy[b]],curve_controls=arc_controls(xy[a],xy[b],region['center'])))
    if red_phase:curves=red_phase['routes']
    reserved=[];reserved_keys=set();assignments=[]
    major={rid for r in backbones['routes'] for rid in r['reaction_ids']}
    ordered=sorted(records,key=lambda g:(g['exchange'] not in fixed,g.get('internal',g['external']) not in fixed,g['reaction_id']))
    xmin,xmax,ymin,ymax=backbones.get('routing_boundary',[0,cols-1,0,rows-1])
    sides=[('left',0,float(xmin),1.),('right',0,float(xmax),-1.),('bottom',1,float(ymin),1.),('top',1,float(ymax),-1.)]
    for group in ordered:
        ex,e=group['exchange'],group['external'];target=group.get('internal',e);base=xy[target]
        path=min(group.get('transport_paths',[]),key=lambda path:(len(path),sum(len(source_reactions[path[i][4:]]['metabolites']) for i in range(1,len(path),2)),path)) if group.get('transport_paths') else []
        primary=path[1][4:] if len(path)==3 else None
        chain_ids=[ex,*path] if path else [ex,e]
        # Major cycle nodes stay fixed; ordinary cytosolic partners may move.
        deltas=[(0,0)] if target in fixed else sorted([(dx,dy) for dx in range(-4,5) for dy in range(-4,5) if abs(dx)+abs(dy)<=4],key=lambda d:(abs(d[0])+abs(d[1]),d))
        candidates=[]
        for dx,dy in deltas:
            c=[base[0]+dx,base[1]+dy]
            if not (xmin+.16<c[0]<xmax-.16 and ymin+.16<c[1]<ymax-.16):continue
            if target not in (region or {}).get('node_ids',[]) and not allowed_points([c],region)[0]:continue
            if target not in reserved_keys and any(_on(c,p['points'][0],p['points'][-1]) for p in reserved):continue
            for face,axis,value,sign in sides:
                if backbones.get('adaptive_grid') and abs(c[1-axis]-round(c[1-axis]))>EPS:continue
                a=list(c);a[axis]=value;length=abs(c[axis]-value)
                if length<(.64 if path else .16):continue
                proposed={ex:a,target:c}
                if path:
                    if backbones.get('adaptive_grid') and length<len(chain_ids)-1-EPS:continue
                    for index,key in enumerate(chain_ids[1:-1],1):
                        point=list(a)
                        distance=float(index) if backbones.get('adaptive_grid') else min(.75*index,length*index/(len(chain_ids)-1))
                        point[axis]+=sign*distance
                        proposed[key]=list(original[key]) if key in fixed else point
                    chain=[proposed[k] for k in chain_ids]
                    if any(abs(p[1-axis]-a[1-axis])>EPS for p in chain):continue
                    if any(sign*(v[axis]-u[axis])<.16-EPS for u,v in zip(chain,chain[1:])):continue
                if any(k in fixed and math.dist(v,original[k])>EPS for k,v in proposed.items()):continue
                if any(not allowed_points([p],region)[0] for k,p in proposed.items() if k not in (region or {}).get('node_ids',[])):continue
                if any(math.dist(p,q)<.16-EPS for i,p in enumerate(proposed.values()) for q in list(proposed.values())[i+1:]):continue
                if any(math.dist(p,xy[k])<.16-EPS for p in proposed.values() for k in reserved_keys|fixed if k not in proposed):continue
                corridor=dict(rid=group['reaction_id'],mid=target,source=ex,target=target,priority=2 if group['reaction_id'] in major else 1,points=[a,c])
                if backbones.get('adaptive_grid'):
                    from adaptive_routing import segment_clearance,incident_exempt
                    if segment_clearance(a,c,region)<1-EPS and not incident_exempt(a,c,corridor,dict(xy,**proposed),region):continue
                if any(_on(xy[k],a,c) for k in fixed|reserved_keys if k not in proposed):continue
                own={(group['reaction_id'],e[4:])}
                if path:own.update(corridor_incidences(dict(group,transport_path=path)))
                obstacles=[p for p in curves if (p['rid'],p['mid']) not in own]
                def clear(other):
                    if backbones.get('adaptive_grid') and not other.get('curve_controls'):
                        start,end=corridor['points']
                        along=0 if abs(start[1]-end[1])<EPS else 1
                        low,high=sorted((start[along],end[along]));level=start[1-along]
                        for u,v in zip(other['points'],other['points'][1:]):
                            axis=0 if abs(u[1]-v[1])<EPS else 1
                            if axis!=along or abs(level-u[1-axis])>=1-EPS:continue
                            lo,hi=sorted((u[axis],v[axis]))
                            if min(high,hi)-max(low,lo)>=2-EPS:return False
                    if measure([corridor,other])['hard_verified']:return True
                    # Adaptive routing audits red/orange point crossings after
                    # placement. Do not reject its straight candidate here;
                    # positive-length overlaps remain forbidden.
                    if not (backbones.get('adaptive_grid') and corridor['priority']==1 and other['priority']==2):return False
                    if other.get('curve_controls'):return False
                    for a,b in intersections(corridor,other):
                        if math.dist(a,b)>=EPS:return False
                        orientations={0 if abs(u[1]-v[1])<EPS else 1
                                      for u,v in zip(other['points'],other['points'][1:]) if _on(a,u,v)}
                        if len(orientations)!=1:return False
                    return True
                if any(not clear(p) for p in [*obstacles,*reserved]):continue
                displacement=abs(dx)+abs(dy)
                cost=10*displacement+length+.001*sum(math.dist(original.get(k,p),p) for k,p in proposed.items())
                candidates.append((cost,face,proposed,corridor,primary))
        if not candidates:raise ValueError('No straight exchange corridor satisfies backbone/boundary constraints: '+group['reaction_id'])
        _,face,proposed,corridor,primary=min(candidates,key=lambda x:(x[0],x[1]))
        xy.update(proposed);reserved_keys.update(proposed);reserved.append(corridor)
        assignments.append(dict(group,face=face,primary_transport=primary,transport_path=path,
                                chain_node_ids=chain_ids,node_ids=list(proposed),corridor=[list(p) for p in corridor['points']]))
    # Clear unrelated ordinary nodes from the straight corridors. Identity and
    # all fixed backbone coordinates are retained; cofactors are attached later.
    occupied={tuple(xy[k]) for k in fixed|reserved_keys}
    for key in sorted(xy):
        if key in fixed|reserved_keys or nodes[key]['kind']=='currency':continue
        point=xy[key]
        if key in exchange_nodes:
            on_boundary=(abs(point[0]-xmin)<EPS or abs(point[0]-xmax)<EPS or abs(point[1]-ymin)<EPS or abs(point[1]-ymax)<EPS)
            inside=xmin-EPS<=point[0]<=xmax+EPS and ymin-EPS<=point[1]<=ymax+EPS
            if not (inside and on_boundary):
                choices=[]
                boundary_points=([(x,float(y)) for x in (xmin,xmax) for y in range(math.floor(ymin)+1,math.ceil(ymax))]+
                                 [(float(x),y) for y in (ymin,ymax) for x in range(math.floor(xmin)+1,math.ceil(xmax))])
                for q in boundary_points:
                    if not allowed_points([q],region)[0] or any(math.dist(q,p)<.16-EPS for p in occupied):continue
                    if any(_on(q,a,b) for p in [*reserved,*curves] if not p.get('curve_controls') for a,b in zip(p['points'],p['points'][1:])):continue
                    choices.append((abs(q[0]-point[0])+abs(q[1]-point[1]),q))
                if not choices:raise ValueError('No boundary position for cytosolic sink: '+key)
                xy[key]=list(min(choices)[1])
            occupied.add(tuple(xy[key]));continue
        conflict=not (xmin+.16<point[0]<xmax-.16 and ymin+.16<point[1]<ymax-.16) or any(_on(point,a,b) for p in [*reserved,*curves] if not p.get('curve_controls') for a,b in zip(p['points'],p['points'][1:])) or any(math.dist(point,q)<.16-EPS for q in occupied)
        if conflict:
            candidates=[]
            for x in range(math.floor(xmin)+1,math.ceil(xmax)):
                for y in range(math.floor(ymin)+1,math.ceil(ymax)):
                    q=[float(x),float(y)]
                    if not allowed_points([q],region)[0] or any(math.dist(q,p)<.16-EPS for p in occupied):continue
                    if any(_on(q,a,b) for p in [*reserved,*curves] if not p.get('curve_controls') for a,b in zip(p['points'],p['points'][1:])):continue
                    candidates.append((abs(q[0]-point[0])+abs(q[1]-point[1]),q))
            if not candidates:raise ValueError('No node position outside reserved exchange corridors: '+key)
            xy[key]=min(candidates)[1]
        occupied.add(tuple(xy[key]))
    changed={k:dict(before=original[k],after=list(xy[k])) for k in xy if math.dist(original[k],xy[k])>EPS}
    return dict(groups=assignments,changed_nodes=changed,straight_node_ids=sorted(reserved_keys),verified=True)


def lock_corridor_paths(paths,specification):
    lookup={(p['rid'],p['mid']):p for p in paths};locked=[]
    for group in specification['groups']:
        keys=corridor_incidences(group)
        for key in keys:
            if key not in lookup:raise ValueError('Source incidence missing from exchange corridor: '+str(key))
            p=lookup[key];a,b=p['points'][0],p['points'][-1]
            if abs(a[0]-b[0])>EPS and abs(a[1]-b[1])>EPS:raise ValueError('Exchange corridor endpoints are not collinear: '+str(key))
            if p.get('red_fixed') and (any(not _on(q,a,b) for q in p['points']) or abs(sum(math.dist(u,v) for u,v in zip(p['points'],p['points'][1:]))-math.dist(a,b))>1e-6):
                raise ValueError('Frozen red route is incompatible with straight exchange corridor: '+str(key))
            p['points']=[a,b];p['straight_locked']=True;locked.append(p)
    return locked


def validate_corridors(paths,specification):
    lookup={(p['rid'],p['mid']):p for p in paths};count=0
    for group in specification['groups']:
        keys=corridor_incidences(group)
        a,b=group['corridor']
        for key in keys:
            p=lookup[key]
            if p.get('curve_controls') or any(not _on(q,a,b) for q in p['points']):raise ValueError('Saved exchange corridor is not straight: '+str(key))
            length=sum(math.dist(u,v) for u,v in zip(p['points'],p['points'][1:]))
            if abs(length-math.dist(p['points'][0],p['points'][-1]))>1e-6:raise ValueError('Saved exchange corridor backtracks: '+str(key))
            count+=1
    return dict(verified=True,groups=len(specification['groups']),triplets=sum('internal' in g for g in specification['groups']),pairs=sum('internal' not in g for g in specification['groups']),
                multistep_groups=sum(len(g.get('transport_path',[]))>3 for g in specification['groups']),straight_incidences=count)
