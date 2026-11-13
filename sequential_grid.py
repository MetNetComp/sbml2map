"""Static primary/auxiliary lattice and reserved geometry for exchange MILPs."""
import math
import numpy as np
from collections import defaultdict
from crossing_routing import EPS,clean
from adaptive_routing import PORT_RADIUS,PARALLEL_LENGTH,integer,incident_exempt

class Grid:
    def __init__(self,xy,boundary,region,auxiliary=False,fanout_nodes=(),local_keys=(),fanout_extents=None,extra_points=()):
        self.xy=xy;self.region=region;self.boundary=boundary;self.auxiliary=auxiliary;self.reserved=[];self.parallel_intervals=defaultdict(list)
        self.local_keys=set(local_keys)
        routing_xy={k:p for k,p in xy.items() if k not in self.local_keys}
        self.copy_h=defaultdict(list);self.copy_v=defaultdict(list)
        for key in self.local_keys:
            x,y=xy[key];self.copy_h[round(y,8)].append(x);self.copy_v[round(x,8)].append(y)
        for values in [*self.copy_h.values(),*self.copy_v.values()]:values.sort()
        self.fanout_nodes=set(fanout_nodes)
        self.fanout_extents={k:(fanout_extents or {}).get(k,1.5) for k in self.fanout_nodes}
        xmin,xmax,ymin,ymax=boundary
        xs=set(float(x) for x in range(math.ceil(xmin),math.floor(xmax)+1));ys=set(float(y) for y in range(math.ceil(ymin),math.floor(ymax)+1))
        if auxiliary:
            for key in self.fanout_nodes:
                cx,cy=xy[key]
                deltas=(-1.25,-.75,-.4,-.3,-.25,-.2,-.1,.1,.2,.25,.3,.4,.75,1.25)
                if self.fanout_extents[key]>1.5:deltas=(*deltas,-2.75,-2.25,-1.75,1.75,2.25,2.75)
                for delta in deltas:
                    if xmin<cx+delta<xmax:xs.add(round(cx+delta,9))
                    if ymin<cy+delta<ymax:ys.add(round(cy+delta,9))
        xs.update(float(p[0]) for p in routing_xy.values());ys.update(float(p[1]) for p in routing_xy.values())
        xs.update(float(p[0]) for p in extra_points);ys.update(float(p[1]) for p in extra_points)
        self.x=sorted(xs);self.y=sorted(ys);self.nx=len(self.x);self.ny=len(self.y)
        self.lookup={k:self.index(p) for k,p in routing_xy.items()};self.owners=defaultdict(set)
        for k,i in self.lookup.items():self.owners[i].add(k)
        self.red_edges=defaultdict(list);self.red_vertices=defaultdict(list);self.orange_edges=defaultdict(list)
        self.clearance={}
        self.hard_h=np.zeros((self.ny,self.nx-1),dtype=bool)
        self.hard_v=np.zeros((self.ny-1,self.nx),dtype=bool)
        self.node_occupied=np.zeros((self.ny,self.nx),dtype=bool)
        self.node_occupied.flat[list(self.lookup.values())]=True
        self._edge_masks={};self._attachment_masks={}
        self.horizontal_clearance,self.vertical_clearance=self._clearance_arrays()
    def _clearance_arrays(self):
        x=np.asarray(self.x);y=np.asarray(self.y)
        if not self.region:return np.full((self.ny,self.nx-1),np.inf),np.full((self.ny-1,self.nx),np.inf)
        def calculate(axis):
            low=x[:-1][None,:] if axis==0 else y[:-1][:,None]
            high=x[1:][None,:] if axis==0 else y[1:][:,None]
            level=y[:,None] if axis==0 else x[None,:]
            r=self.region
            if r['shape']=='circle':
                along=r['center'][axis];across=r['center'][1-axis]
                gap=np.maximum(np.maximum(low-along,along-high),0)
                minimum=np.sqrt(gap**2+(level-across)**2)
                maximum=np.sqrt(np.maximum((low-along)**2,(high-along)**2)+(level-across)**2)
                return np.maximum(np.maximum(minimum-r['radius'],r['radius']-maximum),0)
            l,t,side=r['left'],r['top'],r['side_length']
            along_low,along_high=(l,l+side) if axis==0 else (t-side,t)
            across_low,across_high=(t-side,t) if axis==0 else (l,l+side)
            gap=np.maximum(np.maximum(along_low-high,low-along_high),0)
            distances=[np.sqrt(gap**2+(level-edge)**2) for edge in (across_low,across_high)]
            across_gap=np.maximum(np.maximum(across_low-level,level-across_high),0)
            distances.extend(np.sqrt(np.maximum(np.maximum(low-edge,edge-high),0)**2+across_gap**2) for edge in (along_low,along_high))
            return np.minimum.reduce(distances)
        return calculate(0),calculate(1)

    def route_edges(self,path):
        endpoints={path['source'],path['target']}
        ring=endpoints&set((self.region or {}).get('node_ids',[]))
        fans=endpoints&self.fanout_nodes if self.auxiliary else set()
        key=(frozenset(ring),frozenset(fans))
        if key in self._edge_masks:return self._edge_masks[key]
        x=np.asarray(self.x);y=np.asarray(self.y)
        primary=lambda values:np.isclose(values,np.rint(values),atol=EPS,rtol=0)
        h=np.broadcast_to(primary(y)[:,None],(self.ny,self.nx-1)).copy()
        v=np.broadcast_to(primary(x)[None,:],(self.ny-1,self.nx)).copy()
        for node in fans:
            cx,cy=self.xy[node]
            h |= (np.abs(y[:,None]-cy)<self.fanout_extents[node]+EPS)&(np.maximum(np.abs(x[:-1][None,:]-cx),np.abs(x[1:][None,:]-cx))<self.fanout_extents[node]+EPS)
            v |= (np.abs(x[None,:]-cx)<self.fanout_extents[node]+EPS)&(np.maximum(np.abs(y[:-1][:,None]-cy),np.abs(y[1:][:,None]-cy))<self.fanout_extents[node]+EPS)
        clear_h=self.horizontal_clearance>=1-EPS;clear_v=self.vertical_clearance>=1-EPS
        for node in ring:
            if node not in self._attachment_masks:
                ah=np.zeros_like(clear_h);av=np.zeros_like(clear_v)
                dummy=dict(source=node,target=node)
                for iy,ix in np.argwhere(~clear_h):
                    a=(self.x[ix],self.y[iy]);b=(self.x[ix+1],self.y[iy])
                    ah[iy,ix]=incident_exempt(a,b,dummy,self.xy,self.region)
                for iy,ix in np.argwhere(~clear_v):
                    a=(self.x[ix],self.y[iy]);b=(self.x[ix],self.y[iy+1])
                    av[iy,ix]=incident_exempt(a,b,dummy,self.xy,self.region)
                self._attachment_masks[node]=(ah,av)
            ah,av=self._attachment_masks[node];clear_h=clear_h|ah;clear_v=clear_v|av
        result=(h&clear_h,v&clear_v);self._edge_masks[key]=result
        return result

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
                if abs(step)==1:self.hard_h[v//self.nx,min(v%self.nx,w%self.nx)]=True
                else:self.hard_v[min(v//self.nx,w//self.nx),v%self.nx]=True
                if p['priority']==2:
                    self.red_vertices[v].append((p,1 if abs(step)==1 else 2));self.red_vertices[w].append((p,1 if abs(step)==1 else 2))
    def index(self,p):return min(range(len(self.y)),key=lambda i:abs(self.y[i]-p[1]))*self.nx+min(range(len(self.x)),key=lambda i:abs(self.x[i]-p[0]))
    def point(self,v):return (self.x[v%self.nx],self.y[v//self.nx])
