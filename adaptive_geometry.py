"""Compact TCA anchors at exact intersections with primary grid lines."""
import math
import numpy as np


def compact_tca(fixed, backbones, scale=2/3, grid_point_keys=()):
    region=backbones.get('tca_geometry')
    if not region or region['shape']!='circle':return fixed,backbones
    old_radius=region['radius'];radius=math.sqrt(round((old_radius*scale)**2))
    center=[float(math.ceil(region['center_offset'][0])),float(math.floor(region['center_offset'][1]))]
    candidates=set()
    for axis in (0,1):
        for level in range(math.ceil(center[axis]-radius),math.floor(center[axis]+radius)+1):
            delta=math.sqrt(max(0.,radius*radius-(level-center[axis])**2))
            for sign in (-1,1):
                p=list(center);p[axis]=float(level);p[1-axis]+=sign*delta
                candidates.add(tuple(p))
    phase=3*math.pi/4
    angle=lambda p:(phase-math.atan2(p[1]-center[1],p[0]-center[0]))%(2*math.pi)
    points=sorted(candidates,key=angle);angles=[angle(p) for p in points]
    old_center=region['center_offset']
    keys=sorted(region['node_ids'],key=lambda k:(phase-math.atan2(fixed[k][1]-old_center[1],fixed[k][0]-old_center[0]))%(2*math.pi))
    n,m=len(keys),len(points)
    if m<n:raise ValueError('Too few circle/grid-line intersections for TCA nodes')
    # Ordered matching avoids crossing arcs and places anchors approximately evenly.
    dp=np.full((n+1,m+1),np.inf);dp[0,:]=0.;take=np.zeros((n+1,m+1),dtype=bool)
    for i in range(1,n+1):
        target=(i-.5)*2*math.pi/n
        for j in range(1,m+1):
            dp[i,j]=dp[i,j-1]
            candidate=points[j-1]
            permitted=keys[i-1] not in grid_point_keys or all(abs(v-round(v))<1e-7 for v in candidate)
            value=dp[i-1,j-1]+(angles[j-1]-target)**2 if permitted else math.inf
            if value<dp[i,j]:dp[i,j]=value;take[i,j]=True
    selected=[];i,j=n,m
    if not math.isfinite(dp[n,m]):raise ValueError('Required TCA grid-point ports cannot fit circle')
    while i:
        if take[i,j]:selected.append(points[j-1]);i-=1
        j-=1
    selected.reverse()
    positions=dict(zip(keys,map(list,selected)))
    fixed.update(positions);backbones['fixed_node_offsets'].update(positions)
    for route in backbones['routes']:
        for key in positions:
            if key in route.get('node_offsets',{}):route['node_offsets'][key]=positions[key]
    region.update(center_offset=center,radius=radius,previous_radius=old_radius,diameter_ratio=radius/old_radius,requested_diameter_ratio=scale,
                  anchor_mode='grid_line_intersection',all_tca_nodes_on_grid=False,
                  all_tca_nodes_on_grid_lines=True,node_grid_step=1.,lattice_point_count=m)
    backbones['circle'].update(center_offset=center,radius=radius)
    return fixed,backbones
