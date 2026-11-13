"""Independently validate a saved sequential exchange drawing."""
import argparse
import csv
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from adaptive_routing import validate_paths, validate_node_grid
from crossing_routing import validate_saved_geometry
from escher_grid import model_hash
from export_escher import validate_map
from model_input import load_model
from reaction_network import build_network, validate_network
from red_objective import route_length
from sequential_orange import organic_graph, hop_distances


def validate(folder, source):
    folder=Path(folder)
    read=lambda name:json.loads((folder/name).read_text())
    report=read('report.json');backbones=read('major_backbones.json')
    plan=read('sequential_orange.json');model=read('display_model.json');data=read('escher_map.json')
    original=load_model(source)
    assert read('source_model.json')==original
    assert report['source_model_sha256']==model_hash(original)
    rows=list(csv.DictReader((folder/'coordinates.csv').open()))
    xy={r['id']:[float(r['x']),float(r['y'])] for r in rows}
    nodes,_=build_network(model,{r['original_id'] for r in rows if r['kind']=='currency'})
    edges=list(csv.DictReader((folder/'edges.csv').open()))
    for edge in edges:edge['coefficient']=float(edge['coefficient'])
    validate_network(model,nodes,edges,xy)
    assert validate_map(data,model,nodes,xy)['schema_valid']
    actual=validate_saved_geometry(data)
    assert actual['verified'] and actual['red_preservation']['exactly_preserved']
    paths=backbones['crossing_routing']['routes'];by_key={(p['rid'],p['mid']):p for p in paths}
    validate_paths(paths,xy,backbones['routing_boundary'],backbones['tca_exclusion']['region'],[])
    validate_node_grid(nodes,xy,backbones['tca_exclusion']['region'])
    red=[p for p in paths if p['priority']==2]
    roots={k for p in red for k in (p['source'],p['target'])}
    adjacency,lookup=organic_graph(nodes,edges)
    distances=hop_distances(adjacency,roots)
    connected=set(roots);chosen=set();new_length=0.;order=[];orange_attachments=0
    for record in plan['records']:
        terminal=record['exchange'];order.append((distances[terminal],terminal))
        assert record['red_hops']==distances[terminal]
        if record.get('already_connected'):
            assert terminal in connected;continue
        chain=record['path']
        assert chain[0]==terminal and chain[-1] in connected
        assert all(k not in connected for k in chain[:-1])
        assert all(nodes[k]['kind']!='currency' and not nodes[k].get('shared_cofactor') for k in chain)
        assert record['attachment_node']==chain[-1]
        expected='red' if chain[-1] in roots else 'orange'
        assert record['attachment_color']==expected
        orange_attachments+=expected=='orange'
        length=0.
        for a,b in zip(chain,chain[1:]):
            index=lookup[frozenset((a,b))];edge=edges[index]
            assert index not in chosen
            path=by_key[edge['reaction_id'],edge['metabolite_id']]
            assert path['priority']==1 and path['sequential_fixed']
            length+=route_length(path);chosen.add(index)
        assert math.isclose(length,record['new_length_grid'],abs_tol=1e-6)
        assert record['solver']['lower_bound']<=length+1e-6
        connected.update(chain);new_length+=length
    assert order==sorted(order)
    assert chosen==set(plan['edge_indices'])
    assert math.isclose(new_length,sum(route_length(p) for p in paths if p['priority']==1),abs_tol=1e-6)
    render=read('render_validation.json');assert render['labels']['overlaps']==0
    # Fit the renderer's affine transform from endpoints, then check every SVG
    # line/control point against its saved route (including reverse reactions).
    svg=ET.parse(folder/'layout.svg');namespace='{http://www.w3.org/2000/svg}'
    expected_paths=[];actual_paths=[]
    for i,edge in enumerate(edges):
        path=by_key[edge['reaction_id'],edge['metabolite_id']]
        points=[path['points'][0],*path['curve_controls'],path['points'][-1]] if path.get('curve_controls') else path['points']
        start=xy[edge['target'] if edge['direction']=='reverse' else edge['source']]
        if math.dist(points[0],start)>1e-7:points=list(reversed(points))
        element=svg.find(f'.//*[@id="edge-{i}"]/{namespace}path')
        assert element is not None
        numbers=re.findall(r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?',element.get('d'))
        vertices=np.array(list(map(float,numbers))).reshape(-1,2)
        assert len(vertices)==len(points),(i,len(vertices),len(points))
        expected_paths.append(np.array(points));actual_paths.append(vertices)
    x=np.concatenate(expected_paths);y=np.concatenate(actual_paths)
    transform=np.linalg.lstsq(np.c_[x,np.ones(len(x))],y,rcond=None)[0]
    assert np.max(np.abs(np.c_[x,np.ones(len(x))]@transform-y))<2e-5
    return dict(model='e_coli_core',shape=backbones['tca_exclusion']['region']['shape'],
                drawn_reactions=report['drawn_reactions'],source_model_unchanged=True,
                grid=validate_node_grid(nodes,xy,backbones['tca_exclusion']['region']),
                red_preservation=actual['red_preservation'],parallel_conflicts=0,labels=render['labels'],
                exchanges=len(plan['records']),orange_attachments=orange_attachments,orange_connections=len(chosen),
                orange_length_grid=new_length,coenzyme_transit=False,initial_red_hop_order_verified=True,
                sequential_orange_preserved=True,svg_escher_verified=True,
                boundary=backbones['routing_boundary'],
                global_optimality_proven=False,topology_search=plan['topology_search'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path);parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=validate(args.folder,args.model)
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    print(text)
