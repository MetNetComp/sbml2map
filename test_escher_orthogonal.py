"""Orthogonal Escher routes preserve source chemistry and exported endpoints."""
import copy
import json
from pathlib import Path
import unittest

from export_escher import create_map,validate_map
from reaction_network import build_network
from run_fba_escher import read_fba_layout


class OrthogonalExportTests(unittest.TestCase):
    def fixture(self):
        model=dict(id='routing',metabolites=[dict(id=m) for m in ('a','b','atp')],reactions=[
            dict(id='reverse',metabolites={'a':-1,'b':1,'atp':-1},lower_bound=-10,upper_bound=0),
            dict(id='reversible',metabolites={'b':-1,'a':1},lower_bound=-10,upper_bound=10)])
        nodes,edges=build_network(model,{'atp'})
        xy={key:[float(i),float(i*i)] for i,key in enumerate(nodes)}
        return model,nodes,edges,xy

    def test_every_non_aligned_incidence_gets_one_elbow(self):
        model,nodes,edges,xy=self.fixture();original=copy.deepcopy((model,nodes,xy))
        data=create_map(model,nodes,edges,xy)
        check=validate_map(data,model,nodes,xy)
        expected=sum(xy[e['source']][0]!=xy[e['target']][0] and xy[e['source']][1]!=xy[e['target']][1] for e in edges)
        self.assertEqual(check['routing_nodes'],expected)
        self.assertEqual(check['segments'],len(edges)+expected)
        self.assertEqual((model,nodes,xy),original)
        self.assertTrue(check['orthogonal_segments_verified'])
        self.assertEqual(check['reverse_only_reactions'],['reverse'])

    def test_native_save_sidecar_keeps_routes_verifiable(self):
        model,nodes,edges,xy=self.fixture();data=create_map(model,nodes,edges,xy)
        # The viewer fills null controls with collinear midpoint controls.
        for reaction in data[1]['reactions'].values():
            for segment in reaction['segments'].values():
                a,b=(data[1]['nodes'][segment[k]] for k in ('from_node_id','to_node_id'))
                segment['b1']={axis:(a[axis]+b[axis])/2 for axis in ('x','y')}
                segment['b2']=dict(segment['b1'])
        sidecar=data[0].pop('metadata')
        check=validate_map(data,model,nodes,xy,source_metadata=sidecar)
        self.assertTrue(check['orthogonal_segments_verified'])
        self.assertGreater(check['routing_nodes'],0)

    def test_moved_corner_and_wrong_source_mapping_are_rejected(self):
        model,nodes,edges,xy=self.fixture();data=create_map(model,nodes,edges,xy)
        broken=copy.deepcopy(data);key=broken[0]['metadata']['routing_node_ids'][0]
        broken[1]['nodes'][key]['x']+=1
        with self.assertRaisesRegex(ValueError,'horizontal or vertical'):validate_map(broken,model,nodes,xy)
        broken=copy.deepcopy(data);sid=next(iter(broken[0]['metadata']['segment_source_edges']))
        broken[0]['metadata']['segment_source_edges'][sid]=['wrong','nodes']
        with self.assertRaisesRegex(ValueError,'incidence'):validate_map(broken,model,nodes,xy)

    def test_saved_e_coli_core_map_and_fba_are_orthogonal(self):
        root=Path(__file__).resolve().parent;folder=root/'results/e_coli_core'
        if not folder.is_dir():folder=root/'public/results/e_coli_core'
        model=json.loads((folder/'source_model.json').read_text())
        projected,nodes,edges,xy=read_fba_layout(folder,model)
        data=json.loads((folder/'escher_map.json').read_text())
        check=validate_map(data,projected,nodes,xy)
        self.assertEqual(check['source_nodes'],len(nodes))
        self.assertTrue(check['orthogonal_segments_verified'])
        self.assertGreater(check['routing_nodes'],0)
        report=json.loads((folder/'fba_report.json').read_text())
        self.assertTrue(report['map_validation']['orthogonal_segments_verified'])
        self.assertAlmostEqual(report['objective_value'],0.8739215069684304,places=8)


if __name__=='__main__':unittest.main()
