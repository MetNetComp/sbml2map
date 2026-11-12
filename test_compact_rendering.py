"""Check actual label extents and terminal-only exchange compaction."""
import json
import unittest

from adaptive_routing import compact_terminal_corridors
from test_adaptive_routing import path


class CompactRenderingTests(unittest.TestCase):
    def test_label_packing_has_serializable_zero_overlap_report(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from escher_grid import separate_static_labels
        fig,ax=plt.subplots()
        try:
            for i in range(4):
                text=ax.annotate('overlapping label '+str(i),(.5,.5),xytext=(4,5),textcoords='offset points')
                text.set_gid('label-'+str(i))
            fig.canvas.draw()
            report=separate_static_labels(ax,fig.canvas.get_renderer())
            self.assertEqual(report['overlaps'],0)
            self.assertGreater(report['moved'],0)
            json.dumps(report)
        finally:plt.close(fig)

    def test_exchange_tail_compaction_preserves_internal_nodes(self):
        nodes={k:dict(kind='reaction' if k.startswith('rxn:') else 'metabolite')
               for k in ('rxn:EX_e','met:e','rxn:T','met:c','met:other')}
        xy={'rxn:EX_e':[20,4],'met:e':[19,4],'rxn:T':[18,4],'met:c':[4,4],'met:other':[2,2]}
        paths=[path('EX_e','rxn:EX_e','met:e',xy['rxn:EX_e'],xy['met:e'],1),
               path('T','rxn:T','met:e',xy['rxn:T'],xy['met:e'],1),
               path('T','rxn:T','met:c',xy['rxn:T'],xy['met:c'],1)]
        for p in paths:p['mid']=p['target'][4:];p['straight_locked']=True
        group=dict(reaction_id='EX_e',exchange='rxn:EX_e',external='met:e',internal='met:c',
                   primary_transport='T',face='right',node_ids=['rxn:EX_e','met:e','rxn:T','met:c'],corridor=[[20,4],[4,4]])
        backbones=dict(routing_boundary=[0,20,0,10],fixed_node_offsets={'met:c':[4,4]},adaptive_layout={},
                       exchange_corridors=dict(groups=[group]),
                       crossing_routing=dict(routes=paths,signature='initial',report=dict(conditional_crossings=[])))
        self.assertTrue(compact_terminal_corridors(nodes,xy,backbones,lambda *a,**k:None))
        self.assertEqual(xy['met:c'],[4,4]);self.assertEqual(xy['met:other'],[2,2])
        self.assertEqual(backbones['routing_boundary'][1],7)
        self.assertTrue(backbones['crossing_routing']['report']['verified'])

    def test_moving_endpoint_along_a_short_first_leg_keeps_orthogonal_elbow(self):
        ids=('rxn:EX_e','met:e','rxn:T','met:c','met:other','met:right')
        nodes={k:dict(kind='reaction' if k.startswith('rxn:') else 'metabolite') for k in ids}
        xy={'rxn:EX_e':[20,4],'met:e':[19,4],'rxn:T':[18,4],'met:c':[4,4],'met:other':[2,2],'met:right':[18,8]}
        paths=[path('EX_e','rxn:EX_e','met:e',xy['rxn:EX_e'],xy['met:e'],1),
               path('T','rxn:T','met:e',xy['rxn:T'],xy['met:e'],1),
               path('T','rxn:T','met:c',xy['rxn:T'],xy['met:c'],1)]
        for p in paths:p['mid']=p['target'][4:];p['straight_locked']=True
        branch=path('T','rxn:T','met:other',xy['rxn:T'],xy['met:other'],0)
        branch.update(mid='other',points=[(18,4),(17.75,4),(17.75,2),(2,2)],auxiliary_lane=True)
        paths.append(branch)
        group=dict(reaction_id='EX_e',exchange='rxn:EX_e',external='met:e',internal='met:c',
                   primary_transport='T',face='right',node_ids=['rxn:EX_e','met:e','rxn:T','met:c'],corridor=[[20,4],[4,4]])
        b=dict(routing_boundary=[0,20,0,10],fixed_node_offsets={'met:c':[4,4]},adaptive_layout={},
               exchange_corridors=dict(groups=[group]),
               crossing_routing=dict(routes=paths,signature='initial',report=dict(conditional_crossings=[])))
        self.assertTrue(compact_terminal_corridors(nodes,xy,b,lambda *a,**k:None))
        branch=next(p for p in b['crossing_routing']['routes'] if p['mid']=='other')
        self.assertEqual(branch['points'][0],(17.,4.))
        self.assertEqual(branch['points'][1],(17.75,4))
        self.assertTrue(all(a[0]==c[0] or a[1]==c[1] for a,c in zip(branch['points'],branch['points'][1:])))


if __name__=='__main__':unittest.main()
