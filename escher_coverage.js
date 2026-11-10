// Filter saved reactions by sampled rank without modifying the Escher map data.
(()=>{
const options=mapData[0].metadata.reaction_coverage;
if(!options)return;
const host=document.getElementById('map'),slider=document.getElementById('escher-coverage');
const style=document.createElement('style');style.textContent='[data-escher-coverage-hidden]{display:none !important;}';
document.head.append(style);
document.getElementById('escher-coverage-controls').hidden=false;
slider.value=options.initial_percent??100;
let pending=false;
function schedule(){if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;sync();});}}
function sync(){
 const map=window.builder?.map;if(!map)return;
 const percent=Math.max(0,Math.min(100,Number(slider.value))),order=options.reaction_order.filter(id=>map.reactions[id]);
 const visible=new Set(order.slice(0,Math.ceil(percent*order.length/100)));
 const nodes=new Set(),force=document.getElementById('hl-force')?.checked;
 // New reactions have no sampled rank and are visible at 100%.
 if(percent===100)for(const id of Object.keys(map.reactions))visible.add(id);
 if(force)for(const el of host.querySelectorAll('#reactions .reaction.user-blue'))visible.add(el.id.slice(1));
 for(const id of visible)for(const segment of Object.values(map.reactions[id]?.segments||{})){
  nodes.add(segment.from_node_id);nodes.add(segment.to_node_id);
 }
 if(force)for(const el of host.querySelectorAll('#nodes .node.user-blue'))nodes.add(el.id.slice(1));
 let changed=false;
 for(const [selector,set] of [['#reactions .reaction',visible],['#nodes .node',nodes]]){
  for(const el of host.querySelectorAll(selector)){
   const hidden=!set.has(el.id.slice(1));
   if(el.hasAttribute('data-escher-coverage-hidden')!==hidden){el.toggleAttribute('data-escher-coverage-hidden',hidden);changed=true;}
  }
 }
 options.initial_percent=percent;
 document.getElementById('escher-coverage-value').textContent=percent+'%';
 document.getElementById('escher-coverage-count').textContent=visible.size+' / '+Object.keys(map.reactions).length+' reactions';
 if(changed){window.syncEscherLayers?.();window.syncEscherLabels?.();}
}
window.syncEscherCoverage=sync;
slider.addEventListener('input',()=>{sync();document.dispatchEvent(new Event('escher-coverage-change'));});
new MutationObserver(records=>{
 if(records.some(r=>!r.target.closest?.('#escher-priority-layer')))schedule();
}).observe(host,{subtree:true,childList:true});
sync();
})();
