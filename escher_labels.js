// Resolve metabolite/reaction labels together using their screen-space bounds.
(()=>{
const host=document.getElementById('map');
const options=mapData[0].metadata.label_display;
const style=document.createElement('style');
style.textContent=`[data-escher-label-hidden],[data-escher-label-hidden] *{visibility:hidden!important;}
#map .reaction-label,#map .reaction-label text,#map .node-label,#map .node-label text{
 paint-order:stroke fill;stroke:#fff!important;stroke-width:3px;stroke-linejoin:round;vector-effect:non-scaling-stroke;
}`;
document.head.append(style);
let pending=false;
const controls=[['show-metabolite-labels','show_metabolites'],['show-reaction-labels','show_reactions'],['allow-label-overlap','allow_overlap']];
for(const [id,key] of controls){
 const input=document.getElementById(id);if(!input)continue;input.checked=options[key];
 input.addEventListener('change',()=>{options[key]=input.checked;schedule();});
}
function schedule(){if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;sync();});}}
function sync(){
 const map=window.builder?.map;if(!map)return;
 const labels=[];
 for(const [selector,kind,enabled,priorities] of [
  ['#reactions .reaction-label','reaction',options.show_reactions,options.reaction_priorities],
  ['#nodes .node-label','metabolite',options.show_metabolites,options.metabolite_priorities]
 ]){
  for(const el of host.querySelectorAll(selector)){
   const owner=el.closest(kind==='reaction'?'.reaction':'.node');
   if(!owner)continue;
   const id=owner.id.slice(1),record=(kind==='reaction'?map.reactions:map.nodes)[id];
   labels.push({el,key:owner.id,rank:priorities[id]??(kind==='metabolite'&&!record?.node_is_primary?1:2),
                order:kind+':'+(record?.bigg_id??id)+':'+id,enabled,hidden:!enabled});
  }
 }
 // Read all bounds before applying suppression. Previously hidden labels must be
 // measured again so zooming, edits, and control changes can reveal them.
 for(const item of labels)item.el.removeAttribute('data-escher-label-hidden');
 for(const item of labels){
  if(!item.enabled)continue;
  const css=getComputedStyle(item.el),r=item.el.getBoundingClientRect();
  if(css.display==='none'||css.visibility==='hidden'||!r.width||!r.height)continue;
  item.box={left:r.left-2,right:r.right+2,top:r.top-2,bottom:r.bottom+2};
 }
 labels.sort((a,b)=>b.rank-a.rank||(a.order<b.order?-1:a.order>b.order?1:0));
 const cells=new Map(),cellSize=100;
 function keys(box){
  const result=[];
  for(let x=Math.floor(box.left/cellSize);x<=Math.floor(box.right/cellSize);x++)
   for(let y=Math.floor(box.top/cellSize);y<=Math.floor(box.bottom/cellSize);y++)result.push(x+','+y);
  return result;
 }
 for(const item of labels){
  if(!item.box||options.allow_overlap)continue;
  const b=item.box,slots=keys(b);
  item.hidden=slots.some(key=>(cells.get(key)||[]).some(a=>a.left<b.right&&a.right>b.left&&a.top<b.bottom&&a.bottom>b.top));
  if(!item.hidden)for(const key of slots){if(!cells.has(key))cells.set(key,[]);cells.get(key).push(b);}
 }
 const byKey=new Map();
 for(const item of labels){item.el.toggleAttribute('data-escher-label-hidden',item.hidden);byKey.set(item.key,item.el);}
 // The priority layer contains non-interactive copies of reaction labels.
 for(const clone of host.querySelectorAll('#escher-priority-layer [data-source-label]')){
  const source=byKey.get(clone.dataset.sourceLabel);
  const hidden=!source||getComputedStyle(source).visibility==='hidden'||getComputedStyle(source).display==='none';
  clone.toggleAttribute('data-escher-label-hidden',hidden);
  const originals=source?[source,...source.querySelectorAll('*')]:[];
  [clone,...clone.querySelectorAll('*')].forEach((el,i)=>{
   const visibility=hidden?'hidden':getComputedStyle(originals[i]||source).visibility;
   if(el.style.visibility!==visibility)el.style.visibility=visibility;
  });
 }
}
window.syncEscherLabels=sync;
window.syncEscherLayers?.(true);
const observer=new MutationObserver(records=>{
 if(records.some(r=>!r.target.closest?.('#escher-priority-layer')))schedule();
});
observer.observe(host,{subtree:true,childList:true,characterData:true,attributes:true,
 attributeFilter:['d','class','style','transform','visibility','x','y','dx','dy','font-size']});
new ResizeObserver(schedule).observe(host);
window.addEventListener('resize',schedule);
document.fonts?.ready.then(schedule);
schedule();
})();
