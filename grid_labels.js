// Label controls share SVG coverage/cofactor filtering without moving geometry.
(()=>{
const host=document.getElementById('viewer'),options=window.gridLabelOptions;
const items=registry.filter(item=>item.label_kind).map(item=>({...item,el:document.getElementById(item.id)})).filter(item=>item.el);
for(const item of items){item.el.dataset.labelKind=item.label_kind;item.el.dataset.labelCofactor=String(item.coenzyme);}
const style=document.createElement('style');
style.textContent='[data-grid-label-hidden],[data-grid-label-hidden] *{visibility:hidden!important;}';
document.head.append(style);
let pending=false;
function schedule(){if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;sync();});}}
for(const [id,key] of [['show-metabolite-labels','show_metabolites'],['show-reaction-labels','show_reactions'],['allow-label-overlap','allow_overlap']]){
 const input=document.getElementById(id);input.checked=options[key];
 input.addEventListener('change',()=>{options[key]=input.checked;schedule();});
}
function sync(){
 for(const item of items){item.el.removeAttribute('data-grid-label-hidden');item.box=null;item.hidden=!(item.label_kind==='reaction'?options.show_reactions:options.show_metabolites);}
 for(const item of items){
  if(item.hidden)continue;
  const css=getComputedStyle(item.el),r=item.el.getBoundingClientRect();
  if(css.display==='none'||css.visibility==='hidden'||!r.width||!r.height)continue;
  item.box={left:r.left-1.5,right:r.right+1.5,top:r.top-1.5,bottom:r.bottom+1.5};
 }
 const cells=new Map(),size=100;
 function slots(b){const out=[];for(let x=Math.floor(b.left/size);x<=Math.floor(b.right/size);x++)for(let y=Math.floor(b.top/size);y<=Math.floor(b.bottom/size);y++)out.push(x+','+y);return out;}
 for(const item of [...items].sort((a,b)=>b.label_priority-a.label_priority||a.id.localeCompare(b.id))){
  if(item.box&&!options.allow_overlap){
   const b=item.box,keys=slots(b);
   item.hidden=keys.some(key=>(cells.get(key)||[]).some(a=>a.left<b.right&&a.right>b.left&&a.top<b.bottom&&a.bottom>b.top));
   if(!item.hidden)for(const key of keys){if(!cells.has(key))cells.set(key,[]);cells.get(key).push(b);}
  }
  item.el.toggleAttribute('data-grid-label-hidden',item.hidden);
 }
}
window.syncGridLabels=schedule;
new ResizeObserver(schedule).observe(host.querySelector('svg'));
window.addEventListener('resize',schedule);document.fonts?.ready.then(schedule);schedule();
})();
