// Paint important segments above other reactions without changing Escher's DOM
// bindings. Non-interactive copies have no IDs or Escher segment/node classes.
(()=>{
const NS='http://www.w3.org/2000/svg';
let pending=false,last='';
const props=['stroke','fill','stroke-width','stroke-opacity','fill-opacity','opacity','stroke-dasharray','stroke-linecap','stroke-linejoin','display','visibility','font-family','font-size','font-weight','text-anchor'];
function sync(){
 const map=window.builder?.map,root=document.getElementById('reactions');
 if(!map||!root)return;
 const fixed=document.getElementById('major-pathway-colors');
 const structural=!map.has_data_on_reactions && !!fixed?.sheet && !fixed.sheet.disabled;
 const parent=root.parentNode,items=[];
 for(const path of root.querySelectorAll('.segment')){
  const source=path.parentElement,sid=source.id.replace(/^s/,'');
  const rank=path.closest('.user-blue')?3:structural?(mapData[0].metadata.segment_priorities?.[sid]||0):0;
  if(!rank || !path.getClientRects().length)continue;
  const css=getComputedStyle(path),matrix=parent.getCTM().inverse().multiply(source.getCTM());
  items.push({source,rank,key:[sid,rank,path.getAttribute('d'),css.stroke,css.strokeWidth,source.innerHTML,matrix.toString()].join('|')});
 }
 const labels=[...root.querySelectorAll('.reaction-label')].filter(label=>label.getClientRects().length);
 const key=items.map(x=>x.key).join(';')+'::'+labels.map(label=>[label.closest('.reaction').id,label.textContent,label.getCTM()?.toString()].join('|')).join(';');
 if(key===last && (document.getElementById('escher-priority-layer') || !items.length))return;
 last=key;
 document.getElementById('escher-priority-layer')?.remove();
 if(!items.length)return;
 const layer=document.createElementNS(NS,'g');layer.id='escher-priority-layer';layer.style.pointerEvents='none';
 parent.insertBefore(layer,root.nextSibling);
 function copy(source){
  const clone=source.cloneNode(true),a=[source,...source.querySelectorAll('*')],b=[clone,...clone.querySelectorAll('*')];
  a.forEach((el,i)=>{const style=getComputedStyle(el);b[i].removeAttribute('id');b[i].removeAttribute('class');for(const p of props)b[i].style.setProperty(p,style.getPropertyValue(p));b[i].style.pointerEvents='none';});
  clone.setAttribute('transform',parent.getCTM().inverse().multiply(source.getCTM()).toString());return clone;
 }
 items.sort((a,b)=>a.rank-b.rank);
 for(const item of items){const clone=copy(item.source);clone.dataset.sourceSegment=item.source.id;clone.dataset.priority=item.rank;layer.append(clone);}
 // Keep reaction labels and biological node symbols above the emphasized lines.
 for(const label of labels){const clone=copy(label);clone.dataset.sourceLabel=label.closest('.reaction').id;layer.append(clone);}
 for(const id of ['nodes','text-labels']){const el=document.getElementById(id);if(el?.parentNode===parent)parent.append(el);}
 window.syncEscherLabels?.();
}
window.syncEscherLayers=sync;
const observer=new MutationObserver(records=>{
 if(records.every(r=>r.target.closest?.('#escher-priority-layer')))return;
 if(!pending){pending=true;requestAnimationFrame(()=>{pending=false;sync();});}
});
observer.observe(document.getElementById('map'),{subtree:true,childList:true,attributes:true,attributeFilter:['d','class','style','transform']});
sync();
})();
