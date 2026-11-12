const $ = id => document.getElementById(id);
const selected = {metabolites:new Set(config.initial.metabolites), reactions:new Set(config.initial.reactions)};
const lookup = Object.fromEntries(Object.entries(config.catalog).map(([kind, items]) => [kind, new Map(items.map(x => [x.id,x]))]));
let elements = config.elements, forced = new Set(), blue = new Set();
const originalDisplay = new Map();
const panel = $('highlight-panel');
document.querySelector('header')?.append(panel);
if(config.escher) $('save')?.addEventListener('click',()=>{
  if(typeof mapData!=='undefined') (mapData[0].metadata??={}).highlights=selectionJSON();
},true);
if(config.escher && $('map')) $('map').style.height='75vh';
function selectionJSON(){return Object.fromEntries(Object.entries(selected).map(([k,v])=>[k,[...v].sort()]));}
function download(name, value, type){
  const url=URL.createObjectURL(new Blob([value],{type}));
  const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
}
function isVisible(item){return item.ids.some(id=>{const e=$(id);return e && getComputedStyle(e).display!=='none';});}
function status(kind,id){
  const item=lookup[kind].get(id), matches=elements.filter(e=>e.kind===kind && e.key===id);
  if(!matches.length) return 'Not in saved map ('+(item?.reason || 'not_in_saved_map')+')';
  return matches.some(isVisible)?'Visible':'Hidden by filter';
}
function button(label,handler){const b=document.createElement('button');b.type='button';b.textContent=label;b.onclick=handler;return b;}
function results(){
  const kind=$('hl-kind').value, q=$('hl-search').value.trim().toLowerCase(), out=$('hl-results');out.replaceChildren();
  const matches=config.catalog[kind].filter(x=>(x.id+' '+x.name).toLowerCase().includes(q));
  for(const x of matches.slice(0,40)){
    out.append(button(`${x.id} — ${x.name}${x.compartment?' ['+x.compartment+']':''} — ${status(kind,x.id)}`,()=>{
      selected[kind].add(x.id);apply();
    }));
  }
  if(matches.length>40){const hint=document.createElement('span');hint.textContent=`${matches.length} matches; refine your search`;out.append(hint);}
}
function list(){
  const out=$('hl-selected');out.replaceChildren();
  for(const kind of ['metabolites','reactions']) for(const id of selected[kind]){
    out.append(button(`${kind==='metabolites'?'Compound':'Reaction'}: ${id} × — ${status(kind,id)}`,()=>{selected[kind].delete(id);apply();}));
  }
}
function paint(){
  for(const el of blue) el.classList.remove('user-blue');blue.clear();
  for(const id of forced){const el=$(id);if(el)el.style.display=originalDisplay.get(id) || '';}
  forced.clear();originalDisplay.clear();
  if(!config.escher && typeof update==='function') update();
  const force=$('hl-force').checked;
  function show(id){const el=$(id);if(!el)return;if(!originalDisplay.has(id))originalDisplay.set(id,el.style.display);el.style.display='';forced.add(id);}
  for(const item of elements){
    if(!selected[item.kind].has(item.key))continue;
    for(const id of item.ids){const el=$(id);if(!el)continue;el.classList.add('user-blue');el.dataset.hlRole=item.role;blue.add(el);if(force)show(id);}
    if(force)for(const id of item.endpoints||[])show(id);
  }
  if(config.escher)window.syncEscherCoverage?.();
}
function apply(){paint();results();list();}
function validate(value){
  if(!value || typeof value!=='object' || Array.isArray(value) || Object.keys(value).some(k=>!lookup[k]))throw Error('Invalid highlight JSON');
  const next={metabolites:[],reactions:[]};
  for(const kind of Object.keys(next)){
    const ids=value[kind]??[];
    if(!Array.isArray(ids)||ids.some(id=>typeof id!=='string'))throw Error('Expected lists of IDs');
    const unknown=ids.filter(id=>!lookup[kind].has(id));
    if(unknown.length)throw Error('Unknown IDs: '+unknown.join(', '));
    next[kind]=ids;
  }return next;
}
$('hl-kind').onchange=results;$('hl-search').oninput=results;$('hl-force').onchange=apply;
$('hl-add').onclick=()=>{
  const kind=$('hl-kind').value,ids=$('hl-bulk').value.split(/[\s,]+/).filter(Boolean),unknown=ids.filter(id=>!lookup[kind].has(id));
  for(const id of ids)if(lookup[kind].has(id))selected[kind].add(id);
  $('hl-message').textContent=unknown.length?'Unknown IDs: '+unknown.join(', '):'';apply();
};
$('hl-clear').onclick=()=>{for(const s of Object.values(selected))s.clear();$('hl-message').textContent='';apply();};
$('hl-save').onclick=()=>download('highlights.json',JSON.stringify(selectionJSON(),null,2),'application/json');
$('hl-load').onchange=async event=>{
  try{const file=event.target.files[0];if(!file)return;const next=validate(JSON.parse(await file.text()));
    for(const kind of Object.keys(next))selected[kind]=new Set(next[kind]);$('hl-message').textContent='';apply();
  }catch(e){$('hl-message').textContent=e.message;}event.target.value='';
};
$('hl-svg').onclick=()=>{
  const source=document.querySelector(config.escher?'#map svg':'#viewer svg');if(!source)return;
  const clone=source.cloneNode(true), a=[source,...source.querySelectorAll('*')],b=[clone,...clone.querySelectorAll('*')];
  const props=['fill','stroke','stroke-width','stroke-opacity','fill-opacity','opacity','display','visibility','font-family','font-size','font-weight','text-anchor','paint-order','stroke-linejoin','vector-effect'];
  a.forEach((el,i)=>{const css=getComputedStyle(el);for(const prop of props)b[i].style.setProperty(prop,css.getPropertyValue(prop));});
  clone.setAttribute('xmlns','http://www.w3.org/2000/svg');
  download('highlighted_map.svg',new XMLSerializer().serializeToString(clone),'image/svg+xml');
};
if(!config.escher){
  for(const id of ['coverage','coenzymes','focus']) $(id)?.addEventListener('input',apply);
  for(const id of ['coenzymes','focus']) $(id)?.addEventListener('change',apply);
  apply();
}else{
  if(!$('map')){apply();return;}
  let scheduled=false;
  function sync(){
    const map=window.builder?.map;if(!map)return;
    elements=[];
    for(const [id,node] of Object.entries(map.nodes))if(node.node_type==='metabolite'){
      elements.push({kind:'metabolites',key:node.bigg_id,ids:['n'+id],role:'node'});
      if(!lookup.metabolites.has(node.bigg_id)){const x={id:node.bigg_id,name:node.name||''};lookup.metabolites.set(x.id,x);config.catalog.metabolites.push(x);}
    }
    for(const [id,r] of Object.entries(map.reactions)){
      elements.push({kind:'reactions',key:r.bigg_id,ids:['r'+id],role:'edge',endpoints:[...new Set(Object.values(r.segments).flatMap(s=>['n'+s.from_node_id,'n'+s.to_node_id]))]});
      if(!lookup.reactions.has(r.bigg_id)){const x={id:r.bigg_id,name:r.name||''};lookup.reactions.set(x.id,x);config.catalog.reactions.push(x);}
    }
    apply();
  }
  document.addEventListener('escher-coverage-change',apply);
  const observer=new MutationObserver(records=>{
    if(!records.some(r=>r.type==='childList'))return;
    if(!scheduled){scheduled=true;requestAnimationFrame(()=>{scheduled=false;sync();});}
  });observer.observe($('map'),{childList:true,subtree:true});sync();
  // The viewer may finish loading before the controls are appended.
  window.addEventListener('load',sync,{once:true});
}
