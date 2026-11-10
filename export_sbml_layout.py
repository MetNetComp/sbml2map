#!/usr/bin/env python3
"""Add standard SBML Level 3 Layout to the original model, on explicit request."""
import argparse
from collections import defaultdict
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import re
import tempfile
from urllib.parse import quote
from xml.dom import Node, minidom
import xml.etree.ElementTree as ET

from model_input import load_model

XMLNS = 'http://www.w3.org/2000/xmlns/'
LAYOUT = 'http://www.sbml.org/sbml/level3/version1/layout/version1'


def _library():
    try:
        import libsbml
    except ImportError as exc:
        raise RuntimeError('SBML Layout export requires: python -m pip install -r requirements-sbml-layout.txt') from exc
    return libsbml


def validate_export_request(model_path, layout_id='sbml2map_layout'):
    path=Path(model_path)
    if path.suffix.lower()=='.json':
        raise ValueError('SBML Layout export currently requires an SBML Level 3 input, not JSON')
    if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*',layout_id):
        raise ValueError('Layout ID must be a valid SBML SId')
    try:
        _,root=next(ET.iterparse(path,events=('start',)))
    except (ET.ParseError,StopIteration) as exc:
        raise ValueError('Invalid SBML input') from exc
    if root.tag.rsplit('}',1)[-1]!='sbml' or root.get('level')!='3' or root.get('version') not in ('1','2'):
        raise ValueError('SBML Layout export currently supports SBML Level 3 Version 1 or 2')
    _library()


def _decode(identifier,prefix):
    return re.sub(r'__(\d+)__',lambda m:chr(int(m[1])),identifier.removeprefix(prefix))


def _reverse_ids(elements,prefix):
    result={}
    for element in elements:
        key=_decode(element.getId(),prefix)
        if key in result:raise ValueError('Ambiguous original SBML ID: '+key)
        result[key]=element.getId()
    return result


def _snapshot(node):
    if node.nodeType==Node.ELEMENT_NODE:
        attrs=tuple(sorted((a.namespaceURI or '',a.localName or a.name,a.value) for a in node.attributes.values()))
        return (node.namespaceURI or '',node.localName or node.tagName,attrs,
                tuple(v for child in node.childNodes if (v:=_snapshot(child)) is not None))
    if node.nodeType in (Node.TEXT_NODE,Node.CDATA_SECTION_NODE):
        return ('text',node.data) if node.data.strip() else None
    if node.nodeType in (Node.COMMENT_NODE,Node.PROCESSING_INSTRUCTION_NODE):
        return (node.nodeType,node.toxml())
    return None


def _messages(document):
    messages=[]
    for i in range(document.getNumErrors()):
        error=document.getError(i)
        messages.append(dict(id=error.getErrorId(),severity=error.getSeverity(),message=error.getMessage().strip()))
    return messages


def _require_valid(document,libsbml,context):
    errors=[m for m in _messages(document) if m['severity']>=libsbml.LIBSBML_SEV_ERROR]
    if errors:raise ValueError(context+': '+'; '.join(m['message'] for m in errors[:3]))


def _read_saved_layout(directory,source):
    from export_escher import validate_map
    from reaction_network import build_network
    saved=directory/'source_model.json'
    if saved.exists() and json.loads(saved.read_text())!=source:
        raise ValueError('Input SBML does not match the saved source model')
    model_file=directory/('display_model.json' if (directory/'display_model.json').exists() else 'model.json')
    projected=json.loads(model_file.read_text())
    original={r['id']:r for r in source['reactions']}
    for reaction in projected['reactions']:
        parent=original.get(reaction['id'])
        if parent is None or any(reaction[k]!=parent[k] for k in ('lower_bound','upper_bound')) or any(
                parent['metabolites'].get(mid)!=value for mid,value in reaction['metabolites'].items()):
            raise ValueError('Saved display chemistry or bounds do not match the input SBML')
    with (directory/'coordinates.csv').open() as stream:rows=list(csv.DictReader(stream))
    nodes,edges=build_network(projected,{r['original_id'] for r in rows if r['kind']=='currency'})
    if len(rows)!=len(nodes) or {r['id'] for r in rows}!=set(nodes):
        raise ValueError('Saved coordinates do not match the display model')
    for row in rows:
        if any(row[k]!=nodes[row['id']][k] for k in ('kind','original_id')):
            raise ValueError('Saved coordinate identity mismatch')
    xy={r['id']:[float(r['x']),float(r['y'])] for r in rows}
    data=json.loads((directory/'escher_map.json').read_text())
    if nodes:
        validate_map(data,projected,nodes,xy)
    elif data[1]['nodes'] or data[1]['reactions']:
        raise ValueError('Unexpected nodes in an empty layout')
    return data


def _incidences(data):
    """Recover full marker-to-metabolite polylines, including lane offsets."""
    metadata=data[0]['metadata'];source=metadata['source_node_ids']
    for reaction in data[1]['reactions'].values():
        rid=reaction['bigg_id'];marker=metadata['reactions'][rid]['marker_node_id'];groups=defaultdict(list)
        for sid,segment in reaction['segments'].items():
            pair=metadata.get('segment_source_edges',{}).get(sid)
            if pair is None:pair=[source[segment['from_node_id']],source[segment['to_node_id']]]
            groups[frozenset(pair)].append(segment)
        for pair,segments in groups.items():
            graph=defaultdict(list)
            for segment in segments:
                a,b=segment['from_node_id'],segment['to_node_id'];graph[a].append(b);graph[b].append(a)
            target=next((n for n in graph if n in source and n!=marker),None)
            if target is None:raise ValueError('Missing metabolite endpoint')
            path=[marker];previous=None
            while path[-1]!=target:
                choices=[n for n in graph[path[-1]] if n!=previous]
                if len(choices)!=1 or choices[0] in path:raise ValueError('Invalid incidence path')
                previous=path[-1];path.append(choices[0])
            yield rid,target,path


def export_sbml_layout(model_path,layout_dir,output=None,layout_id='sbml2map_layout',verify_fba=False):
    model_path,layout_dir=Path(model_path),Path(layout_dir)
    output=Path(output) if output is not None else layout_dir/'layout.sbml'
    if output.resolve()==model_path.resolve():raise ValueError('Output must not overwrite the input SBML')
    validate_export_request(model_path,layout_id);lib=_library()
    raw=model_path.read_bytes();original_dom=minidom.parseString(raw)
    document=lib.readSBMLFromString(raw.decode('utf-8-sig'))
    _require_valid(document,lib,'Cannot read original SBML')
    document.checkConsistency();_require_valid(document,lib,'Original SBML consistency failure')
    original_warnings=_messages(document)
    source=load_model(model_path);data=_read_saved_layout(layout_dir,source)
    original_model=document.getModel()
    species_ids=_reverse_ids(original_model.getListOfSpecies(),'M_')
    reaction_ids=_reverse_ids(original_model.getListOfReactions(),'R_')
    original_reactions={r['id']:r for r in source['reactions']}
    used={a.value for e in original_dom.getElementsByTagName('*') for a in e.attributes.values() if a.localName=='id'}
    def identifier(preferred):
        preferred=re.sub(r'[^A-Za-z0-9_]','_',preferred)
        if not preferred or preferred[0].isdigit():preferred='g_'+preferred
        result=preferred;i=2
        while result in used:result=f'{preferred}_{i}';i+=1
        used.add(result);return result
    layout_id=identifier(layout_id)
    root=original_dom.documentElement
    prefix=next((a.localName for a in root.attributes.values() if a.namespaceURI==XMLNS and a.value==LAYOUT),None)
    if not prefix:
        prefix='layout';i=2
        while root.hasAttribute('xmlns:'+prefix):prefix=f'layout{i}';i+=1
    document.enablePackage(LAYOUT,prefix,True)
    document.setPackageRequired('layout',False)
    plugin=document.getModel().getPlugin('layout');existing_layouts=plugin.getNumLayouts()
    layout=plugin.createLayout();layout.setId(layout_id)
    layout.setName('sbml2map saved diagram (displayed subset)')
    body=data[1];metadata=data[0].get('metadata',{});canvas=body['canvas'];margin=20.0
    offset=(margin-canvas['x'],margin-canvas['y']);width=canvas['width']+2*margin;height=canvas['height']+2*margin
    glyph_ids={};coordinates={};reaction_glyphs={};curve_coordinates={};text_count=0;line_count=0;reference_count=0
    def position(node):return node['x']+offset[0],node['y']+offset[1]
    def box(glyph,x,y,w,h):
        nonlocal width,height
        b=glyph.getBoundingBox();b.setX(x);b.setY(y);b.setZ(0);b.setWidth(w);b.setHeight(h);b.setDepth(0)
        width=max(width,x+w+margin);height=max(height,y+h+margin)
    def text(glyph_id,label,x,y):
        nonlocal text_count
        t=layout.createTextGlyph();t.setId(identifier(layout_id+'_text_'+str(text_count)))
        t.setGraphicalObjectId(glyph_id);t.setText(label);box(t,x,y-12,max(20,len(label)*7),14);text_count+=1
    for node_id,source_id in metadata.get('source_node_ids',{}).items():
        node=body['nodes'][node_id];x,y=position(node)
        if node['node_type']=='metabolite':
            mid=node['bigg_id']
            if mid not in species_ids:raise ValueError('Unknown original SBML species: '+mid)
            glyph=layout.createSpeciesGlyph();gid=identifier(layout_id+'_species_'+node_id)
            glyph.setId(gid);glyph.setSpeciesId(species_ids[mid]);radius=10 if node['node_is_primary'] else 3
            box(glyph,x-radius,y-radius,2*radius,2*radius)
            text(gid,mid,node['label_x']+offset[0],node['label_y']+offset[1])
            glyph_ids[node_id]=gid;coordinates[gid]=[x,y]
    for map_id,reaction in body['reactions'].items():
        rid=reaction['bigg_id']
        if rid not in reaction_ids:raise ValueError('Unknown original SBML reaction: '+rid)
        glyph=layout.createReactionGlyph();gid=identifier(layout_id+'_reaction_'+map_id)
        glyph.setId(gid);glyph.setReactionId(reaction_ids[rid]);reaction_glyphs[rid]=glyph
        marker=metadata['reactions'][rid]['marker_node_id'];x,y=position(body['nodes'][marker])
        box(glyph,x-4,y-4,8,8);coordinates[gid]=[x,y]
        text(gid,rid,reaction['label_x']+offset[0],reaction['label_y']+offset[1])
    if body['reactions']:
        for rid,target,path in _incidences(data):
            mid=body['nodes'][target]['bigg_id'];coefficient=original_reactions[rid]['metabolites'][mid]
            glyph=reaction_glyphs[rid].createSpeciesReferenceGlyph()
            glyph.setId(identifier(layout_id+'_incidence_'+str(reference_count)))
            glyph.setSpeciesGlyphId(glyph_ids[target])
            glyph.setRole(lib.SPECIES_ROLE_SUBSTRATE if coefficient<0 else lib.SPECIES_ROLE_PRODUCT)
            reaction=original_model.getReaction(reaction_ids[rid])
            references=reaction.getListOfReactants() if coefficient<0 else reaction.getListOfProducts()
            matches=[ref for ref in references if ref.getSpecies()==species_ids[mid]]
            if len(matches)==1 and matches[0].isSetId():glyph.setSpeciesReferenceId(matches[0].getId())
            points=[position(body['nodes'][n]) for n in path]
            if coefficient<0:points.reverse()
            curve_coordinates[glyph.getId()]=points
            for a,b in zip(points,points[1:]):
                if not (math.isclose(a[0],b[0],abs_tol=1e-6,rel_tol=0) or math.isclose(a[1],b[1],abs_tol=1e-6,rel_tol=0)):
                    raise ValueError('SBML Layout export requires orthogonal saved connections')
                segment=glyph.createLineSegment();segment.setStart(a[0],a[1],0);segment.setEnd(b[0],b[1],0);line_count+=1
            reference_count+=1
    dimensions=layout.getDimensions();dimensions.setWidth(width);dimensions.setHeight(height);dimensions.setDepth(0)
    # Transfer only the new Layout XML into the original DOM. This avoids
    # libSBML normalizing or dropping original core/FBC/annotation content.
    generated=minidom.parseString(lib.writeSBMLToString(document))
    new_layout=next(e for e in generated.getElementsByTagNameNS(LAYOUT,'layout') if e.getAttributeNS(LAYOUT,'id')==layout_id or e.getAttribute('id')==layout_id)
    added_attributes=[]
    for a in generated.documentElement.attributes.values():
        if a.namespaceURI==XMLNS and not root.hasAttribute(a.name):
            root.setAttributeNS(XMLNS,a.name,a.value);added_attributes.append(a.name)
    required=prefix+':required'
    if not root.hasAttributeNS(LAYOUT,'required'):
        root.setAttributeNS(LAYOUT,required,'false');added_attributes.append(required)
    model_element=next(e for e in root.childNodes if e.nodeType==Node.ELEMENT_NODE and e.localName=='model')
    lists=[e for e in model_element.childNodes if e.nodeType==Node.ELEMENT_NODE and e.namespaceURI==LAYOUT and e.localName=='listOfLayouts']
    created_list=not lists
    container=lists[0] if lists else original_dom.createElementNS(LAYOUT,prefix+':listOfLayouts')
    if created_list:model_element.appendChild(container)
    container.appendChild(original_dom.importNode(new_layout,True))
    candidate=original_dom.toxml(encoding='utf-8')
    verification=minidom.parseString(candidate)
    attached=next(e for e in verification.getElementsByTagNameNS(LAYOUT,'layout') if e.getAttributeNS(LAYOUT,'id')==layout_id or e.getAttribute('id')==layout_id)
    parent=attached.parentNode;parent.removeChild(attached)
    if created_list:parent.parentNode.removeChild(parent)
    for name in added_attributes:verification.documentElement.removeAttribute(name)
    if _snapshot(verification.documentElement)!=_snapshot(minidom.parseString(raw).documentElement):
        raise RuntimeError('Original SBML XML semantics changed during export')
    checked=lib.readSBMLFromString(candidate.decode());_require_valid(checked,lib,'Cannot reload exported SBML')
    checked.checkConsistency();_require_valid(checked,lib,'Exported SBML consistency failure')
    loaded=checked.getModel().getPlugin('layout').getLayout(layout_id)
    if loaded is None or loaded.getNumSpeciesGlyphs()!=len(glyph_ids) or loaded.getNumReactionGlyphs()!=len(reaction_glyphs):
        raise RuntimeError('Glyphs were lost during SBML round trip')
    for glyph in list(loaded.getListOfSpeciesGlyphs())+list(loaded.getListOfReactionGlyphs()):
        b=glyph.getBoundingBox();center=[b.getX()+b.getWidth()/2,b.getY()+b.getHeight()/2]
        if any(not math.isclose(a,b,abs_tol=1e-6,rel_tol=0) for a,b in zip(center,coordinates[glyph.getId()])):
            raise RuntimeError('Glyph coordinates changed during SBML round trip')
    reloaded_segments=0
    for reaction in loaded.getListOfReactionGlyphs():
        for reference in reaction.getListOfSpeciesReferenceGlyphs():
            expected=curve_coordinates[reference.getId()];curve=reference.getCurve()
            if curve.getNumCurveSegments()!=len(expected)-1:raise RuntimeError('Curve segments lost during SBML round trip')
            for i in range(curve.getNumCurveSegments()):
                segment=curve.getCurveSegment(i)
                for actual,point in ((segment.getStart(),expected[i]),(segment.getEnd(),expected[i+1])):
                    if not (math.isclose(actual.getXOffset(),point[0],abs_tol=1e-6,rel_tol=0) and
                            math.isclose(actual.getYOffset(),point[1],abs_tol=1e-6,rel_tol=0)):
                        raise RuntimeError('Curve coordinates changed during SBML round trip')
                reloaded_segments+=1
    if reloaded_segments!=line_count or loaded.getNumTextGlyphs()!=text_count:
        raise RuntimeError('Layout geometry or text lost during SBML round trip')
    output.parent.mkdir(parents=True,exist_ok=True)
    temporary=None
    try:
        with tempfile.NamedTemporaryFile(dir=output.parent,suffix='.sbml',delete=False) as stream:
            temporary=Path(stream.name);stream.write(candidate)
        if load_model(temporary)!=source:raise RuntimeError('Source model changed after SBML export')
        fba=None
        if verify_fba:
            from run_fba_escher import read_objective,solve_fba
            ids=set(reaction_ids);objective,sense=read_objective(model_path,ids)
            _,before=solve_fba(source,objective,sense);other_objective,other_sense=read_objective(temporary,ids)
            _,after=solve_fba(load_model(temporary),other_objective,other_sense)
            if objective!=other_objective or sense!=other_sense or not math.isclose(before['objective_value'],after['objective_value'],rel_tol=1e-9,abs_tol=1e-9):
                raise RuntimeError('FBA objective changed after SBML export')
            fba=dict(before=before['objective_value'],after=after['objective_value'],unchanged=True)
        temporary.replace(output);temporary=None
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)
    displayed_species={body['nodes'][n]['bigg_id'] for n in glyph_ids}
    warnings=[m for m in _messages(checked) if m['severity']<lib.LIBSBML_SEV_ERROR]
    report=dict(status='valid',format='SBML Level 3 Layout Version 1',render_exported=False,
                source_sbml=str(model_path),source_sha256=hashlib.sha256(raw).hexdigest(),
                output_file=output.name,output_sha256=hashlib.sha256(candidate).hexdigest(),
                level=document.getLevel(),version=document.getVersion(),libsbml_version=lib.getLibSBMLDottedVersion(),
                layout_id=layout_id,existing_layouts_preserved=existing_layouts,
                source_reactions=len(reaction_ids),source_species=len(species_ids),
                species_glyphs=len(glyph_ids),reaction_glyphs=len(reaction_glyphs),text_glyphs=text_count,
                species_reference_glyphs=reference_count,line_segments=line_count,
                undisplayed_reaction_ids=sorted(set(reaction_ids)-set(reaction_glyphs)),
                undisplayed_species_ids=sorted(set(species_ids)-displayed_species),
                coordinate_offset=list(offset),layout_dimensions=[width,height],
                validation=dict(original_xml_semantics_preserved=True,source_model_preserved=True,
                                sbml_consistency_errors=0,glyph_references_and_coordinates_verified=True,
                                curves_horizontal_or_vertical=True),
                warning_count=len(warnings),source_warning_count=len(original_warnings),
                additional_warning_count=max(0,len(warnings)-len(original_warnings)),warnings=warnings[:20],fba_verification=fba,
                limitation='Layout stores geometry only. Colors and arrow styling require Render support in a consuming application.')
    report_path=output.parent/'sbml_layout_report.json'
    report_path.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    if output.parent.resolve()==layout_dir.resolve():
        link=f'<p id="sbml-layout-download"><a href="{html.escape(quote(output.name))}">Download SBML with layout</a></p>'
        for filename in ('index.html','escher_preview.html','escher_fba.html'):
            path=layout_dir/filename
            if path.exists():
                page=path.read_text()
                if 'id="sbml-layout-download"' not in page:
                    path.write_text(page.replace('<header>','<header>'+link,1))
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,required=True)
    parser.add_argument('--layout-dir',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--layout-id',default='sbml2map_layout')
    parser.add_argument('--verify-fba',action='store_true',help='Also compare the full-model FBA objective before and after export')
    args=parser.parse_args()
    try:report=export_sbml_layout(args.model,args.layout_dir,args.output,args.layout_id,args.verify_fba)
    except (ValueError,RuntimeError,OSError) as exc:parser.error(str(exc))
    print(json.dumps({key:report[key] for key in ('output_file','layout_id','species_glyphs','reaction_glyphs','line_segments','status')},indent=2))


if __name__=='__main__':main()
