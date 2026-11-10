"""Read identity hints without changing source chemistry or using drawing coordinates."""
import re
from urllib.parse import unquote, urlparse
from chemical_identity import source_values, uri_identifier, metanetx_id

CD = '{http://www.sbml.org/2001/ns/celldesigner}'
FIELDS = re.compile(r'(?:^|\s)(Name|Abbreviation|ChargedFormula|Formula|Charge):\s*')


def local(tag):
    return tag.rsplit('}', 1)[-1]


def decode(identifier, prefix):
    return re.sub(r'__(\d+)__', lambda m: chr(int(m[1])), identifier.removeprefix(prefix))


def hidden_metabolite(record):
    return record['kind'] == 'inorganic' or record.get('boundary_placeholder', False)


def note_fields(text):
    """Read labeled fields, retaining repeated/conflicting values as candidates."""
    matches = list(FIELDS.finditer(text)); result = {}
    for i, match in enumerate(matches):
        value = text[match.end():matches[i + 1].start() if i + 1 < len(matches) else len(text)].strip()
        # Link captions follow the last notes field in the flattened SBML reader.
        value = re.split(r'\s+VMH (?:Metabolite|Reaction)\b', value)[0].strip()
        if value:
            result.setdefault(match[1], []).append(value)
    return result


def reaction_labels(reaction):
    labels = [('SBML reaction ID', reaction.get('id', '')),
              ('SBML reaction name', reaction.get('name', ''))]
    labels.extend(reaction.get('_reaction_labels', []))
    labels.extend(('SBML notes:SOURCE',value.removeprefix('R_')) for value in source_values(reaction.get('notes','')))
    for item in reaction.get('annotations',[]):
        qualifier=local(item['qualifier']);key=uri_identifier(item['uri'])
        if qualifier in ('is','isRelatedTo') and key and key[0]=='bigg.reaction':
            labels.append(('annotation:bigg.reaction:'+qualifier,key[1].removeprefix('R_')))
    return [(source, value.strip()) for source, value in labels if value.strip()]


def boundary_labels(reaction):
    """Related aliases must not replace an explicit EX/DM/SK designation."""
    labels=[(source,value) for source,value in reaction_labels(reaction)
            if value.startswith(('EX_','DM_','SK_'))]
    explicit=[(source,value) for source,value in labels if not source.endswith(':isRelatedTo')]
    return explicit or labels


def enrich_context(context, root):
    """Enrich an isolated discovery copy; never edit the parsed source model."""
    metabolites = {m['id']: m for m in context['metabolites']}
    reactions = {r['id']: r for r in context['reactions']}
    for entity in root.iter():
        tag = local(entity.tag)
        if tag not in ('species', 'reaction'):
            continue
        identifier = decode(entity.get('id', ''), 'M_' if tag == 'species' else 'R_')
        record = (metabolites if tag == 'species' else reactions).get(identifier)
        if record is None:
            continue
        annotations = [c for c in entity if local(c.tag) == 'annotation']
        if tag == 'reaction':
            labels = []
            for annotation in annotations:
                for extension in annotation.findall(CD + 'extension'):
                    value = extension.findtext(CD + 'name')
                    if value:
                        labels.append(('CellDesigner reaction name', value))
                    kinds = extension.findtext(CD + 'reactionType')
                    if kinds:
                        record['_celldesigner_reaction_type'] = kinds
                    refs = [n.get('species') for n in extension.iter()
                            if n.tag in (CD+'baseReactant', CD+'baseProduct', CD+'reactantLink', CD+'productLink')]
                    if refs:
                        record['_celldesigner_connections_verified'] = all(
                            decode(v, 'M_') in record['metabolites'] for v in refs if v)
            if labels:
                record['_reaction_labels'] = labels
            continue
        hints = []
        for value in source_values(record.get('notes','')):
            source_id=decode(value,'M_')
            source_id=re.sub(r'_[a-z][0-9]?$', '', source_id)
            hints.append(dict(source='SBML notes:SOURCE',value=source_id,kind='identifier'))
        fields = note_fields(record.get('notes', ''))
        for field, kind in [('Abbreviation', 'identifier'), ('Name', 'name')]:
            for value in fields.get(field, []):
                hints.append(dict(source='SBML notes:' + field, value=value, kind=kind))
        for notes in (c for c in entity if local(c.tag) == 'notes'):
            for node in notes.iter():
                href = node.get('href', '')
                url = urlparse(href)
                if url.hostname in ('vmh.life', 'www.vmh.life', 'delta.vmh.life'):
                    match = re.fullmatch(r'/metabolite/([^/]+)/?', url.path)
                    if match:
                        hints.append(dict(source='SBML notes:VMH link', value=unquote(match[1]),
                                          kind='identifier', uri=href))
        for annotation in annotations:
            for identity in annotation.iter(CD + 'speciesIdentity'):
                value = identity.findtext(CD + 'name')
                if value:
                    record['_celldesigner_name'] = value
                    hints.append(dict(source='CellDesigner speciesIdentity/name', value=value, kind='label'))
                record['_celldesigner_class'] = identity.findtext(CD + 'class', '')
        if hints:
            if record.get('name'):
                hints.insert(0, dict(source='SBML species name', value=record['name'], kind='label'))
            record['_identity_hints'] = hints
        formulas = set(fields.get('ChargedFormula', fields.get('Formula', [])))
        if not record.get('chemicalFormula', record.get('formula')) and len(formulas) == 1:
            formula = next(iter(formulas))
            if re.fullmatch(r'(?:[A-Z][a-z]?\d*)+', formula) and not ({'R', 'X', 'Z'} & set(re.findall(r'[A-Z][a-z]?', formula))):
                record['chemicalFormula'] = formula
                record['_formula_evidence'] = dict(source='SBML notes:formula', value=formula)
    # UNKNOWN '.' in the default compartment is a drawing source/sink symbol,
    # not a metabolite. Require explicit exchange/demand labels on every use.
    for mid, metabolite in metabolites.items():
        if (metabolite.get('_celldesigner_class') != 'UNKNOWN' or (metabolite.get('name') or metabolite.get('_celldesigner_name')) != '.'
                or metabolite.get('compartment') != 'default'):
            continue
        uses = [r for r in context['reactions'] if mid in r['metabolites']]
        if uses and all(any(re.match(r'^(EX_|DM_|SK_)', value) for _, value in reaction_labels(r))
                        and len(r['metabolites']) == 2 for r in uses):
            metabolite['_boundary_placeholder'] = True
            metabolite['_boundary_evidence'] = 'CellDesigner UNKNOWN dot in default compartment; all uses explicitly labeled EX/DM/SK'
    # MetaNetX represents exchange as physical species <-> same chemical in
    # an explicit model-boundary compartment. It is hidden only in the view.
    compartments=context.get('compartments',{})
    def keys(m):
        result={uri_identifier(a['uri']) for a in m.get('annotations',[])
                if local(a['qualifier']) in ('is','isRelatedTo')}
        result.discard(None)
        if mid:=metanetx_id(m['id']):result.add(('metanetx.chemical',mid))
        if '@' in m['id']:result.add(('model-local chemical ID',m['id'].split('@',1)[0]))
        return result
    def same_chemical(a,b):
        left,right=metanetx_id(a['id']),metanetx_id(b['id'])
        if left and right:return left==right
        return bool(keys(a)&keys(b))
    for mid,metabolite in metabolites.items():
        comp=metabolite.get('compartment','')
        if comp.upper()!='BOUNDARY' and str(compartments.get(comp,'')).lower()!='model boundary':continue
        uses=[r for r in context['reactions'] if mid in r['metabolites']]
        verified=True
        for reaction in uses:
            others=[m for m in reaction['metabolites'] if m!=mid]
            labels=reaction_labels(reaction)
            explicit=(reaction.get('sbo_term')=='SBO:0000627' or
                      any(v.startswith(('EX_','DM_','SK_')) for _,v in labels))
            if (len(others)!=1 or not explicit or
                    metabolites[others[0]].get('compartment')==comp or
                    reaction['metabolites'][mid]*reaction['metabolites'][others[0]]>=0 or
                    not same_chemical(metabolite,metabolites[others[0]])):
                verified=False;break
        if uses and verified:
            metabolite['_boundary_placeholder']=True
            metabolite['_boundary_evidence']='Explicit model-boundary compartment; same chemical and exchange/SBO evidence on every use'
    return context


def fallback_role_evidence(metabolite, reference):
    """Exact catalog identifiers/names only; no fuzzy or formula-only identity."""
    if metabolite.get('_boundary_placeholder') or metabolite.get('_celldesigner_class') in ('PROTEIN', 'GENE', 'RNA'):
        return []
    identifiers = {}; names = {}
    for role, record in reference['roles'].items():
        for value in set(record['aliases'].get('bigg.metabolite', [])) | {role}:
            identifiers.setdefault(value, set()).add(role)
            legacy = re.sub(r'__(D|L|R|S)$', r'_\1', value)
            identifiers.setdefault(legacy, set()).add(role)
        for value in record['names']:
            name = re.sub(r'[^a-z0-9]', '', value.lower())
            names.setdefault(name, set()).add(role)
    evidence = []
    for hint in metabolite.get('_identity_hints', []):
        value = hint['value'].strip()
        label = re.sub(r'\[[^]]+\]$', '', value).strip()
        candidates = identifiers.get(label, set()) if hint['kind'] in ('label', 'identifier') else set()
        tier = 80
        if not candidates and hint['kind'] in ('label', 'name'):
            candidates = names.get(re.sub(r'[^a-z0-9]', '', label.lower()), set())
            tier = 60
        for role in sorted(candidates):
            evidence.append(dict(role=role, tier=tier, source=hint['source'], value=value,
                                 match='exact catalog identifier/legacy stereo alias' if tier == 80 else 'exact catalog name',
                                 **({'uri': hint['uri']} if 'uri' in hint else {})))
    return evidence
