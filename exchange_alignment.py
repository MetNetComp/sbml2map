"""Evidence-based EX/extracellular/cytosolic triplets and collinearity checks."""

EPS = 1e-7


def _compartment(metabolite):
    compartment = str(metabolite.get('compartment', '')).lower()
    if compartment in ('e', 'extracellular', 'external', 'extracellular_space'):
        return 'e'
    if compartment in ('c', 'cytosol', 'cytoplasm', 'cytosolic'):
        return 'c'
    return None


def _base(metabolite, compartment):
    identifier = metabolite['id']
    raw = str(metabolite.get('compartment', compartment))
    for suffix in ('_'+raw, '['+raw+']', '_'+compartment):
        if identifier.endswith(suffix):
            return identifier[:-len(suffix)]
    return None


def _identities(metabolite):
    identities = set()
    # Only exact chemical identities are used; broad isRelatedTo links do not
    # establish that two compartments contain the same chemical.
    for item in metabolite.get('annotations', ()):
        if not isinstance(item, dict):
            continue
        qualifier = item.get('qualifier', item.get('relation', '')).rsplit('}',1)[-1]
        if qualifier not in ('is', 'bqbiol:is', 'BQB_IS'):
            continue
        value = item.get('resource', item.get('uri', ''))
        if value:
            from annotation_context import uri_identifier
            identity = uri_identifier(value)
            if identity and identity[0] in ('bigg.metabolite','kegg.compound','chebi','metanetx.chemical'):
                identities.add(identity)
    return identities


def discover_alignment(model, nodes, exchange_nodes):
    metabolites = {m['id']: m for m in model['metabolites']}
    reactions = {r['id']: r for r in model['reactions']}
    cytosolic = [m for m in metabolites.values() if _compartment(m) == 'c' and 'met:'+m['id'] in nodes]
    groups, skipped = [], []
    for key in sorted(exchange_nodes):
        rid = nodes[key]['original_id']
        reaction = reactions[rid]
        extracellular = [metabolites[mid] for mid in reaction['metabolites']
                         if mid in metabolites and _compartment(metabolites[mid]) == 'e' and 'met:'+mid in nodes]
        if len(extracellular) != 1:
            skipped.append(dict(reaction_id=rid, reason='No unique displayed extracellular participant'))
            continue
        external = extracellular[0]
        identities = _identities(external)
        base = _base(external, 'e')
        exact = [m for m in cytosolic if identities and identities & _identities(m)]
        named = [m for m in cytosolic if base is not None and _base(m, 'c') == base]
        matches = exact if exact else named
        if len(matches) != 1:
            skipped.append(dict(reaction_id=rid, external_metabolite_id=external['id'],
                                reason='No unique displayed cytosolic chemical match', candidates=[m['id'] for m in matches]))
            continue
        groups.append(dict(reaction_id=rid, exchange=key, external='met:'+external['id'],
                           internal='met:'+matches[0]['id'],
                           evidence='Exact chemical identity annotation' if exact else 'Matching compartment-normalized source IDs'))
    # An external species shared by multiple EX reactions or multiple matched
    # intracellular identities has no unique three-point layout specification.
    duplicate = {g['external'] for g in groups if sum(h['external'] == g['external'] for h in groups) > 1}
    for group in groups:
        if group['external'] in duplicate:
            skipped.append(dict(reaction_id=group['reaction_id'], reason='Extracellular species has multiple exchange reactions'))
    return [g for g in groups if g['external'] not in duplicate], skipped


def triplet_face(group, xy, cols, rows):
    a,b,c = (xy[group[k]] for k in ('exchange','external','internal'))
    if abs(a[1]-b[1]) < EPS and abs(b[1]-c[1]) < EPS:
        if abs(a[0]) < EPS and a[0]+.16-EPS <= b[0] <= c[0]-.16+EPS:
            return 'left'
        if abs(a[0]-(cols-1)) < EPS and c[0]+.16-EPS <= b[0] <= a[0]-.16+EPS:
            return 'right'
    if abs(a[0]-b[0]) < EPS and abs(b[0]-c[0]) < EPS:
        if abs(a[1]) < EPS and a[1]+.16-EPS <= b[1] <= c[1]-.16+EPS:
            return 'bottom'
        if abs(a[1]-(rows-1)) < EPS and c[1]+.16-EPS <= b[1] <= a[1]-.16+EPS:
            return 'top'
    return None


def validate_alignment(groups, xy, cols, rows):
    results = [dict(group, face=triplet_face(group, xy, cols, rows)) for group in groups]
    failures = [r for r in results if r['face'] is None]
    if failures:
        raise RuntimeError('Exchange alignment constraints violated: '+str([r['reaction_id'] for r in failures]))
    return dict(verified=True, checked_triplets=len(groups), aligned_triplets=results,
                ordering='EX -> extracellular -> cytosolic, from outer boundary inward')
