"""Normalize chemical identifiers and read optional, versioned MNXref tables."""
from collections import defaultdict
import hashlib
from pathlib import Path
import re
from urllib.parse import unquote, urlparse


def normalize_identifier(namespace, value):
    namespace=namespace.lower().strip()
    namespace={'inchikey':'inchi_key', 'biggm':'bigg.metabolite',
               'keggc':'kegg.compound'}.get(namespace,namespace)
    value=unquote(value).strip()
    if namespace=='chebi':
        value='CHEBI:'+re.sub(r'^CHEBI:', '', value, flags=re.I)
    return namespace,value


def uri_identifier(uri):
    uri=unquote(uri)
    if uri.startswith('urn:miriam:'):
        namespace,_,value=uri[len('urn:miriam:'):].partition(':')
    else:
        parsed=urlparse(uri)
        if parsed.hostname not in ('identifiers.org','www.identifiers.org'):
            return None
        part=parsed.path.strip('/')
        namespace,_,value=part.partition('/' if '/' in part else ':')
    return normalize_identifier(namespace,value) if value else None


def source_values(notes):
    match=re.search(r'(?:^|\s)SOURCE:\s*(.*?)(?=\s+[A-Z][A-Z _-]*:\s|$)', notes)
    return [s.strip() for s in match[1].split(';') if s.strip()] if match else []


def metanetx_id(identifier):
    # model_input has already decoded __64__ to @; strip only the compartment.
    match=re.fullmatch(r'(MNXM\d+)(?:@[^@]+)?',identifier)
    return match[1] if match else None


class MetaNetXReference:
    """Load local xrefs and deprecated IDs; no implicit network requests."""
    def __init__(self,reference,directory=None):
        self.roles=defaultdict(set);self.replacements=defaultdict(set)
        self.provenance=None
        external=defaultdict(set)
        for role,record in reference['roles'].items():
            for namespace,values in record['aliases'].items():
                for value in values:
                    key=normalize_identifier(namespace,value)
                    external[key].add(role)
                    if key[0]=='metanetx.chemical':self.roles[key[1]].add(role)
        if directory is None:return
        directory=Path(directory)
        paths=[directory/name for name in ('chem_xref.tsv','chem_depr.tsv')]
        for path in paths:
            if not path.is_file():raise ValueError('Missing MetaNetX reference table: '+str(path))
        files={}
        def rows(path):
            digest=hashlib.sha256();headers=[]
            with path.open('rb') as stream:
                for raw in stream:
                    digest.update(raw);line=raw.decode('utf-8').rstrip('\r\n')
                    if line.startswith('#'):
                        if len(headers)<30:headers.append(line)
                        continue
                    if line:yield line.split('\t')
            files[path.name]=dict(sha256=digest.hexdigest(),headers=headers)
        for row in rows(paths[0]):
            if len(row)<2:raise ValueError('Invalid MetaNetX chem_xref.tsv row')
            namespace,_,value=row[0].partition(':')
            for role in external.get(normalize_identifier(namespace,value),()):
                self.roles[row[1]].add(role)
        for row in rows(paths[1]):
            if len(row)<2:raise ValueError('Invalid MetaNetX chem_depr.tsv row')
            if row[0]!=row[1]:self.replacements[row[0]].add(row[1])
        self.provenance=dict(directory=str(directory.resolve()),files=files)
        expanded=defaultdict(set)
        for mid,roles in self.roles.items():
            terminals=self.terminals(mid)
            # A split must not copy an obsolete compound's identity onto all
            # of its chemically distinct replacements.
            if len(terminals)==1:
                expanded[next(iter(terminals))].update(roles)
        self.roles=expanded

    def terminals(self,mid):
        # Explicit stack handles long replacement chains; cycles are invalid.
        result=set();pending=[(mid,frozenset())]
        while pending:
            current,seen=pending.pop()
            if current in seen:raise ValueError('Cycle in MetaNetX deprecated identifiers: '+current)
            children=self.replacements.get(current)
            if not children:result.add(current)
            else:pending.extend((child,seen|{current}) for child in children)
        return result

    def candidates(self,mid):
        return set().union(*(self.roles.get(t,set()) for t in self.terminals(mid)))
