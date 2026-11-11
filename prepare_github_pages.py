#!/usr/bin/env python3
"""Prepare the E. coli core static site for GitHub Pages (main, /docs)."""
from pathlib import Path
import csv
import json
import shutil
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit

from highlighting import grid_elements, inject_controls

ROOT = Path(__file__).resolve().parent
FILES = (
    'index.html', 'escher_preview.html', 'layout.png', 'layout.svg',
    'orange_routes.png', 'orange_routes.svg', 'report.json',
    'reaction_ranking.csv', 'escher_map.json', 'escher_flux_medians.json',
)

FBA_FILES = ('escher_fba.html','escher_fba_map.json','fba_report.json','fba_fluxes.csv','fba_fluxes.json','escher_fluxes.json')
SBML_FILES = ('layout.sbml','sbml_layout_report.json')



class _ModelLinks(HTMLParser):
    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in ('href', 'src', 'download') and value:
                path = unquote(urlsplit(value).path).lower()
                if path.endswith(('.sbml', '.xml')) or '.sbml.' in path:
                    raise ValueError('iCH360 must not link to SBML/XML model downloads: ' + value)


def check_ich360_publication(docs):
    """Keep iCH360 source/positioned SBML files and model links off the site."""
    for path in docs.rglob('*'):
        if 'ich' not in str(path.relative_to(docs)).lower() or not path.is_file():
            continue
        name = path.name.lower()
        if name.endswith(('.sbml', '.xml')) or '.sbml.' in name:
            raise ValueError('iCH360 SBML/XML model must not be published: ' + str(path))
        if path.suffix.lower() == '.html':
            _ModelLinks().feed(path.read_text())


def main():
    check_ich360_publication(ROOT / 'docs')
    source = ROOT / 'results/e_coli_core'
    target = ROOT / 'docs/e_coli_core'
    for name in (*FILES, 'source_model.json', 'coordinates.csv', 'edges.csv', 'reaction_ranking.json'):
        if not (source / name).is_file():
            raise SystemExit(f'Missing drawing file: {name}')
    model = json.loads((source / 'source_model.json').read_text())
    records = json.loads((source / 'reaction_ranking.json').read_text())
    with (source / 'coordinates.csv').open() as stream:
        nodes = {row['id']: row for row in csv.DictReader(stream)}
    with (source / 'edges.csv').open() as stream:
        edges = list(csv.DictReader(stream))
    elements = grid_elements(nodes, edges)
    target.mkdir(parents=True, exist_ok=True)
    files = FILES + tuple(name for name in FBA_FILES + SBML_FILES + ('parallel_routing.json',) if (source/name).is_file())
    for name in files:
        if name.endswith('.html'):
            page = inject_controls((source / name).read_text(), model, elements,
                                   records=records, escher=name in ('escher_preview.html','escher_fba.html'))
            (target / name).write_text(page)
        else:
            shutil.copyfile(source / name, target / name)
    (ROOT / 'docs/.nojekyll').touch()
    print(f'Prepared {len(files)} drawing files in {target}')


if __name__ == '__main__':
    main()
