#!/usr/bin/env python3
"""Prepare the E. coli core static site for GitHub Pages (main, /docs)."""
from pathlib import Path
import csv
import json
import shutil

from highlighting import grid_elements, inject_controls

ROOT = Path(__file__).resolve().parent
FILES = (
    'index.html', 'escher_preview.html', 'layout.png', 'layout.svg',
    'orange_routes.png', 'orange_routes.svg', 'report.json',
    'reaction_ranking.csv', 'escher_map.json', 'escher_flux_medians.json',
)

FBA_FILES = ('escher_fba.html','fba_report.json','fba_fluxes.csv','fba_fluxes.json','escher_fluxes.json')
SBML_FILES = ('layout.sbml','sbml_layout_report.json')


def main():
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
    files = FILES + tuple(name for name in FBA_FILES + SBML_FILES if (source/name).is_file())
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
