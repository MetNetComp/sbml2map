#!/usr/bin/env python3
"""Build stoichiometric and adjacency matrices from BiGG JSON or SBML.

Uses only the Python standard library. Rows of S are metabolites; columns are
reactions. Negative values are substrates; positive values are products.
For B = (S != 0), undirected adjacency is (B B^T > 0) for metabolites and
(B^T B > 0) for reactions. Both have zero diagonals.
Directed metabolite edges link substrates to products in each allowed direction.
Directed reaction edges link a producer to a consumer of the same metabolite.
A[i][j] = 1 means i -> j.

These edges show structure, not steady-state flux feasibility or atom transfer.
ATP, water, and other compounds are included as ordinary metabolites.
"""

import argparse
import csv
from pathlib import Path
from model_input import load_model


def stoichiometric_matrix(model):
    """Return (S, metabolite_ids, reaction_ids) from a model dictionary."""
    metabolites = [m['id'] for m in model['metabolites']]
    reactions = [r['id'] for r in model['reactions']]
    index = {name: i for i, name in enumerate(metabolites)}
    S = [[0.0] * len(reactions) for _ in metabolites]
    for j, reaction in enumerate(model['reactions']):
        for name, coefficient in reaction['metabolites'].items():
            S[index[name]][j] = float(coefficient)
    return S, metabolites, reactions


def steady_state_matrix(model):
    """Balance internal species only; retain every original reaction column.

    Missing boundaryCondition means false, preserving legacy JSON models.
    Drawing-only placeholders do not imply exemption from mass balance.
    """
    S, mids, rids = stoichiometric_matrix(model)
    balanced = []
    for i, metabolite in enumerate(model['metabolites']):
        boundary = metabolite.get('boundaryCondition', False)
        if not isinstance(boundary, bool):
            raise ValueError('boundaryCondition must be boolean: ' + metabolite['id'])
        if not boundary:
            balanced.append(i)
    return [S[i] for i in balanced], [mids[i] for i in balanced], rids


def adjacency_from_stoichiometry(S, network='metabolite', directed=False,
                                 bounds=None):
    """Return a 0/1 adjacency matrix from S, with no self-loops.

    network: 'metabolite' or 'reaction'.
    bounds: (lower_bound, upper_bound) pairs in reaction order.
      Used only for directed graphs. If omitted, all reactions are forward.
      An upper bound > 0 allows forward flow; a lower bound < 0 allows reverse.
      Reactions with both bounds at zero are isolated in directed graphs.
      Undirected graphs use nonzero entries of S, regardless of bounds.
    """
    if network not in ('metabolite', 'reaction'):
        raise ValueError('network must be metabolite or reaction')
    if not S or not S[0] or any(len(row) != len(S[0]) for row in S):
        raise ValueError('S must be a nonempty rectangular matrix')
    n_met, n_rxn = len(S), len(S[0])
    size = n_met if network == 'metabolite' else n_rxn
    A = [[0] * size for _ in range(size)]

    if not directed:
        # Use nonzero sets to avoid cancellation from signed matrix products.
        groups = (
            [{i for i in range(n_met) if S[i][j] != 0}
             for j in range(n_rxn)] if network == 'metabolite' else
            [{j for j in range(n_rxn) if S[i][j] != 0}
             for i in range(n_met)]
        )
        for group in groups:
            for source in group:
                for target in group:
                    if source != target:
                        A[source][target] = 1
        return A

    if bounds is None:
        bounds = [(0, 1)] * n_rxn
    if len(bounds) != n_rxn or any(lo > hi for lo, hi in bounds):
        raise ValueError('bounds must contain one valid (lower, upper) per reaction')

    consumers = [set() for _ in range(n_met)]
    producers = [set() for _ in range(n_met)]
    for j, (lower, upper) in enumerate(bounds):
        substrates = {i for i in range(n_met) if S[i][j] < 0}
        products = {i for i in range(n_met) if S[i][j] > 0}
        directions = []
        if upper > 0:
            directions.append((substrates, products))
        if lower < 0:
            directions.append((products, substrates))
        for sources, targets in directions:
            if network == 'metabolite':
                for source in sources:
                    for target in targets:
                        if source != target:
                            A[source][target] = 1
            else:
                for i in sources:
                    consumers[i].add(j)
                for i in targets:
                    producers[i].add(j)

    if network == 'reaction':
        for producing, consuming in zip(producers, consumers):
            for source in producing:
                for target in consuming:
                    if source != target:
                        A[source][target] = 1
    return A


def write_matrix(path, matrix, row_ids, column_ids):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['id', *column_ids])
        for name, row in zip(row_ids, matrix):
            writer.writerow([name, *row])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path,
                        default=Path(__file__).with_name('e_coli_core.json'))
    parser.add_argument('--network', choices=['metabolite', 'reaction'],
                        default='metabolite')
    parser.add_argument('--directed', action='store_true',
                        help='Build a directed adjacency matrix using reaction bounds')
    parser.add_argument('--output-dir', type=Path, default=Path('network_output'))
    args = parser.parse_args()
    model = load_model(args.model)
    S, metabolites, reactions = stoichiometric_matrix(model)
    bounds = [(r['lower_bound'], r['upper_bound']) for r in model['reactions']]
    A = adjacency_from_stoichiometry(S, args.network, args.directed, bounds)
    ids = metabolites if args.network == 'metabolite' else reactions
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.model.stem
    s_path = args.output_dir / f'{prefix}_stoichiometric.csv'
    mode = 'directed' if args.directed else 'undirected'
    a_path = args.output_dir / f'{prefix}_{args.network}_{mode}_adjacency.csv'
    write_matrix(s_path, S, metabolites, reactions)
    write_matrix(a_path, A, ids, ids)
    edges = sum(map(sum, A)) // (1 if args.directed else 2)
    print(f'S: {len(metabolites)} x {len(reactions)} -> {s_path}')
    print(f'A: {len(ids)} x {len(ids)}, edges: {edges} -> {a_path}')


if __name__ == '__main__':
    main()
