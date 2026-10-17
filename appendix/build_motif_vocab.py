"""Build frozen motif/V dictionaries for the library's statistical feature bases.

Full builds belong on a scheduled compute node. This builder fits reference
parameters only; no repertoire/sample data is read and no distance matrix is cached.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import squareform

from vdjtools.signature.layout import AMINO_ACIDS


def motif_distances(triples: list[str], *, threads: int = 1) -> np.ndarray:
    """One native batch of equal-length, fixed-position BLOSUM62 penalties."""
    from seqtree import SubstitutionMatrix
    from seqtree.gapblock import score_matrix
    from vdjtools.cores import kernel_threads

    budget = kernel_threads(threads)
    return np.asarray(score_matrix(triples, triples, matrix=SubstitutionMatrix.blosum62(), threads=budget))


def cluster_map(names, distances, *, n_clusters: int, prefix: str) -> dict[str, str]:
    """Average-linkage dissimilarity clusters, canonically labelled by first member."""
    names = list(names)
    d = np.asarray(distances)
    n = len(names)
    if len(set(names)) != n or n == 0 or not isinstance(n_clusters, int) or not 1 <= n_clusters <= n:
        raise ValueError('unique nonempty names and 1 <= n_clusters <= vocabulary size required')
    if d.shape != (n, n) or not np.isfinite(d).all() or (d < 0).any():
        raise ValueError('distance table must be square, finite and nonnegative')
    if not np.allclose(d, d.T, rtol=0, atol=1e-6) or not np.allclose(np.diag(d), 0, rtol=0, atol=1e-6):
        raise ValueError('distance table must be symmetric with a zero diagonal')
    order = np.argsort(names)
    sorted_names = [names[i] for i in order]
    d = d[np.ix_(order, order)]
    if n == 1:
        labels = np.zeros(1, dtype=int)
    else:
        condensed = squareform(d, checks=False)
        del d
        labels = cut_tree(linkage(condensed, method='average'), n_clusters=[n_clusters]).ravel()
    groups = sorted(set(labels), key=lambda k: min(name for name, lab in zip(sorted_names, labels) if lab == k))
    canonical = {lab: f'{prefix}{i:03d}' for i, lab in enumerate(groups)}
    return {name: canonical[int(lab)] for name, lab in zip(sorted_names, labels)}


def build_motif_map(*, n_clusters: int = 100, threads: int = 1) -> dict:
    """Cluster all 8,000 triples; retain the mapping, not the pairwise workspace."""
    import scipy
    import seqtree
    triples = [''.join(t) for t in itertools.product(AMINO_ACIDS, repeat=3)]
    distance = motif_distances(triples, threads=threads)
    mapping = cluster_map(triples, distance, n_clusters=n_clusters, prefix='m')
    return {'version': 1, 'mapping': mapping, 'n_clusters': n_clusters,
            'distance': 'seqtree_blosum62_fixed_position_gram_penalty', 'linkage': 'average',
            'amino_acids': AMINO_ACIDS, 'scipy': scipy.__version__, 'seqtree': seqtree.__version__}


def build_v_map(matrix_path, genes, *, organism, locus, n_clusters=16, unclustered_genes=()) -> dict:
    """Cluster genes through representative alleles in a supplied full-V artifact."""
    path = Path(matrix_path)
    if path.stem != f'{organism}_{locus}':
        raise ValueError('V matrix filename must identify the requested organism and locus')
    with np.load(path, allow_pickle=False) as z:
        names = [str(x) for x in z['V_alleles']]
        distance = np.asarray(z['V_dist'])
    # Validate the entire source before selecting a submatrix.
    if distance.shape != (len(names), len(names)) or not np.isfinite(distance).all() or (distance < 0).any():
        raise ValueError('invalid V source distance table')
    if not np.allclose(distance, distance.T, atol=1e-6, rtol=0) or not np.allclose(np.diag(distance), 0, atol=1e-6, rtol=0):
        raise ValueError('V source must be symmetric with zero diagonal')
    if len(set(names)) != len(names):
        raise ValueError('duplicate source allele names')
    positions = {name: i for i, name in enumerate(names)}
    for i, name in enumerate(names):
        for member in name.split(','):
            if member in positions and positions[member] != i:
                if not np.array_equal(distance[positions[member]], distance[i]):
                    raise ValueError('an ambiguity member has nonidentical distance rows')
            positions.setdefault(member, i)
    genes = sorted(set(genes))
    allowed_unclustered = set(unclustered_genes)
    unclustered = []
    reps = {}
    for gene in genes:
        candidates = sorted(a for a in positions if ',' not in a and a.split('*')[0] == gene)
        if not candidates:
            if gene in allowed_unclustered:
                unclustered.append(gene)
                continue
            raise ValueError(f'V matrix has no representative for {gene}')
        reps[gene] = gene + '*01' if gene + '*01' in positions else candidates[0]
    clustered = sorted(reps)
    actual = min(n_clusters, len(clustered))
    indices = [positions[reps[g]] for g in clustered]
    mapping = cluster_map(clustered, distance[np.ix_(indices, indices)], n_clusters=actual, prefix='v')
    mapping.update({gene: '_unknown' for gene in unclustered})
    return {'version': 1, 'mapping': mapping, 'representatives': reps, 'organism': organism,
            'locus': locus, 'n_clusters': actual, 'requested_clusters': n_clusters,
            'distance': 'full_v_self_score_minus_twice_global_blosum62_alignment',
            'linkage': 'average', 'unclustered_genes': unclustered, 'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--clusters', type=int, default=100)
    parser.add_argument('--threads', type=int, default=1, help='Native batch compute threads; 0 uses allocated cores')
    parser.add_argument('--v-distances', type=Path)
    parser.add_argument('--v-clusters', type=int, default=16)
    parser.add_argument('--organism', default='human')
    parser.add_argument('--loci', default='TRA,TRB,TRG,TRD,IGH,IGK,IGL')
    args = parser.parse_args()
    result = {'motif_map': build_motif_map(n_clusters=args.clusters, threads=args.threads), 'v_maps': {}}
    if args.v_distances:
        from vdjtools.signature.features import gene_vocab
        from vdjtools.model.reference import load_germline
        for loc in args.loci.split(','):
            genes = gene_vocab(loc, args.organism)['vus']
            germline = load_germline(loc, args.organism)
            functional = set(germline.filter((germline['segment'] == 'V') &
                                             (germline['functionality'] == 'F'))['gene'])
            # No distance is invented for a nonfunctional gene without a full-V template.
            unclustered = sorted(set(genes) - functional)
            result['v_maps'][loc] = build_v_map(args.v_distances / f'{args.organism}_{loc}.npz', genes,
                                               organism=args.organism, locus=loc, n_clusters=args.v_clusters,
                                               unclustered_genes=unclustered)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, sort_keys=True, indent=1))


if __name__ == '__main__':
    main()
