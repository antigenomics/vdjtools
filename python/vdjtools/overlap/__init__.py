"""vdjtools.overlap — Sample overlap, TCRnet, and sample clustering.

Exact-match overlap metrics live here (:func:`overlap_metrics`, :func:`overlap_pair`);
fuzzy / e-value overlap and TCRnet are **delegated** to the vdjmatch + seqtree engine:

- :func:`fuzzy_overlap` / :func:`fuzzy_overlap_metrics` — edit-distance clonotype
  overlap (``vdjmatch.cluster.overlap``).
- :func:`similarity_overlap` / :func:`similarity_matrix` — sequence-similarity-weighted
  overlap (TINA / Leinster-Cobbold; ``Z=I`` recovers exact overlap, ``Z=1[P≤θ]`` the
  fuzzy one).
- :func:`tcrnet` — neighbourhood-enrichment / convergence test (control-repertoire null)
- :func:`alice` — the same test against a V(D)J generation model (Pgen null)
  (``vdjmatch.evalue.query_evalues``).
- :func:`pairwise_distances` / :func:`cluster_samples` — all-pairs distance matrix and
  MDS / hierarchical embedding.
- :func:`track_clonotypes` — per-sample frequency trajectories.
"""
from .cluster import cluster_samples, pairwise_distances
from .fuzzy import fuzzy_overlap, fuzzy_overlap_metrics
from .metrics import DEFAULT_KEY, overlap_metrics, overlap_pair
from .similarity import SimilarityMatrices, similarity_matrix, similarity_overlap
from .alice import alice
from .tcrnet import tcrnet
from .track import track_clonotypes

__all__ = [
    "overlap_metrics", "overlap_pair", "DEFAULT_KEY",
    "fuzzy_overlap", "fuzzy_overlap_metrics",
    "similarity_overlap", "similarity_matrix", "SimilarityMatrices",
    "alice",
    "tcrnet",
    "pairwise_distances", "cluster_samples",
    "track_clonotypes",
]
