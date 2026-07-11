"""vdjtools.overlap — Sample overlap and TCRnet.

Exact-match overlap metrics live here (:func:`overlap_metrics`, :func:`overlap_pair`);
fuzzy / e-value overlap and TCRnet are **delegated** to the vdjmatch + seqtree engine:

- :func:`fuzzy_overlap` / :func:`fuzzy_overlap_metrics` — edit-distance clonotype
  overlap (``vdjmatch.cluster.overlap``).
- :func:`tcrnet` — neighbourhood-enrichment / convergence test
  (``vdjmatch.evalue.query_evalues``).
"""
from .fuzzy import fuzzy_overlap, fuzzy_overlap_metrics
from .metrics import DEFAULT_KEY, overlap_metrics, overlap_pair
from .tcrnet import tcrnet

__all__ = [
    "overlap_metrics", "overlap_pair", "DEFAULT_KEY",
    "fuzzy_overlap", "fuzzy_overlap_metrics",
    "tcrnet",
]
