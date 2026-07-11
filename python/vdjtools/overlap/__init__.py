"""vdjtools.overlap — Sample overlap and TCRnet.

Exact-match overlap metrics live here (:func:`overlap_metrics`, :func:`overlap_pair`);
fuzzy / e-value overlap and TCRnet are delegated to vdjmatch (``cluster.overlap`` /
``evalue.query_evalues``).
"""
from .metrics import DEFAULT_KEY, overlap_metrics, overlap_pair

__all__ = ["overlap_metrics", "overlap_pair", "DEFAULT_KEY"]
