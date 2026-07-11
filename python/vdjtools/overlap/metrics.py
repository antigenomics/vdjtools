"""Exact-match pairwise repertoire-overlap metrics (pure polars + numpy).

Implements the four legacy vdjtools overlap metrics (``OverlapEvaluator``) on an
exact clonotype match key. This is exact-match only; fuzzy / e-value overlap and
TCRnet are delegated to vdjmatch (``cluster.overlap`` / ``evalue.query_evalues``).
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import CDR3_AA, COUNT, J_CALL, V_CALL

#: Default exact match key (legacy "Strict" overlap: CDR3 aa + V + J).
DEFAULT_KEY = (CDR3_AA, V_CALL, J_CALL)


def _aggregate(df: pl.DataFrame, key: list[str]) -> pl.DataFrame:
    """Collapse to unique clonotype keys, summing counts and recomputing frequency."""
    agg = df.group_by(key, maintain_order=True).agg(pl.col(COUNT).sum().alias("_count"))
    total = agg["_count"].sum()
    total = total if total else 1
    return agg.with_columns((pl.col("_count") / total).alias("_freq"))


def overlap_pair(a: pl.DataFrame, b: pl.DataFrame,
                 key: "tuple[str, ...]" = DEFAULT_KEY) -> "tuple[pl.DataFrame, dict]":
    """Compute the shared-clonotype table and overlap metrics for two samples.

    Both frames are first collapsed to unique clonotype keys (summing counts;
    frequencies are recomputed within each sample). The four legacy metrics are:

    - **D** (diversity): ``d12 / (d1 * d2)``.
    - **F** (frequency): ``sqrt(Σ_shared f_a · Σ_shared f_b)``.
    - **F2**: ``Σ_shared sqrt(f_a · f_b)``.
    - **R**: Pearson correlation of ``log10`` shared frequencies. ``log10`` is used
      because clonotype frequencies span orders of magnitude (legacy plots the scatter
      on log axes); requires ``d12 > 2`` else ``R = 0`` (matching the legacy guard).

    Args:
        a: First clonotype frame.
        b: Second clonotype frame.
        key: Columns forming the exact match key (default
            ``("cdr3_aa", "v_call", "j_call")``; use ``("cdr3_aa",)`` for CDR3-only
            or add ``"cdr3_nt"`` for nucleotide-level matching).

    Returns:
        A tuple ``(shared, metrics)`` where ``shared`` is a ``pl.DataFrame`` of the
        joined shared clonotypes (key columns plus ``count_a, count_b, freq_a,
        freq_b``) and ``metrics`` is a dict with keys ``D, F, F2, R, d1, d2, d12``.
    """
    key = list(key)
    a_agg = _aggregate(a, key)
    b_agg = _aggregate(b, key)
    d1, d2 = a_agg.height, b_agg.height

    shared = a_agg.join(b_agg, on=key, how="inner", suffix="_b").rename({
        "_count": "count_a", "_freq": "freq_a", "_count_b": "count_b", "_freq_b": "freq_b",
    })
    d12 = shared.height

    fa = shared["freq_a"].to_numpy()
    fb = shared["freq_b"].to_numpy()
    div = float(d12) / (d1 * d2) if d1 and d2 else 0.0
    freq = float(np.sqrt(fa.sum() * fb.sum())) if d12 else 0.0
    freq2 = float(np.sum(np.sqrt(fa * fb))) if d12 else 0.0
    if d12 > 2:
        r = float(np.corrcoef(np.log10(fa), np.log10(fb))[0, 1])
        if np.isnan(r):
            r = 0.0
    else:
        r = 0.0

    metrics = {"D": div, "F": freq, "F2": freq2, "R": r, "d1": d1, "d2": d2, "d12": d12}
    return shared.select([*key, "count_a", "count_b", "freq_a", "freq_b"]), metrics


def overlap_metrics(a: pl.DataFrame, b: pl.DataFrame,
                    key: "tuple[str, ...]" = DEFAULT_KEY) -> dict:
    """Compute the four exact-match overlap metrics (D, F, F2, R) for two samples.

    Args:
        a: First clonotype frame.
        b: Second clonotype frame.
        key: Exact match key (see :func:`overlap_pair`).

    Returns:
        Dict with keys ``D, F, F2, R, d1, d2, d12``.
    """
    return overlap_pair(a, b, key)[1]
