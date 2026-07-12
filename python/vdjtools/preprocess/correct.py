"""Frequency-based sequencing-error correction (polars; neighbour search via seqtree).

Reimplements the legacy vdjtools ``Corrector``. Clonotypes whose ``junction_nt`` sit
within a few substitutions of a much more abundant clonotype are treated as PCR /
sequencing errors of that parent and merged into it.

The ``<= max_mismatches`` substitution neighbour search is delegated to **seqtree**
(the "delegate search to seqtree" convention), which builds an immutable fuzzy
index and returns, per query, every within-budget neighbour and its substitution
count. Only the abundance-ratio decision is done here.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

from ..io.schema import (
    JUNCTION_NT,
    COUNT,
    J_CALL,
    V_CALL,
    normalize,
    recompute_frequency,
)

_SEQTREE_HINT = (
    "seqtree is required for vdjtools.preprocess.correct; install the extra with "
    "`pip install 'vdjtools[preprocess]'` (or `pip install seqtree>=0.3`)."
)


def _corrected_counts(counts: np.ndarray, neighbours: "list[list[tuple[int, int]]]",
                      log_ratio_threshold: float) -> np.ndarray:
    """Compute per-clonotype corrected counts from neighbour lists (legacy core).

    For each clonotype ``i`` with count ``c_i`` this mirrors legacy
    ``Corrector.computeCorrectedCount`` over the original counts:

    - it **absorbs** a neighbour ``j`` (``total += c_j``) when
      ``log10(c_i / c_j) > m * log_ratio_threshold`` (i.e. ``c_j`` is more than
      ``ratio**m`` times smaller — a child error), and
    - it is itself **removed** (result ``0``) as soon as a neighbour ``j`` satisfies
      ``log10(c_i / c_j) < -m * log_ratio_threshold`` (``i`` is the child of a
      bigger parent), returning immediately as the legacy loop does.

    Args:
        counts: Original per-clonotype counts, indexed like ``neighbours``.
        neighbours: For each clonotype, a list of ``(j, m)`` pairs — neighbour
            index ``j`` (``!= i``) and its substitution distance ``m``.
        log_ratio_threshold: ``-log10(ratio)`` from :func:`correct`.

    Returns:
        A new int64 array of corrected counts (``0`` marks a removed error).
    """
    out = counts.astype(np.int64).copy()
    for i, nbrs in enumerate(neighbours):
        c_i = int(counts[i])
        total = c_i
        removed = False
        for j, m in nbrs:
            c_j = int(counts[j])
            log_ratio = math.log10(c_i / c_j)
            if log_ratio > m * log_ratio_threshold:
                total += c_j
            elif log_ratio < -m * log_ratio_threshold:
                removed = True
                break
        out[i] = 0 if removed else total
    return out


def _neighbours(seqs: "list[str]", max_mismatches: int):
    """Return per-sequence ``(j, m)`` substitution neighbours via seqtree."""
    import seqtree as st

    params = st.SearchParams(max_subs=max_mismatches, max_ins=0, max_dels=0)
    hits = st.pairwise_batch(seqs, seqs, params, alphabet="nt")
    return [[(h.ref_id, h.n_subs) for h in row if h.ref_id != i]
            for i, row in enumerate(hits)]


def _correct_block(block: pl.DataFrame, max_mismatches: int,
                   log_ratio_threshold: float) -> pl.DataFrame:
    """Correct one block (all rows share V/J when ``same_vj``); returns kept rows."""
    seqs = [s.upper() for s in block[JUNCTION_NT].to_list()]
    neighbours = _neighbours(seqs, max_mismatches)
    counts = block[COUNT].to_numpy()
    corrected = _corrected_counts(counts, neighbours, log_ratio_threshold)
    block = block.with_columns(pl.Series(COUNT, corrected, dtype=pl.Int64))
    return block.filter(pl.col(COUNT) > 0)


def correct(df: pl.DataFrame, max_mismatches: int = 2, ratio: float = 0.05,
            same_vj: bool = False) -> pl.DataFrame:
    """Merge low-frequency sequencing-error clonotypes into their parents.

    Reimplements ``Corrector``. Clonotype pairs whose ``junction_nt`` are within
    ``max_mismatches`` substitutions are compared by abundance: a smaller clonotype
    is merged into a larger one when its count is below ``ratio ** m`` times the
    larger's (``m`` = number of substitutions). All decisions use the *original*
    counts (a single, order-independent pass, as in the legacy parallel stream).

    Args:
        df: A clonotype frame with ``junction_nt`` and ``duplicate_count`` columns.
        max_mismatches: Maximum substitutions for two clonotypes to be neighbours
            (legacy default 2); insertions/deletions are not considered.
        ratio: Per-mismatch parent/child count ratio (legacy default 0.05). A child
            is merged when ``child_count < ratio ** m * parent_count``.
        same_vj: If ``True`` (the opt-in "match-segment" mode, legacy ``-a`` /
            ``--match-segment``) only clonotypes sharing the exact ``v_call`` and
            ``j_call`` are compared; if ``False`` (default, matching legacy
            fidelity — legacy ``Corrector`` is segment-agnostic by default) all
            clonotypes are compared regardless of segment.

    Returns:
        The corrected frame (errors dropped, parents' counts increased), sorted by
        descending ``duplicate_count`` with ``frequency`` recomputed. Rows with a
        null ``junction_nt`` pass through uncorrected.

    Raises:
        ImportError: If seqtree is not installed (see the ``preprocess`` extra).
    """
    try:
        import seqtree  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without seqtree
        raise ImportError(_SEQTREE_HINT) from exc

    log_ratio_threshold = -math.log10(ratio)
    d = normalize(df)

    # Only ACGT nucleotide junctions can be searched; the rest pass through.
    valid = pl.col(JUNCTION_NT).is_not_null() & pl.col(JUNCTION_NT).str.to_uppercase().str.contains(
        r"^[ACGT]+$"
    )
    passthrough = d.filter(~valid)
    work = d.filter(valid)

    blocks = []
    if same_vj:
        for _, block in work.group_by([V_CALL, J_CALL], maintain_order=True):
            blocks.append(_correct_block(block, max_mismatches, log_ratio_threshold))
    elif work.height:
        blocks.append(_correct_block(work, max_mismatches, log_ratio_threshold))

    kept = [b for b in blocks if b.height]
    out = pl.concat([*kept, passthrough], how="vertical_relaxed") if kept else passthrough
    return recompute_frequency(out.sort(COUNT, descending=True, maintain_order=True))
