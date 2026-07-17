"""Similarity-aware (functional) diversity of a single repertoire.

Classical diversity treats every clonotype as equally distinct: two repertoires of 1000 clones
score the same whether those clones are near-identical variants of one motif or 1000 unrelated
rearrangements. That is the wrong null for TCRs, where convergent recombination and clonal
expansion both manufacture near-neighbours.

Leinster & Cobbold (2012, *Ecology* 93(3):477, doi:10.1890/10-2402.1) fix this by folding a
**similarity matrix** ``Z`` into the Hill numbers. With ``p`` the clonotype abundance vector and
``(Zp)_i = Σ_j Z_ij p_j`` the *ordinariness* of clonotype ``i`` (how well-represented its
neighbourhood is):

    ᑫD^Z(p) = ( Σ_i p_i (Zp)_i^(q-1) )^(1/(1-q))      q ≠ 1
    ¹D^Z(p) = exp( − Σ_i p_i ln (Zp)_i )              the q→1 limit
    Rao's Q = Σ_ij p_i p_j (1 − Z_ij) = 1 − pᵀZp      expected dissimilarity of two draws

``Z = I`` (nothing resembles anything but itself) recovers the plain Hill numbers exactly —
richness, exp(Shannon), inverse Simpson at q = 0, 1, 2. That identity is the defining special
case and is what ``test_functional.py`` pins against :mod:`vdjtools.stats.inext`.

This is the single-community counterpart of :func:`vdjtools.overlap.similarity_overlap`, which
applies the same kernel as a *two*-sample bilinear form ``pᵀZq``. Both take ``Z`` from seqtree,
so the kernel is defined in exactly one place.

``q`` orders the profile by how much it cares about rare clonotypes: ``q=0`` counts them fully
(and is the most sensitive to sequencing depth), ``q=2`` is dominated by the expanded ones. Read
the profile, not a single number — that is the whole point of a Hill profile.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..io.schema import JUNCTION_AA


def _profile(p: np.ndarray, zp: np.ndarray, q: float) -> float:
    """One ``ᑫD^Z`` value from abundances ``p`` and ordinariness ``Zp``.

    Zero-abundance clonotypes are dropped rather than special-cased: they contribute nothing at
    any ``q`` and would put ``0·log 0`` / ``0^(-1)`` in the sum. ``(Zp)_i ≥ Z_ii·p_i = p_i > 0``
    on what survives (the kernel diagonal is 1), so no guard against ``log 0`` is needed.
    """
    nz = p > 0
    p, zp = p[nz], zp[nz]
    if p.size == 0:
        return 0.0
    if q == 1:
        return float(np.exp(-np.sum(p * np.log(zp))))
    return float(np.sum(p * zp ** (q - 1)) ** (1.0 / (1.0 - q)))


def functional_diversity(df: pl.DataFrame, *, q=(0, 1, 2), key: "tuple[str, ...]" = (JUNCTION_AA,),
                         kernel: str = "exp", tau: float | None = None, matrix=None,
                         max_penalty: int | None = None, gap_prior="central",
                         gap_open: int | None = None, dense: bool | None = None,
                         weight: str = "freq", threads: int = 0) -> pl.DataFrame:
    """Leinster-Cobbold similarity-aware diversity profile of one repertoire.

    Args:
        df: Clonotype frame (canonical schema), one row per clonotype.
        q: Diversity orders to report. ``q=0`` weights rare clonotypes fully, ``q=1`` is the
            Shannon-equivalent, ``q=2`` is dominated by expanded clones. Non-integer ``q`` works.
        key: Columns forming the clonotype identity; **must include** ``junction_aa``.
        kernel: ``"exp"`` (``Zᵢⱼ = exp(−Pᵢⱼ/τ)``), ``"step"`` (``1[P ≤ max_penalty]``), or
            ``"identity"`` (``Z = I`` — recovers the plain Hill numbers).
        tau: Kernel bandwidth for ``"exp"``.
        matrix: Substitution matrix (default BLOSUM62, via seqtree).
        max_penalty: Alignment-penalty cutoff; the ``"step"`` kernel threshold.
        gap_prior: Gap-placement prior passed to seqtree.
        gap_open: Gap-open penalty passed to seqtree.
        dense: Force the dense (``True``) or sparse (``False``) kernel path.
        weight: ``"freq"`` (relative abundance) or ``"presence"`` (uniform — the unweighted,
            incidence-style profile).
        threads: Worker threads for the seqtree kernel build (``0`` = engine default).

    Returns:
        One row per requested ``q`` with columns ``q``, ``diversity`` (``ᑫD^Z``) and ``rao``
        (Rao's quadratic entropy — constant across ``q``, carried for convenience).

    Raises:
        ImportError: If seqtree (or, for the sparse path, scipy) is missing.
        ValueError: On an unknown ``kernel``/``weight``/``gap_prior``, or a ``key`` without
            ``junction_aa``.

    Example:
        >>> functional_diversity(sample, q=(0, 1, 2))["diversity"].to_list()
        [812.4, 402.1, 233.7]
        >>> # Z = I is the sanity anchor: plain Hill numbers.
        >>> functional_diversity(sample, q=(0,), kernel="identity")["diversity"][0]
        4000.0
    """
    from ..overlap.similarity import _weights, similarity_matrix

    qs = [q] if np.isscalar(q) else list(q)
    # ponytail: similarity_matrix is a pair API, so (df, df) builds z_ab/z_aa/z_bb -- three copies
    #           of one block, i.e. ~3x the seqtree work. Correct, and the shortest diff that
    #           reuses the kernel rather than forking a second definition of it. If profiling a
    #           large repertoire makes this hurt, factor a _self_block(a) out of similarity_matrix
    #           and share it with overlap.pairwise_distances, which re-derives the same block per
    #           pair and wants it too.
    sm = similarity_matrix(df, df, key=key, kernel=kernel, tau=tau, matrix=matrix,
                           max_penalty=max_penalty, gap_prior=gap_prior, gap_open=gap_open,
                           dense=dense, threads=threads)
    p = _weights(sm.freq_a, weight)
    zp = np.asarray(sm.z_aa @ p).ravel()          # ordinariness; .ravel() flattens the sparse path
    rao = 1.0 - float(p @ zp)                     # 1 - pᵀZp, no second kernel pass
    return pl.DataFrame({
        "q": [float(x) for x in qs],
        "diversity": [_profile(p, zp, float(x)) for x in qs],
        "rao": [rao] * len(qs),
    })
