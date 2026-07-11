"""Diversity estimators on a clonotype count vector.

Every estimator takes the per-clonotype ``duplicate_count`` vector of a single
clonotype frame (one row per clonotype). Definitions follow legacy vdjtools
(``com.antigenomics.vdjtools.diversity``); each deviation is noted in its docstring.

Notation: ``n`` = total reads, ``Sobs`` = observed richness (# clonotypes),
``f_i = count_i / n``, ``F1`` / ``F2`` = number of singleton / doubleton clonotypes.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

from ..io.schema import COUNT


def _counts(df: pl.DataFrame) -> np.ndarray:
    """Extract the positive integer clonotype count vector from a frame."""
    c = df[COUNT].drop_nulls().to_numpy().astype(np.int64)
    return c[c > 0]


def observed_richness(counts: np.ndarray) -> int:
    """Observed species richness ``Sobs`` — the number of clonotypes.

    Args:
        counts: Per-clonotype count vector.

    Returns:
        Number of clonotypes.
    """
    return int(counts.size)


def _f1_f2(counts: np.ndarray) -> tuple[int, int]:
    """Return (singletons F1, doubletons F2)."""
    return int(np.count_nonzero(counts == 1)), int(np.count_nonzero(counts == 2))


def chao1(counts: np.ndarray) -> float:
    """Chao1 lower-bound richness estimate.

    ``Chao1 = Sobs + F1*(F1-1) / (2*(F2+1))`` (bias-corrected form, legacy default).

    Args:
        counts: Per-clonotype count vector.

    Returns:
        The Chao1 richness estimate.
    """
    sobs = observed_richness(counts)
    f1, f2 = _f1_f2(counts)
    return sobs + f1 * (f1 - 1) / (2 * (f2 + 1))


def chao_e(counts: np.ndarray, extrapolate_to: int | None = None) -> float:
    """Chao extrapolated richness at a larger sampling depth.

    Follows legacy ``ChaoEstimator.chaoE``:
    ``S = Sobs + F0 * [1 - (1 - F1/(n*F0))^(m* )]`` where ``F0 = F1*(F1-1)/(2*(F2+1))``
    and ``m* = extrapolate_to - n``.

    Args:
        counts: Per-clonotype count vector.
        extrapolate_to: Target read depth (``>= n``). Defaults to ``2*n`` — there is
            no cohort here to borrow "largest sample" from as legacy does, so the
            sample is extrapolated to double its size.

    Returns:
        The extrapolated richness. Falls back to ``Sobs`` when extrapolation is
        undefined (``F0 == 0``, i.e. too few singletons/doubletons).

    Raises:
        ValueError: If ``extrapolate_to < n``.
    """
    n = int(counts.sum())
    sobs = observed_richness(counts)
    f1, f2 = _f1_f2(counts)
    f0 = f1 * (f1 - 1) / (2 * (f2 + 1))
    target = 2 * n if extrapolate_to is None else int(extrapolate_to)
    if target < n:
        raise ValueError(f"extrapolate_to ({target}) must be >= n ({n})")
    if f0 == 0 or n == 0:
        return float(sobs)
    m_star = target - n
    brackets = 1.0 - (1.0 - f1 / (n * f0)) ** m_star
    return sobs + f0 * brackets


def efron_thisted(counts: np.ndarray, max_depth: int = 20,
                  cv_threshold: float = 0.05) -> float:
    """Efron–Thisted total-diversity lower-bound estimate.

    Direct port of legacy ``ExactEstimator.getEfronThisted``: increases the
    truncation depth until the coefficient of variation ``D/S`` reaches
    ``cv_threshold``, accumulating alternating binomial coefficients over the
    count-frequency spectrum ``f_x`` (number of clonotypes seen exactly ``x`` times).

    Args:
        counts: Per-clonotype count vector.
        max_depth: Maximum truncation depth.
        cv_threshold: Stop once ``std/estimate`` reaches this value.

    Returns:
        The Efron–Thisted richness estimate.
    """
    sobs = observed_richness(counts)
    max_count = int(counts.max()) if counts.size else 0
    fx = {int(k): int(v) for k, v in zip(*np.unique(counts, return_counts=True))}
    s = float(sobs)
    for depth in range(1, max_depth + 1):
        h = [0.0] * depth
        nx = [float(fx.get(y, 0)) for y in range(1, depth + 1)]
        for y in range(1, depth + 1):
            for x in range(1, y + 1):
                coef = math.comb(y - 1, x - 1)
                h[x - 1] += coef if x % 2 == 1 else -coef
        s = sobs + sum(h[i] * nx[i] for i in range(depth))
        d = math.sqrt(sum(h[i] * h[i] * nx[i] for i in range(depth)))
        if depth >= max_count:
            break
        if s != 0 and d / s >= cv_threshold:
            break
    return float(s)


def _shannon_entropy(counts: np.ndarray) -> float:
    """Shannon entropy ``H = -Σ f_i ln f_i`` in nats."""
    n = counts.sum()
    if n == 0:
        return 0.0
    f = counts / n
    return float(-np.sum(f * np.log(f)))


def shannon_wiener(counts: np.ndarray) -> float:
    """Shannon–Wiener diversity ``exp(H)`` (effective number of clonotypes).

    Matches legacy ``shannonWienerIndex`` = ``exp(-Σ f_i ln f_i)`` (the Hill number
    of order 1), i.e. the exponential of Shannon entropy rather than the entropy
    itself.

    Args:
        counts: Per-clonotype count vector.

    Returns:
        ``exp(H)``.
    """
    return float(math.exp(_shannon_entropy(counts)))


def normalized_shannon_wiener(counts: np.ndarray) -> float:
    """Normalised Shannon–Wiener index ``H / ln(Sobs)`` (Pielou evenness).

    Matches legacy ``normalizedShannonWienerIndex``. This normalises the *entropy*
    ``H`` (not ``exp(H)``) by ``ln(Sobs)``, so it lies in ``[0, 1]``.

    Args:
        counts: Per-clonotype count vector.

    Returns:
        ``H / ln(Sobs)``; ``0.0`` when ``Sobs <= 1`` (no diversity to normalise).
    """
    sobs = observed_richness(counts)
    if sobs <= 1:
        return 0.0
    return _shannon_entropy(counts) / math.log(sobs)


def inverse_simpson(counts: np.ndarray) -> float:
    """Inverse Simpson index ``1 / Σ f_i^2`` (Hill number of order 2).

    Args:
        counts: Per-clonotype count vector.

    Returns:
        The inverse Simpson index.
    """
    n = counts.sum()
    if n == 0:
        return 0.0
    f = counts / n
    return float(1.0 / np.sum(f * f))


def d50(counts: np.ndarray, fraction: float = 0.5) -> float:
    """D50 — dominance fraction: clones covering ``fraction`` of reads, over ``Sobs``.

    Ranks clonotypes by descending frequency and finds the minimum number ``k`` whose
    cumulative frequency reaches ``fraction``; returns ``k / Sobs``.

    Note:
        Legacy ``getDxxIndex`` returns the complement ``1 - k/Sobs``. This function
        follows the task's explicit definition (the covering fraction itself), which
        is the more common reading of "D50".

    Args:
        counts: Per-clonotype count vector.
        fraction: Cumulative-frequency target in ``[0, 1]`` (default ``0.5``).

    Returns:
        ``k / Sobs`` in ``(0, 1]``.
    """
    sobs = observed_richness(counts)
    if sobs == 0:
        return 0.0
    n = counts.sum()
    order = np.sort(counts)[::-1]
    cum = np.cumsum(order) / n
    k = int(np.searchsorted(cum, fraction, side="left")) + 1
    k = min(k, sobs)
    return k / sobs


def diversity_stats(df: pl.DataFrame, extrapolate_to: int | None = None) -> pl.DataFrame:
    """Compute all diversity estimators for a clonotype frame as a one-row frame.

    Args:
        df: A clonotype frame (one row per clonotype) with ``duplicate_count``.
        extrapolate_to: Target depth for ``chaoE`` (see :func:`chao_e`; defaults
            to ``2*n``).

    Returns:
        A single-row ``pl.DataFrame`` with columns ``reads, observed_diversity,
        chao1, chaoE, efron_thisted, shannon_wiener, normalized_shannon_wiener,
        inverse_simpson, d50``.
    """
    counts = _counts(df)
    return pl.DataFrame({
        "reads": [int(counts.sum())],
        "observed_diversity": [observed_richness(counts)],
        "chao1": [chao1(counts)],
        "chaoE": [chao_e(counts, extrapolate_to)],
        "efron_thisted": [efron_thisted(counts)],
        "shannon_wiener": [shannon_wiener(counts)],
        "normalized_shannon_wiener": [normalized_shannon_wiener(counts)],
        "inverse_simpson": [inverse_simpson(counts)],
        "d50": [d50(counts)],
    })
