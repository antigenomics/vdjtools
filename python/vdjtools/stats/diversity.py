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
    """Extract the positive integer clonotype count vector from a frame (sorted).

    Returned ascending so every estimator's floating-point reduction (Σf·ln f, Σf²) is
    summed in a fixed, numerically-stable order — making :func:`diversity_stats`
    independent of the input's clonotype row order, so the streamed cohort path
    (:func:`diversity_cohort`, which reconstructs the vector from the count spectrum)
    is bit-identical to the per-sample path.
    """
    c = df[COUNT].drop_nulls().to_numpy().astype(np.int64)
    return np.sort(c[c > 0])


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
    truncation depth (up to ``max_depth``) until the coefficient of variation
    ``D/S`` reaches ``cv_threshold``, accumulating alternating binomial coefficients
    over the count-frequency spectrum ``f_x`` (number of clonotypes seen exactly
    ``x`` times). The CV stopping rule is the *only* stopping rule — legacy has no
    "stop at the max observed count" cap.

    Args:
        counts: Per-clonotype count vector.
        max_depth: Maximum truncation depth.
        cv_threshold: Stop once ``std/estimate`` reaches this value.

    Returns:
        The Efron–Thisted richness estimate.
    """
    sobs = observed_richness(counts)
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
    """D50 — the legacy ``getDxxIndex`` dominance index ``1 - k/Sobs``.

    Ranks clonotypes by descending count and finds the minimum number ``k`` whose
    cumulative read fraction reaches ``fraction``; returns ``1 - k/Sobs`` — the
    fraction of clonotypes *not* needed to cover ``fraction`` of the reads. This is
    the exact legacy ``ExactEstimator.getDxxIndex`` definition
    (``1.0 - div / frequencyTable.diversity``).

    Args:
        counts: Per-clonotype count vector.
        fraction: Cumulative-frequency target in ``[0, 1]`` (default ``0.5``).

    Returns:
        ``1 - k/Sobs`` in ``[0, 1)``; ``0.0`` when ``Sobs == 0``.
    """
    sobs = observed_richness(counts)
    if sobs == 0:
        return 0.0
    n = counts.sum()
    order = np.sort(counts)[::-1]
    cum = np.cumsum(order) / n
    k = int(np.searchsorted(cum, fraction, side="left")) + 1
    k = min(k, sobs)
    return 1.0 - k / sobs


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


def diversity_cohort(cohort, extrapolate_to: int | None = None, *,
                     sample_col: str = "sample_id") -> pl.DataFrame:
    """Per-sample :func:`diversity_stats` for a whole cohort, low peak memory + exact.

    Collapses the cohort to each sample's **count-frequency spectrum** in one streamed
    pass — ``group_by([sample_id, duplicate_count]).agg(len)`` — which is bounded by the
    number of distinct counts per sample (tiny), never the number of clonotypes. Each
    sample's count vector is then reconstructed from its spectrum and fed to the exact
    per-sample estimators, so every value is **bit-identical** to running
    :func:`diversity_stats` on the full sample (including Efron–Thisted / chaoE / d50,
    which need the spectrum/sorted form). Peak memory is ``O(Σ distinct counts)``, so a
    cohort far larger than RAM is handled with one :meth:`~polars.LazyFrame.collect`.

    Args:
        cohort: A clonotype cohort with a ``sample_id`` column — a lazy
            :func:`vdjtools.io.scan_cohort` frame (streamed) or an eager
            ``pl.DataFrame``.
        extrapolate_to: Target depth for ``chaoE`` (see :func:`diversity_stats`).
        sample_col: The per-sample id column (default ``"sample_id"``).

    Returns:
        One row per sample: ``sample_id`` followed by the :func:`diversity_stats`
        columns, in first-appearance order of ``sample_id``.
    """
    lf = cohort if isinstance(cohort, pl.LazyFrame) else cohort.lazy()
    # Sort the (tiny) spectrum by sample_col so partition order — and thus output row
    # order — is deterministic (sample_id-sorted), independent of scan/group_by order.
    spectrum = (lf.filter(pl.col(COUNT) > 0)
                  .group_by([sample_col, COUNT]).agg(pl.len().alias("_fx"))
                  .collect(engine="streaming").sort(sample_col))
    rows = []
    for sub in spectrum.partition_by(sample_col, maintain_order=True):
        counts = np.repeat(sub[COUNT].to_numpy(), sub["_fx"].to_numpy())
        row = diversity_stats(pl.DataFrame({COUNT: counts}), extrapolate_to)
        rows.append(row.select(pl.lit(sub[sample_col][0]).alias(sample_col), pl.all()))
    if not rows:
        base = diversity_stats(pl.DataFrame({COUNT: np.array([], dtype=np.int64)}))
        return base.select(pl.lit(None, dtype=pl.Utf8).alias(sample_col), pl.all()).clear()
    return pl.concat(rows, how="vertical_relaxed")
