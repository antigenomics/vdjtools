"""Rarefaction: analytic interpolation + Chao extrapolation of clonotype richness.

Interpolation is the exact Hurlbert/Coleman expectation (no random draws) — the
same estimator legacy vdjtools computes in ``ChaoEstimator.chaoI`` — and
extrapolation follows ``ChaoEstimator.chaoE``. Both are deterministic.

:func:`rarefaction` is the legacy richness-only (``q = 0``) curve. For the full
iNEXT framework — size- and coverage-based R/E of Hill numbers of orders 0/1/2
with bootstrap confidence intervals — use :func:`inext`, :func:`inext_coverage`,
:func:`asymptotic_diversity`, :func:`sample_coverage`, and :func:`estimate_d`
(implemented in :mod:`vdjtools.stats.inext` and re-exported here).
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.special import gammaln

from ..io.schema import COUNT
from .diversity import observed_richness
from .inext import (  # noqa: F401  re-exported iNEXT public API
    asymptotic_diversity,
    estimate_d,
    inext,
    inext_coverage,
    sample_coverage,
)


def _spectrum(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return distinct counts ``k`` and their multiplicities ``f_k``."""
    k, fk = np.unique(counts, return_counts=True)
    return k.astype(np.int64), fk.astype(np.float64)


def _f0(counts: np.ndarray) -> float:
    """Chao ``F0`` = estimated number of unseen clonotypes."""
    f1 = float(np.count_nonzero(counts == 1))
    f2 = float(np.count_nonzero(counts == 2))
    return f1 * (f1 - 1) / (2 * (f2 + 1))


def _interpolate(k: np.ndarray, fk: np.ndarray, n: int, sobs: int, f0: float,
                 m: int) -> tuple[float, float]:
    """Hurlbert expected richness (and its std) at subsample size ``m`` (0<=m<=n)."""
    if m <= 0:
        return 0.0, 0.0
    if m >= n:
        return float(sobs), 0.0
    seen = k <= (n - m)  # clonotypes that can be missed at depth m
    # alpha = C(n-k, m) / C(n, m) — probability clonotype with count k is unseen.
    log_denom = gammaln(n + 1) - gammaln(n - m + 1)
    alpha = np.zeros_like(fk)
    ks = k[seen]
    alpha[seen] = np.exp(gammaln(n - ks + 1) - gammaln(n - ks - m + 1) - log_denom)
    sum1 = float(np.sum(fk[seen] * alpha[seen]))
    sum2 = float(np.sum(fk[seen] * (1 - alpha[seen]) ** 2) + np.sum(fk[~seen]))
    s_ind = sobs - sum1
    var = sum2 - (s_ind * s_ind) / (sobs + f0) if (sobs + f0) > 0 else sum2
    return s_ind, float(np.sqrt(max(var, 0.0)))


def _extrapolate(counts: np.ndarray, n: int, sobs: int, f0: float,
                 target: int) -> tuple[float, float]:
    """Chao extrapolated richness (and delta-method std) at depth ``target`` (>n)."""
    f1 = float(np.count_nonzero(counts == 1))
    f2 = float(np.count_nonzero(counts == 2))
    if f0 == 0 or n == 0:
        return float(sobs), 0.0
    m_star = target - n
    r = 1.0 - f1 / (n * f0)
    brackets = 1.0 - r ** m_star
    d_f0_d_f1 = (2 * f1 - 1) / (2 * (f2 + 1))
    d_f0_d_f2 = -f1 * (f1 - 1) / (2 * (f2 + 1) ** 2)
    d_brackets = -m_star * r ** (m_star - 1)
    d_brackets_d_f1 = d_brackets * (-1 / (n * f0) + f1 / (n * f0 * f0) * d_f0_d_f1)
    d_brackets_d_f2 = d_brackets * (f1 / (n * f0 * f0) * d_f0_d_f2)
    d_s_d_f1 = d_f0_d_f1 * brackets + f0 * d_brackets_d_f1
    d_s_d_f2 = d_f0_d_f2 * brackets + f0 * d_brackets_d_f2
    denom = sobs + f0
    cov11 = f1 * (1 - f1 / denom)
    cov22 = f2 * (1 - f2 / denom)
    cov12 = -f1 * f2 / denom
    var = (sobs * (1 - sobs / denom) + d_s_d_f1 ** 2 * cov11
           + d_s_d_f2 ** 2 * cov22 + 2 * d_s_d_f1 * d_s_d_f2 * cov12)
    return sobs + f0 * brackets, float(np.sqrt(max(var, 0.0)))


def rarefaction(df: pl.DataFrame, steps: int = 40,
                extrapolate_to: int | None = None) -> pl.DataFrame:
    """Rarefaction curve: interpolated, observed, and extrapolated richness.

    Builds ``steps`` sampling depths from 0 to ``extrapolate_to``. Points below the
    observed depth ``n`` use the exact Hurlbert interpolation; the point at ``n`` is
    the observed richness; points above ``n`` use Chao extrapolation. Confidence
    bands are ``mean ± 1.96 * std``.

    Args:
        df: A clonotype frame (one row per clonotype) with ``duplicate_count``.
        steps: Number of depths sampled across ``[0, extrapolate_to]``.
        extrapolate_to: Maximum depth. Defaults to ``2*n``.

    Returns:
        A ``pl.DataFrame`` with columns ``x`` (Int64 depth), ``mean``, ``ci_lo``,
        ``ci_hi`` (Float64), and ``kind`` (``interpolated`` | ``observed`` |
        ``extrapolated``), ordered by ``x``.
    """
    counts = df[COUNT].drop_nulls().to_numpy().astype(np.int64)
    counts = counts[counts > 0]
    n = int(counts.sum())
    sobs = observed_richness(counts)
    f0 = _f0(counts)
    k, fk = _spectrum(counts)
    target = 2 * n if extrapolate_to is None else int(extrapolate_to)
    if target < n:
        raise ValueError(f"extrapolate_to ({target}) must be >= n ({n})")

    xs = np.unique(np.concatenate([
        np.linspace(0, target, steps), np.array([n]),
    ]).astype(np.int64))

    rows = []
    for x in xs:
        x = int(x)
        if x < n:
            mean, std = _interpolate(k, fk, n, sobs, f0, x)
            kind = "interpolated"
        elif x == n:
            mean, std, kind = float(sobs), 0.0, "observed"
        else:
            mean, std = _extrapolate(counts, n, sobs, f0, x)
            kind = "extrapolated"
        rows.append((x, mean, mean - 1.96 * std, mean + 1.96 * std, kind))

    return pl.DataFrame(
        rows, schema=["x", "mean", "ci_lo", "ci_hi", "kind"], orient="row",
    ).with_columns(pl.col("x").cast(pl.Int64))
