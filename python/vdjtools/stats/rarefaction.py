"""Rarefaction / extrapolation of repertoire diversity — the iNEXT framework.

:func:`rarefaction` is the single canonical entry point for
rarefaction/extrapolation (R/E) of clonotype diversity. Its default
(``q=0, base="size"``) is the **vdjtools-original richness R/E curve** — the
classic clonotype-accumulation curve, computed by the validated iNEXT ``q = 0``
estimator. The same call also gives the full Hill-number profile (``q``) and
coverage-based R/E (``base="coverage"``).

The estimators themselves live in :mod:`vdjtools.stats.inext` (transcribed from
the iNEXT R package and numerically validated against it). :func:`inext` and
:func:`inext_coverage` are kept as aliases of the size- and coverage-based
engines for discoverability; :func:`rarefaction` dispatches to them.
"""
from __future__ import annotations

import polars as pl

from .inext import (  # noqa: F401  re-exported iNEXT public API
    asymptotic_diversity,
    coverage,
    estimate_d,
    inext,
    inext_batch,
    inext_coverage,
    rarefaction_batch,
    sample_coverage,
)


def rarefaction(data, q=0, base="size", *, sizes=None, coverages=None,
                endpoint=None, knots=40, se=True, nboot=50, conf=0.95,
                seed=0) -> pl.DataFrame:
    """Rarefaction/extrapolation of Hill-number diversity (canonical entry point).

    Single entry point for repertoire R/E, dispatching to the size- or
    coverage-based iNEXT engine (Chao et al. 2014; Hsieh et al. 2016). The
    default ``q=0, base="size"`` is the **vdjtools-original richness rarefaction /
    extrapolation curve** — the classic clonotype-accumulation curve, computed by
    the validated iNEXT ``q = 0`` estimator. Pass a tuple of orders for the full
    Hill-number profile, or ``base="coverage"`` for coverage-based R/E.

    Args:
        data: A 1-D count vector (list/``np.ndarray``/``pl.Series``) of clonotype
            abundances, or a clonotype ``pl.DataFrame`` (``duplicate_count`` used).
        q: Hill order(s) — a scalar (default ``0`` = richness, the vdjtools-original
            curve) or a tuple such as ``(0, 1, 2)`` for the richness/Shannon/Simpson
            profile.
        base: ``"size"`` for size-based R/E (uses ``sizes``/``endpoint``/``knots``)
            or ``"coverage"`` for coverage-based R/E (uses ``coverages``/``knots``).
        sizes: Explicit sampling depths for ``base="size"``; see :func:`inext`.
        coverages: Explicit target coverages for ``base="coverage"``; see
            :func:`inext_coverage`.
        endpoint: Maximum depth for the default size grid (``base="size"``).
        knots: Number of grid points when ``sizes``/``coverages`` are ``None``.
        se: If ``True``, compute bootstrap confidence intervals.
        nboot: Number of bootstrap replicates.
        conf: Confidence level for the intervals.
        seed: Seed for the bootstrap RNG.

    Returns:
        A tidy ``pl.DataFrame``. For ``base="size"`` (see :func:`inext`): columns
        ``order_q``, ``m``, ``method``, ``sample_coverage``, ``qD``, ``qD_lo``,
        ``qD_hi``. For ``base="coverage"`` (see :func:`inext_coverage`): columns
        ``order_q``, ``sample_coverage``, ``m``, ``method``, ``qD``, ``qD_lo``,
        ``qD_hi``.

    Raises:
        ValueError: If ``base`` is not ``"size"`` or ``"coverage"``.
    """
    if base == "size":
        return inext(data, q, sizes=sizes, endpoint=endpoint, knots=knots,
                     se=se, nboot=nboot, conf=conf, seed=seed)
    if base == "coverage":
        return inext_coverage(data, q, coverages=coverages, knots=knots,
                              se=se, nboot=nboot, conf=conf, seed=seed)
    raise ValueError(f"base must be 'size' or 'coverage', got {base!r}")
