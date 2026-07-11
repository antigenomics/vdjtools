"""iNEXT interpolation/extrapolation of Hill-number diversity for a repertoire.

Size- and coverage-based rarefaction/extrapolation (R/E) plus asymptotic
diversity estimators for Hill numbers of order ``q`` (0 = richness, 1 = Shannon,
2 = Simpson). The estimators are transcribed from the ``iNEXT`` R package
(v3.0.2 internals ``TD.m.est``, ``Chat.Ind``, ``Diversity_profile``,
``EstiBootComm.Ind``, ``invChat.Ind``) and the underlying papers, and are
numerically validated against that package in ``tests/python/test_inext.py``.

For extrapolation of orders ``q != 2`` this uses iNEXT's "beta" method
(``qD(n+m*) = D_obs + (D_asy - D_obs) * (1 - (1-beta)^m*)``), which matches the R
package rather than Eq. 10c of Chao et al. (2014). Order ``q = 2`` uses the exact
closed form for both interpolation and extrapolation. Only non-negative integer
orders are supported (iNEXT's default profile is ``q = (0, 1, 2)``).

References:
    Chao, A., Gotelli, N. J., Hsieh, T. C., Sander, E. L., Ma, K. H., Colwell,
    R. K., & Ellison, A. M. (2014). Rarefaction and extrapolation with Hill
    numbers: a framework for sampling and estimation in species diversity
    studies. *Ecological Monographs*, 84(1), 45-67. doi:10.1890/13-0133.1

    Hsieh, T. C., Ma, K. H., & Chao, A. (2016). iNEXT: an R package for
    rarefaction and extrapolation of species diversity (Hill numbers). *Methods
    in Ecology and Evolution*, 7(12), 1451-1456. doi:10.1111/2041-210X.12613
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import polars as pl
from scipy.optimize import brentq
from scipy.special import digamma, gammaln
from scipy.stats import norm

from ..io.schema import COUNT

# Output schema for the size-based / estimate_d rows.
_RE_SCHEMA = {
    "order_q": pl.Int64,
    "m": pl.Int64,
    "method": pl.Utf8,
    "sample_coverage": pl.Float64,
    "qD": pl.Float64,
    "qD_lo": pl.Float64,
    "qD_hi": pl.Float64,
}


# --------------------------------------------------------------------------- #
# input coercion
# --------------------------------------------------------------------------- #
def _as_counts(data) -> np.ndarray:
    """Coerce ``data`` to a positive ``float64`` clonotype-count vector."""
    if isinstance(data, pl.DataFrame):
        arr = data[COUNT].drop_nulls().to_numpy()
    elif isinstance(data, pl.Series):
        arr = data.drop_nulls().to_numpy()
    else:
        arr = np.asarray(data)
    arr = arr.astype(np.float64)
    arr = arr[arr > 0]
    if arr.size == 0:
        raise ValueError("iNEXT requires a non-empty vector of positive counts")
    return arr


def _as_orders(q) -> list[int]:
    """Coerce ``q`` to a list of non-negative integer Hill orders."""
    seq = [q] if isinstance(q, (int, float)) else list(q)
    out: list[int] = []
    for v in seq:
        if v < 0 or v != int(v):
            raise ValueError(f"Hill orders must be non-negative integers, got {v!r}")
        out.append(int(v))
    return out


def _spectrum(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Distinct counts and their multiplicities ``f_k`` (both ``float64``)."""
    c, f = np.unique(x, return_counts=True)
    return c.astype(np.float64), f.astype(np.float64)


def _logbinom(a, b):
    """Elementwise ``log C(a, b)`` (``-inf`` when ``b < 0`` or ``a < b``).

    Uses ``gammaln`` so ``a`` may be non-integer; broadcasts and works for
    scalars.
    """
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    valid = (b >= -1e-12) & (a - b >= -1e-12)
    aa = np.where(valid, a, 0.0)
    bb = np.where(valid, b, 0.0)
    val = gammaln(aa + 1) - gammaln(bb + 1) - gammaln(aa - bb + 1)
    return np.where(valid, val, -np.inf)


# --------------------------------------------------------------------------- #
# coverage (Chat.Ind)
# --------------------------------------------------------------------------- #
def _chat(x: np.ndarray, m: float) -> float:
    """Sample coverage estimate ``Ĉ(m)`` at (possibly non-integer) size ``m``.

    Chao et al. (2014) Eq. (4a-4c): interpolation for ``m < n`` (linear blend of
    floor/ceil for non-integer ``m``), the observed coverage at ``m = n``, and
    ``1 - (f1/n) * A^(m-n+1)`` for extrapolation.
    """
    n = x.sum()
    f1 = float(np.count_nonzero(x == 1))
    f2 = float(np.count_nonzero(x == 2))
    if f2 == 0:
        f0 = (n - 1) / n * f1 * (f1 - 1) / 2
    else:
        f0 = (n - 1) / n * f1**2 / 2 / f2
    a_ext = n * f0 / (n * f0 + f1) if f1 > 0 else 1.0

    def at_int(k: float) -> float:
        if k == n:
            return 1 - f1 / n * a_ext
        xx = x[(n - x) >= k]
        return 1 - float(np.sum(
            xx / n * np.exp(gammaln(n - xx + 1) - gammaln(n - xx - k + 1)
                            - gammaln(n) + gammaln(n - k))))

    if m == n:
        return 1 - f1 / n * a_ext
    if m > n:
        return 1 - f1 / n * a_ext ** (m - n + 1)
    if m == round(m):
        return at_int(float(round(m)))
    lo, hi = np.floor(m), np.ceil(m)
    return (hi - m) * at_int(lo) + (m - lo) * at_int(hi)


# --------------------------------------------------------------------------- #
# diversity: plug-in, asymptotic, interpolation, extrapolation
# --------------------------------------------------------------------------- #
def _plugin(x: np.ndarray, q: int) -> float:
    """Empirical (MLE plug-in) Hill number of order ``q`` on ``X_i/n``."""
    p = x / x.sum()
    if q == 0:
        return float(np.count_nonzero(p > 0))
    if q == 1:
        return float(np.exp(-np.sum(p * np.log(p))))
    return float(np.exp(np.log(np.sum(p**q)) / (1 - q)))


def _asymptotic(x: np.ndarray, q: int) -> float:
    """Asymptotic Hill-number estimator of order ``q`` (``iNEXT`` Diversity_profile).

    - ``q = 0``: Chao1 richness ``S_obs + f0``.
    - ``q = 1``: Chao-Wang-Jost Shannon Hill number ``exp(H_hat)``.
    - integer ``q >= 2``: MVUE estimator ``(Σ_{X_i>=q} C(X_i,q)/C(n,q))^(1/(1-q))``
      (``q = 2`` gives ``n(n-1)/Σ X_i(X_i-1)``).
    """
    n = x.sum()
    f1 = float(np.count_nonzero(x == 1))
    f2 = float(np.count_nonzero(x == 2))
    if q == 0:
        f0 = f1**2 / 2 / f2 if f2 > 0 else f1 * (f1 - 1) / 2
        return float(len(x) + (n - 1) / n * f0)
    if q == 1:
        c, f = _spectrum(x)
        first = float(np.sum(f * c / n * (digamma(n) - digamma(c))))
        if f2 > 0:
            a_cwj = 2 * f2 / ((n - 1) * f1 + 2 * f2)
        elif f1 > 0:
            a_cwj = 2 / ((n - 1) * (f1 - 1) + 2)
        else:
            a_cwj = 1.0
        if f1 == 0 or a_cwj == 1.0:
            second = 0.0
        else:
            r = np.arange(1, int(n))
            second = (f1 / n * (1 - a_cwj) ** (-n + 1)
                      * (-np.log(a_cwj) - np.sum(1 / r * (1 - a_cwj) ** r)))
        return float(np.exp(first + second))
    c, f = _spectrum(x)
    mask = c >= q
    moment = float(np.sum(f[mask] * np.exp(_logbinom(c[mask], q) - _logbinom(n, q))))
    return float(moment ** (1 / (1 - q)))


def _rtd_moment(cs: np.ndarray, fs: np.ndarray, n: float, m: int, gfun) -> float:
    """``Σ_k g(k) * fhat_k(m)`` via the MVUE frequency counts at integer ``m < n``.

    Computed by summing over the abundance *spectrum* rather than over all
    ``(k, species)`` pairs: for each distinct count ``c`` the inner ``k`` range is
    ``1..min(c, m)``, so the total work is ``O(Σ_distinct min(c, m))`` — bounded by
    ``O(n)`` per size, keeping large repertoires tractable (rather than ``O(m·S)``).
    """
    total = 0.0
    log_cnm = _logbinom(n, m)
    for cval, fcount in zip(cs, fs):
        kmax = int(min(cval, m))
        k = np.arange(1, kmax + 1)
        logterm = _logbinom(cval, k) + _logbinom(n - cval, m - k) - log_cnm
        total += fcount * float(np.sum(gfun(k) * np.exp(logterm)))
    return total


def _d2(x: np.ndarray, m: float) -> float:
    """Closed-form order-2 Hill number at size ``m`` (interpolation and extrapolation)."""
    n = x.sum()
    s = float(np.sum(x * (x - 1)) / (n * (n - 1)))
    return float(1.0 / (1.0 / m + (m - 1.0) / m * s))


def _rtd(x: np.ndarray, m: int, q: int) -> float:
    """Rarefied (interpolated) Hill number of order ``q`` at integer size ``m``."""
    cs, fs = _spectrum(x)
    n = x.sum()
    if q == 0:
        return _rtd_moment(cs, fs, n, m, lambda k: np.ones_like(k, dtype=np.float64))
    if q == 1:
        s = _rtd_moment(cs, fs, n, m, lambda k: -(k / m) * np.log(k / m))
        return float(np.exp(s))
    if q == 2:
        return _d2(x, m)
    s = _rtd_moment(cs, fs, n, m, lambda k: (k / m) ** q)
    return float(s ** (1 / (1 - q)))


def _diversity_at(x: np.ndarray, m: float, q: int) -> float:
    """Hill number of order ``q`` at size ``m`` (``iNEXT`` ``TD.m.est``).

    ``m < n`` interpolates (MVUE, with a linear floor/ceil blend for non-integer
    ``m``); ``m == n`` returns the plug-in value; ``m > n`` extrapolates via the
    beta method (``q != 2``) or the closed form (``q == 2``).
    """
    n = x.sum()
    if q == 2:
        return _plugin(x, q) if m == n else _d2(x, m)
    obs = _plugin(x, q)
    if m < n:
        if m == round(m):
            return _rtd(x, int(round(m)), q)
        lo, hi = int(np.floor(m)), int(np.ceil(m))
        return (hi - m) * _rtd(x, lo, q) + (m - lo) * _rtd(x, hi, q)
    if m == n:
        return obs
    asy = _asymptotic(x, q)
    rfd = _rtd(x, int(n - 1), q)
    beta = 0.0 if asy == rfd else (obs - rfd) / (asy - rfd)
    mstar = m - n
    return obs + (asy - obs) * (1 - (1 - beta) ** mstar)


def _invert_coverage(x: np.ndarray, cvrg: float) -> float:
    """Size ``m`` achieving target coverage ``cvrg`` (``iNEXT`` ``invChat.Ind``)."""
    n = x.sum()
    ref = _chat(x, n)
    if abs(cvrg - ref) <= 1e-12:
        return float(n)
    if cvrg < ref:
        return float(brentq(lambda m: _chat(x, m) - cvrg, 1e-9, n))
    f1 = float(np.count_nonzero(x == 1))
    f2 = float(np.count_nonzero(x == 2))
    if f1 > 0 and f2 > 0:
        a_ext = (n - 1) * f1 / ((n - 1) * f1 + 2 * f2)
    elif f1 > 1 and f2 == 0:
        a_ext = (n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2)
    else:
        a_ext = 0.0
    mstar = 0.0 if a_ext == 0 else (np.log(n / f1) + np.log(1 - cvrg)) / np.log(a_ext) - 1
    return float(max(n + mstar, 1.0))


# --------------------------------------------------------------------------- #
# bootstrap (EstiBootComm.Ind)
# --------------------------------------------------------------------------- #
def _bootstrap_probs(x: np.ndarray) -> np.ndarray:
    """Augmented-assemblage detection probabilities (``iNEXT`` ``EstiBootComm.Ind``)."""
    n = x.sum()
    f1 = float(np.count_nonzero(x == 1))
    f2 = float(np.count_nonzero(x == 2))
    if f2 == 0:
        f0 = (n - 1) / n * f1 * (f1 - 1) / 2
    else:
        f0 = (n - 1) / n * f1**2 / 2 / f2
    a_ext = n * f0 / (n * f0 + f1) if f1 > 0 else 1.0
    a = f1 / n * a_ext
    b = float(np.sum(x / n * (1 - x / n) ** n))
    w = 0.0 if (f0 == 0 or b == 0) else a / b
    p_obs = x / n * (1 - w * (1 - x / n) ** n)
    k = int(np.ceil(f0))
    if k > 0:
        return np.concatenate([p_obs, np.full(k, a / k)])
    return p_obs


def _bootstrap_se(x: np.ndarray, ms, qs, nboot: int, seed: int) -> np.ndarray:
    """Bootstrap standard errors of ``qD(m)`` — shape ``(len(qs), len(ms))``."""
    p = _bootstrap_probs(x)
    p = p / p.sum()
    n = int(x.sum())
    rng = np.random.default_rng(seed)
    reps = rng.multinomial(n, p, size=nboot)
    vals = np.empty((nboot, len(qs), len(ms)))
    for r in range(nboot):
        xb = reps[r].astype(np.float64)
        xb = xb[xb > 0]
        for i, q in enumerate(qs):
            for j, m in enumerate(ms):
                vals[r, i, j] = _diversity_at(xb, float(m), q)
    return vals.std(axis=0, ddof=1)


def _z(conf: float) -> float:
    return float(norm.ppf(1 - (1 - conf) / 2))


# --------------------------------------------------------------------------- #
# public API
# --------------------------------------------------------------------------- #
def inext(data, q=(0, 1, 2), *, sizes=None, endpoint=None, knots=40, se=True,
          nboot=50, conf=0.95, seed=0) -> pl.DataFrame:
    """Size-based rarefaction/extrapolation of Hill-number diversity.

    Interpolates (``m < n``) and extrapolates (``m > n``) the Hill number of each
    order in ``q`` across a grid of sampling depths, following Chao et al. (2014)
    and the ``iNEXT`` R package (Hsieh et al. 2016). Point estimates are computed
    on the original data; confidence bands come from an augmented-assemblage
    bootstrap.

    Args:
        data: A 1-D count vector (list/``np.ndarray``/``pl.Series``) of clonotype
            abundances, or a clonotype ``pl.DataFrame`` (``duplicate_count`` used).
        q: Hill orders (non-negative integers); default ``(0, 1, 2)``.
        sizes: Explicit sampling depths. If ``None``, ``knots`` depths from 1 to
            ``endpoint`` are used, always including the observed depth ``n``.
        endpoint: Maximum depth when ``sizes`` is ``None``. Defaults to ``2*n``.
        knots: Number of depths when ``sizes`` is ``None``.
        se: If ``True``, compute bootstrap confidence intervals.
        nboot: Number of bootstrap replicates.
        conf: Confidence level for the intervals.
        seed: Seed for the bootstrap RNG.

    Returns:
        A tidy ``pl.DataFrame`` with columns ``order_q`` (Int64), ``m`` (Int64),
        ``method`` (``rarefaction`` | ``observed`` | ``extrapolation``),
        ``sample_coverage``, ``qD``, ``qD_lo``, ``qD_hi`` (the CI columns are null
        when ``se=False``).
    """
    x = _as_counts(data)
    n = int(x.sum())
    qs = _as_orders(q)
    if sizes is None:
        end = 2 * n if endpoint is None else int(endpoint)
        grid = np.floor(np.linspace(1, end, knots)).astype(np.int64)
        grid = np.unique(np.concatenate([grid, np.array([n], dtype=np.int64)]))
    else:
        grid = np.unique(np.asarray(sizes, dtype=np.int64))
    ms = [int(m) for m in grid if m >= 1]

    se_arr = _bootstrap_se(x, ms, qs, nboot, seed) if se else None
    z = _z(conf) if se else 0.0
    rows = []
    for i, qq in enumerate(qs):
        for j, m in enumerate(ms):
            qd = _diversity_at(x, float(m), qq)
            cov = min(max(_chat(x, float(m)), 0.0), 1.0)
            method = ("observed" if m == n
                      else "rarefaction" if m < n else "extrapolation")
            if se:
                s = se_arr[i, j]
                lo, hi = max(0.0, qd - z * s), qd + z * s
            else:
                lo = hi = None
            rows.append((qq, m, method, cov, qd, lo, hi))
    return pl.DataFrame(rows, schema=_RE_SCHEMA, orient="row")


def inext_coverage(data, q=(0, 1, 2), *, coverages=None, knots=40, se=True,
                   nboot=50, conf=0.95, seed=0) -> pl.DataFrame:
    """Coverage-based rarefaction/extrapolation of Hill-number diversity.

    Each target sample coverage is inverted to a sampling depth ``m`` (``iNEXT``
    ``invChat.Ind``) and the Hill number is then evaluated at that depth.

    Args:
        data: Count vector or clonotype ``pl.DataFrame`` (see :func:`inext`).
        q: Hill orders (non-negative integers); default ``(0, 1, 2)``.
        coverages: Explicit target coverages in ``(0, 1)``. If ``None``, a grid of
            the coverages attained by the default size grid (1..``2*n``) is used.
        knots: Number of coverages when ``coverages`` is ``None``.
        se: If ``True``, compute bootstrap confidence intervals.
        nboot: Number of bootstrap replicates.
        conf: Confidence level.
        seed: Seed for the bootstrap RNG.

    Returns:
        A tidy ``pl.DataFrame`` with columns ``order_q``, ``sample_coverage``,
        ``m`` (Float64), ``method``, ``qD``, ``qD_lo``, ``qD_hi``.
    """
    x = _as_counts(data)
    n = int(x.sum())
    qs = _as_orders(q)
    if coverages is None:
        end = 2 * n
        grid = np.unique(np.concatenate([
            np.floor(np.linspace(1, end, knots)).astype(np.int64),
            np.array([n], dtype=np.int64)]))
        covs = np.array([_chat(x, float(m)) for m in grid if m >= 1])
        cvals = np.unique(np.clip(covs, 0.0, 1.0))
    else:
        cvals = np.asarray(list(coverages) if isinstance(coverages, Iterable)
                           else [coverages], dtype=np.float64)

    ms = [_invert_coverage(x, float(c)) for c in cvals]
    se_arr = _bootstrap_se(x, ms, qs, nboot, seed) if se else None
    z = _z(conf) if se else 0.0
    schema = {
        "order_q": pl.Int64, "sample_coverage": pl.Float64, "m": pl.Float64,
        "method": pl.Utf8, "qD": pl.Float64, "qD_lo": pl.Float64, "qD_hi": pl.Float64,
    }
    rows = []
    for i, qq in enumerate(qs):
        for j, (c, m) in enumerate(zip(cvals, ms)):
            qd = _diversity_at(x, m, qq)
            method = ("observed" if m == n
                      else "rarefaction" if m < n else "extrapolation")
            if se:
                s = se_arr[i, j]
                lo, hi = max(0.0, qd - z * s), qd + z * s
            else:
                lo = hi = None
            rows.append((qq, float(c), float(m), method, qd, lo, hi))
    return pl.DataFrame(rows, schema=schema, orient="row")


def asymptotic_diversity(data, q=(0, 1, 2), *, se=True, nboot=50, conf=0.95,
                         seed=0) -> pl.DataFrame:
    """Asymptotic (estimated true) Hill-number diversity of each order.

    Uses Chao1 (``q = 0``), the Chao-Wang-Jost Shannon Hill number (``q = 1``) and
    the MVUE Simpson Hill number (``q = 2``); see :func:`_asymptotic`.

    Args:
        data: Count vector or clonotype ``pl.DataFrame`` (see :func:`inext`).
        q: Hill orders (non-negative integers); default ``(0, 1, 2)``.
        se: If ``True``, compute bootstrap standard errors and intervals.
        nboot: Number of bootstrap replicates.
        conf: Confidence level.
        seed: Seed for the bootstrap RNG.

    Returns:
        A ``pl.DataFrame`` with columns ``order_q``, ``observed`` (plug-in),
        ``estimator`` (asymptotic), ``se``, ``lo``, ``hi``.
    """
    x = _as_counts(data)
    qs = _as_orders(q)
    z = _z(conf)
    reps = None
    if se:
        p = _bootstrap_probs(x)
        p = p / p.sum()
        rng = np.random.default_rng(seed)
        reps = rng.multinomial(int(x.sum()), p, size=nboot)
    schema = {
        "order_q": pl.Int64, "observed": pl.Float64, "estimator": pl.Float64,
        "se": pl.Float64, "lo": pl.Float64, "hi": pl.Float64,
    }
    rows = []
    for qq in qs:
        obs = _plugin(x, qq)
        est = _asymptotic(x, qq)
        if se:
            vals = np.array([
                _asymptotic(xb[xb > 0], qq)
                for xb in reps.astype(np.float64)])
            s = float(vals.std(ddof=1))
            lo, hi = max(0.0, est - z * s), est + z * s
        else:
            s = lo = hi = None
        rows.append((qq, obs, est, s, lo, hi))
    return pl.DataFrame(rows, schema=schema, orient="row")


def sample_coverage(data, m=None):
    """Estimated sample coverage ``Ĉ`` (Chao et al. 2014).

    Args:
        data: Count vector or clonotype ``pl.DataFrame`` (see :func:`inext`).
        m: A single depth, an iterable of depths, or ``None`` (the observed depth
            ``n``).

    Returns:
        ``Ĉ(n)`` (float) when ``m`` is ``None``, ``Ĉ(m)`` (float) for a scalar
        ``m``, or a ``pl.DataFrame`` with columns ``m`` and ``sample_coverage``
        for an iterable ``m``.
    """
    x = _as_counts(data)
    n = int(x.sum())
    if m is None:
        return min(max(_chat(x, float(n)), 0.0), 1.0)
    if isinstance(m, Iterable):
        ms = [float(v) for v in m]
        return pl.DataFrame({
            "m": [int(v) for v in ms],
            "sample_coverage": [min(max(_chat(x, v), 0.0), 1.0) for v in ms],
        })
    return min(max(_chat(x, float(m)), 0.0), 1.0)


def estimate_d(data, base="size", level=None, q=(0, 1, 2), *, se=True, nboot=50,
               conf=0.95, seed=0) -> pl.DataFrame:
    """Diversity at a target sample size or coverage (``iNEXT`` ``estimateD`` analog).

    Args:
        data: Count vector or clonotype ``pl.DataFrame`` (see :func:`inext`).
        base: ``"size"`` to fix a sampling depth, ``"coverage"`` to fix a coverage.
        level: The depth (``base="size"``; defaults to the observed depth ``n``) or
            the coverage (``base="coverage"``; required, in ``(0, 1)``).
        q: Hill orders (non-negative integers); default ``(0, 1, 2)``.
        se: If ``True``, compute bootstrap confidence intervals.
        nboot: Number of bootstrap replicates.
        conf: Confidence level.
        seed: Seed for the bootstrap RNG.

    Returns:
        A ``pl.DataFrame`` with columns ``order_q``, ``m`` (Float64), ``method``,
        ``sample_coverage``, ``qD``, ``qD_lo``, ``qD_hi``.
    """
    x = _as_counts(data)
    n = int(x.sum())
    qs = _as_orders(q)
    if base == "size":
        m = float(n if level is None else level)
        cov = min(max(_chat(x, m), 0.0), 1.0)
    elif base == "coverage":
        if level is None:
            raise ValueError("estimate_d(base='coverage') requires an explicit level")
        cov = float(level)
        m = _invert_coverage(x, cov)
    else:
        raise ValueError(f"base must be 'size' or 'coverage', got {base!r}")

    se_arr = _bootstrap_se(x, [m], qs, nboot, seed) if se else None
    z = _z(conf) if se else 0.0
    method = "observed" if m == n else "rarefaction" if m < n else "extrapolation"
    schema = {
        "order_q": pl.Int64, "m": pl.Float64, "method": pl.Utf8,
        "sample_coverage": pl.Float64, "qD": pl.Float64,
        "qD_lo": pl.Float64, "qD_hi": pl.Float64,
    }
    rows = []
    for i, qq in enumerate(qs):
        qd = _diversity_at(x, m, qq)
        if se:
            s = se_arr[i, 0]
            lo, hi = max(0.0, qd - z * s), qd + z * s
        else:
            lo = hi = None
        rows.append((qq, float(m), method, cov, qd, lo, hi))
    return pl.DataFrame(rows, schema=schema, orient="row")
