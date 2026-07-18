"""edgeR-style negative-binomial exact test for clonotype expansion between two timepoints.

The complementary per-clone caller of Pavlova, Zvyagin & Shugay 2024 (§2.5): the classic
edgeR pipeline — **TMM** library normalization, a **common dispersion** by quantile-adjusted
conditional maximum likelihood (qCML), and the **negative-binomial exact test** — used to call
vaccine-associated clonotypes (``|log2 FC| >= 5`` & ``p <= 0.01`` in the paper).

For two libraries with one replicate each, the NB exact test collapses to a closed form: after
equalizing library sizes, each count is ``NB(mean=λ, size=r=1/φ)`` and the conditional law of
``y_a`` given the total ``t = y_a + y_b`` is **Beta-Binomial(t, r, r)** — λ cancels — so both the
qCML dispersion fit and the two-sided ("minimum-likelihood") exact p-value are beta-binomial
computations (``scipy.stats.betabinom``). This is the same reduction edgeR's ``exactTest`` uses.

The per-clonotype :func:`vdjtools.dynamics.test_pair` (Ayestaran N_eff + Fisher) is the default
caller; this one models over-dispersion explicitly instead of correcting the sample size.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.optimize import minimize_scalar
from scipy.stats import betabinom

from ..biomarker.stats import fdr_bh
from ..io.schema import COUNT
from .paired import DEFAULT_KEY, _joined


def _tmm_factor(c_a: np.ndarray, c_b: np.ndarray, *, log_ratio_trim: float = 0.3,
                sum_trim: float = 0.05) -> float:
    """TMM normalization factor of library *b* relative to *a* (Robinson & Oshlack 2010).

    Returns the factor ``f`` such that b's effective library size is ``sum(c_b) * f`` (and a's
    is ``sum(c_a) / f`` after the pair is symmetrized by the caller). Genes present in both
    libraries contribute; M (log-fold-change) and A (mean log-expression) are trimmed, then a
    variance-weighted mean of M gives ``log2 f``.
    """
    na, nb = c_a.sum(), c_b.sum()
    keep = (c_a > 0) & (c_b > 0)
    if keep.sum() < 2:
        return 1.0
    a, b = c_a[keep].astype(float), c_b[keep].astype(float)
    m = np.log2((b / nb) / (a / na))                       # log fold change (b over a)
    aval = 0.5 * (np.log2(b / nb) + np.log2(a / na))       # mean log expression
    w = (nb - b) / (nb * b) + (na - a) / (na * a)          # asymptotic variance of M
    # double-trim on M and A
    def _trim_mask(x, frac):
        lo, hi = np.quantile(x, [frac, 1 - frac])
        return (x >= lo) & (x <= hi)
    sel = _trim_mask(m, log_ratio_trim) & _trim_mask(aval, sum_trim)
    if sel.sum() == 0 or not np.isfinite(w[sel]).any() or w[sel].sum() == 0:
        return 1.0
    log2f = np.average(m[sel], weights=1.0 / w[sel])
    return float(2.0 ** log2f)


def _equalize(c_a: np.ndarray, c_b: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """TMM-normalize, then scale both libraries to their common (geometric-mean) size."""
    na, nb = c_a.sum(), c_b.sum()
    f = _tmm_factor(c_a, c_b)
    eff_a, eff_b = na / np.sqrt(f), nb * np.sqrt(f)        # symmetric split of the factor
    common = np.sqrt(eff_a * eff_b)
    q_a = c_a * (common / eff_a)
    q_b = c_b * (common / eff_b)
    return q_a, q_b, common


def _common_dispersion(q_a: np.ndarray, q_b: np.ndarray) -> float:
    """qCML common dispersion φ: maximize the beta-binomial conditional likelihood over r=1/φ."""
    t = np.rint(q_a + q_b).astype(np.int64)
    ya = np.rint(q_a).astype(np.int64)
    use = t > 0
    if use.sum() < 2:
        return 0.1
    t, ya = t[use], ya[use]

    def _negll(log_phi):
        r = 1.0 / np.exp(log_phi)
        return -np.sum(betabinom.logpmf(ya, t, r, r))

    res = minimize_scalar(_negll, bounds=(np.log(1e-4), np.log(1e3)), method="bounded")
    return float(np.exp(res.x))


def _exact_p(ya: np.ndarray, t: np.ndarray, r: float) -> np.ndarray:
    """Two-sided minimum-likelihood beta-binomial exact p-value per clonotype."""
    p = np.ones(ya.shape, dtype=float)
    for i in range(ya.size):
        ti = int(t[i])
        if ti == 0:
            continue
        pmf = betabinom.pmf(np.arange(ti + 1), ti, r, r)
        obs = pmf[int(ya[i])]
        p[i] = min(1.0, pmf[pmf <= obs * (1 + 1e-9)].sum())   # small-p two-sided
    return p


def expansion_test(a: pl.DataFrame, b: pl.DataFrame, *, key=DEFAULT_KEY,
                   dispersion: float | str = "auto", min_total: int = 6,
                   log2fc: float = 1.0, alpha: float = 0.05) -> pl.DataFrame:
    """Call clonotypes expanded/contracted between two samples by the edgeR NB exact test.

    Args:
        a: The earlier sample (canonical clonotype frame).
        b: The later sample (direction is b relative to a).
        key: Clonotype match key (default CDR3 aa + V + J).
        dispersion: ``"auto"`` fits the common φ by qCML; a float pins it.
        min_total: Clonotypes with a combined equalized count below this are left ``untested``
            (the discrete tail cannot reach a small p).
        log2fc: A clonotype is called only if ``|log2 FC| >= log2fc`` **and** ``q < alpha``.
            The paper used ``5`` for vaccine-association; ``1`` (a doubling) is a gentler default.
        alpha: BH-FDR threshold.

    Returns:
        One row per clonotype: the ``key`` columns, ``count_a``/``count_b`` (raw), ``q_a``/``q_b``
        (equalized), ``log2fc``, ``p_value``, ``q_value``, ``dispersion``, and ``call`` — one of
        ``expanded``/``contracted``/``unchanged``/``untested``. Sorted by ``q_value``.
    """
    key = list(key)
    j = _joined(a, b, key)
    c_a = j["count_a"].to_numpy().astype(float)
    c_b = j["count_b"].to_numpy().astype(float)
    q_a, q_b, _ = _equalize(c_a, c_b)

    phi = _common_dispersion(q_a, q_b) if dispersion == "auto" else float(dispersion)
    r = 1.0 / max(phi, 1e-6)

    ya = np.rint(q_a).astype(np.int64)
    yb = np.rint(q_b).astype(np.int64)
    t = ya + yb
    tested = t >= min_total

    p = np.ones(t.shape, dtype=float)
    p[tested] = _exact_p(ya[tested], t[tested], r)
    q = np.ones(t.shape, dtype=float)
    if tested.any():
        q[tested] = fdr_bh(p[tested])

    lfc = np.log2((q_b + 0.5) / (q_a + 0.5))
    sig = tested & (q < alpha) & (np.abs(lfc) >= log2fc)
    call = np.where(~tested, "untested",
                    np.where(sig & (lfc > 0), "expanded",
                             np.where(sig, "contracted", "unchanged")))

    return (j.select(key)
            .with_columns(pl.Series("count_a", c_a.astype(np.int64)),
                          pl.Series("count_b", c_b.astype(np.int64)),
                          pl.Series("q_a", q_a), pl.Series("q_b", q_b),
                          pl.Series("log2fc", lfc), pl.Series("p_value", p),
                          pl.Series("q_value", q), pl.lit(phi).alias("dispersion"),
                          pl.Series("call", call))
            .sort("q_value"))


def _demo() -> None:
    """Self-check: a planted 100x-expanded clone is called, an unchanged one is not."""
    rng = np.random.default_rng(0)
    from ..io.schema import J_CALL, JUNCTION_AA, V_CALL

    n = 500
    cdrs = ["CASS" + "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), 5)) + "F"
            for _ in range(n)]
    base = rng.integers(1, 40, n)

    def frame(counts):
        return pl.DataFrame({JUNCTION_AA: cdrs, V_CALL: ["TRBV9"] * n,
                             J_CALL: ["TRBJ2-3"] * n, COUNT: [int(x) for x in counts]})

    post = base.copy()
    post[0] = base[0] * 100          # clone 0 expands 100x
    res = expansion_test(frame(base), frame(post), log2fc=2.0)
    row0 = res.filter((pl.col(JUNCTION_AA) == cdrs[0]))
    assert row0["call"][0] == "expanded", row0
    # a randomly-picked untouched clone should not be called expanded
    other = res.filter((pl.col(JUNCTION_AA) == cdrs[250]))
    assert other["call"][0] in ("unchanged", "untested"), other
    print("expansion._demo OK: dispersion=%.3f\n" % res["dispersion"][0],
          res.head(3).select(JUNCTION_AA, "count_a", "count_b", "log2fc", "q_value", "call"))


if __name__ == "__main__":
    _demo()
