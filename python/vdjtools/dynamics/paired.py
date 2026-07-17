"""Paired within-donor clonotype test: per-pair ``N_eff``, downscale, two-tailed Fisher.

The method is Ayestaran (PhD, Cambridge 2024), Ch. 2. Sequencing is a **two-step** sampling
process — a true frequency ``f`` is sampled into ``N_S1`` molecules, which are sampled into
``N_seq`` reads — so the sample size driving the noise is the harmonic sum (Eq. 2.5–2.6)::

    1 / N_eff = 1 / N_S1 + 1 / N_seq

and is dominated by the *smaller* step. Where ``N_S1 << N_seq`` (the usual case: ~1e6 PBMCs
into ~2e7 reads) the library is heavily oversampled, the observed counts are far noisier than
``N_seq`` implies, and a Fisher test that assumes ``N_seq`` is the sample size calls enormous
numbers of clones changed when nothing changed (thesis Fig. 2.13).

``N_eff`` is a property of the **pair**, not of a sample or a cohort: it is read off the
mean–variance relationship *between two samples*, so a global ``N_eff`` is not merely a bad
idea, it is undefined. That is why :func:`_downscale` is private — a public "rescale to N"
would invite exactly the cohort-wide normalisation this method exists to avoid.

Estimation fits the thesis's Eq. 2.3, ``log σ² = log f + log(1/N)``, with the **slope fixed at
1** — which is just a weighted mean of ``log σ² − log f``, so it needs numpy, not a regression
library. It differs from the thesis in *which* variance it fits, and that choice is
load-bearing enough to record here.

The thesis (p. 19) bins clones by their frequency in sample **B** and takes the mean and
variance of the **A** frequencies in each bin, then **multiplies by 2**: binning on the noisy
``f_B`` instead of the true ``f`` mixes a range of true frequencies into each bin, roughly
doubling the apparent variance, so the fitted value is half the real ``N_eff``.

We instead bin on the mean ``f̄ = (f_A + f_B)/2`` and fit ``var(f_A − f_B) = 2f/N``. The
difference **cancels the true frequency**, so there is no mis-binning to correct and no
correction factor to mis-set — arguably closer to Eq. 2.3's own derivation, which assumed
``f`` known.

Why, measured on simulated two-step data against a *planted* ``N_eff`` (see
``tests/python/test_dynamics_paired.py``): the ``×2`` is exact only when the clone-size
distribution is flat, because only then is the bin contamination equal to the sampling
variance. Real repertoires are heavy-tailed, where the contamination is smaller and a fixed
``×2`` overshoots by ~25%. That is not cosmetic — it inflates the false-positive rate ~1.27×
and **fails the thesis's own p-value-uniformity acceptance test (Fig. 2.13) in 3 of 4 regimes**,
including precisely the oversampled ones (``N_S1 << N_seq``) the method exists to handle::

    regime (N_S1/N_seq)   planted    thesis x2   gate      difference   gate
    200k / 2M             181,818    1.25x       FAIL      1.07x        pass
     50k / 2M              48,780    1.22x       FAIL      1.02x        pass
    200k / 1M             166,667    1.23x       FAIL      0.99x        pass
      1M / 300k           230,769    0.97x       pass      0.95x        pass

The test itself is *not* at fault: pinned to the planted ``N_eff`` it is perfectly calibrated
(conservative, as an exact test on discrete counts should be). The bias was entirely the
estimator's.

The test itself is Eq. 2.4: conditioning on the total ``r1 + r2`` cancels the unknown true
frequency exactly, leaving a hypergeometric — i.e. a two-tailed Fisher exact test on
``[[c_a, c_b], [R_a - c_a, R_b - c_b]]``. It uses the **minimum-likelihood** two-sided
convention (R / scipy), because rounding leaves the two library sizes near-equal rather than
equal, and on near-equal margins the doubling convention differs by up to 2x.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import poisson

from ..biomarker.stats import fdr_bh, fisher_p
from ..io.schema import COUNT, J_CALL, JUNCTION_AA, V_CALL

#: Default clonotype match key (CDR3 aa + V + J), as in :mod:`vdjtools.overlap.track`.
DEFAULT_KEY = (JUNCTION_AA, V_CALL, J_CALL)

#: The five dynamics classes, plus ``untested`` for clones below the testability floor.
CLASSES = ("emergent", "expanded", "persistent", "contracted", "vanishing", "untested")


def _agg(df: pl.DataFrame, key: list[str]) -> pl.DataFrame:
    """Collapse to unique ``key`` with summed counts and within-sample frequency."""
    g = df.group_by(key).agg(pl.col(COUNT).sum().alias("_c"))
    total = g["_c"].sum() or 1
    return g.with_columns((pl.col("_c") / total).alias("_f"))


def _joined(a: pl.DataFrame, b: pl.DataFrame, key: list[str]) -> pl.DataFrame:
    """Full outer join of two samples on ``key``; absent clonotypes get count 0, freq 0.

    A clonotype missing from one sample is a real zero for this test — it was not seen — as
    opposed to "not sampled", which is what a narrower key would silently create.
    """
    ja = _agg(a, key).rename({"_c": "count_a", "_f": "f_a"})
    jb = _agg(b, key).rename({"_c": "count_b", "_f": "f_b"})
    return (ja.join(jb, on=key, how="full", coalesce=True)
              .with_columns(pl.col("count_a", "count_b").fill_null(0),
                            pl.col("f_a", "f_b").fill_null(0.0)))


def _bin_index(x: np.ndarray, bins: int) -> np.ndarray:
    """Log-spaced bin index for strictly positive ``x``."""
    edges = np.geomspace(x.min(), x.max(), bins + 1)
    return np.clip(np.digitize(x, edges) - 1, 0, bins - 1)


def _neff_fit(f_a: np.ndarray, f_b: np.ndarray, bins: int, min_bin: int) -> float:
    """Weighted slope-1 fit of ``log var(f_a - f_b) = log 2f̄ + log(1/N)``; returns ``N``.

    Two independent observations of the same true ``f`` each carry variance ``f/N``, so their
    difference carries ``2f/N`` — and the true ``f`` cancels, which is what removes the thesis's
    ``×2`` correction rather than merely re-tuning it. The weight is the bin's clone count: the
    sampling variance of a log-variance estimate scales as ``1/n``.
    """
    fbar = 0.5 * (f_a + f_b)
    idx = _bin_index(fbar, bins)
    betas, weights = [], []
    for g in range(bins):
        sel = idx == g
        n = int(sel.sum())
        if n < min_bin:
            continue
        var, mean = (f_a[sel] - f_b[sel]).var(ddof=1), fbar[sel].mean()
        if var <= 0 or mean <= 0:
            continue
        betas.append(np.log(var) - np.log(2.0 * mean))   # = log(1/N) for this bin
        weights.append(n)
    if len(betas) < 2:
        raise ValueError(
            f"N_eff fit needs >=2 usable frequency bins, got {len(betas)} "
            f"(from {f_a.size} clonotypes) — the pair is too shallow or shares too few clones"
        )
    return float(np.exp(-np.average(betas, weights=weights)))


def _outlier_mask(f_a: np.ndarray, f_b: np.ndarray, ref_n: float,
                  bins: int) -> np.ndarray:
    """Keep clones whose two frequencies are consistent with pure noise at ``ref_n``.

    Thesis p. 29 / Fig. A.2: real expansions inflate the binned variance and bias ``N_eff``
    **down**, so they are removed before the fit. Bin by the mean frequency across both
    samples; within a bin of ``n`` clones keep only the band beyond which fewer than one clone
    is expected under Poisson noise around that mean. This is a no-op on true replicates (no
    real dynamics to remove) and matters only on a contrast like pre-vs-post.
    """
    fbar = 0.5 * (f_a + f_b)
    idx = _bin_index(fbar, bins)
    keep = np.ones(fbar.shape, dtype=bool)
    for g in range(bins):
        sel = idx == g
        n = int(sel.sum())
        if n < 2:
            continue
        lam = ref_n * fbar[sel].mean()
        lo, hi = poisson.ppf([1 / (2 * n), 1 - 1 / (2 * n)], lam) / ref_n
        keep[sel] = (f_a[sel] >= lo) & (f_a[sel] <= hi) & (f_b[sel] >= lo) & (f_b[sel] <= hi)
    return keep


def estimate_neff(a: pl.DataFrame, b: pl.DataFrame, *, key=DEFAULT_KEY,
                  min_count: int = 2, bins: int = 25, min_bin: int = 10,
                  ref_n: float | None = None) -> float:
    """Estimate the pair's effective sample size from its own mean–variance scaling.

    Args:
        a: First sample of the pair (canonical clonotype frame).
        b: Second sample of the pair.
        key: Clonotype match key (default CDR3 aa + V + J).
        min_count: Discreteness floor for the **fit** — clonotypes whose mean frequency across
            the pair is below this many counts are excluded from the mean–variance fit, because
            at 1–2 counts the frequency is too discrete to estimate a variance from. This is
            *not* the testability floor used by :func:`test_pair` (see its ``min_total``); it is
            the most sensitive knob here, so the benchmark sweeps it rather than trusting a
            default.
        bins: Number of log-spaced frequency bins.
        min_bin: Minimum clonotypes in a bin for it to contribute to the fit.
        ref_n: Reference sample size for the outlier pre-filter. ``None`` (default) derives it
            from the data with an unfiltered first pass — the thesis hard-codes 200,000, which
            is a property of *its* cohort (2x one of our datasets' whole libraries and 33x
            another's), not a constant. Pass a number to pin it.

    Returns:
        The estimated ``N_eff``.

    Raises:
        ValueError: If fewer than two frequency bins are usable — a pair too shallow, or with
            too few shared clonotypes, to fit. Never returns a silent fallback: a wrong
            ``N_eff`` silently mis-scales every downstream p-value.
    """
    j = _joined(a, b, list(key))
    f_a, f_b = j["f_a"].to_numpy(), j["f_b"].to_numpy()
    n_seq = min(int(j["count_a"].sum()), int(j["count_b"].sum()))
    sel = 0.5 * (f_a + f_b) > (min_count / max(n_seq, 1))
    f_a, f_b = f_a[sel], f_b[sel]
    if f_a.size == 0:
        raise ValueError("no clonotypes clear min_count; the pair is too shallow")

    if ref_n is None:                      # first pass, unfiltered, only to scale the filter
        ref_n = _neff_fit(f_a, f_b, bins, min_bin)
    keep = _outlier_mask(f_a, f_b, ref_n, bins)
    if keep.sum() >= min_bin * 2:          # else the filter left too little to fit — use all
        f_a, f_b = f_a[keep], f_b[keep]
    return _neff_fit(f_a, f_b, bins, min_bin)


def _downscale(f: np.ndarray, neff: float) -> np.ndarray:
    """Thesis Eq. 2.7: ``c_i = round(f_i * N_eff)``.

    Deliberately deterministic and deliberately private. This is *not*
    :func:`vdjtools.preprocess.downsample`, which draws a fresh multivariate-hypergeometric
    sample: drawing again would add a second layer of sampling noise on top of the one
    ``N_eff`` already encodes, and the Fisher null below assumes the noise is exactly the one
    it models.
    """
    return np.rint(f * neff).astype(np.int64)


def test_pair(a: pl.DataFrame, b: pl.DataFrame, *, key=DEFAULT_KEY,
              neff: float | str = "auto", min_total: int = 6, alpha: float = 0.01,
              **neff_kw) -> pl.DataFrame:
    """Test every clonotype for a within-donor frequency change between two samples.

    Downscales both samples to the pair's ``N_eff`` and applies a two-tailed Fisher exact test
    to ``[[c_a, c_b], [R_a - c_a, R_b - c_b]]``. Conditioning on ``c_a + c_b`` cancels the
    unknown true frequency exactly (thesis Eq. 2.4), so the test needs no estimate of it.

    Args:
        a: The earlier sample (e.g. pre-vaccination).
        b: The later sample (e.g. post-vaccination). Direction is reported b-relative-to-a.
        key: Clonotype match key.
        neff: ``"auto"`` estimates it from the pair via :func:`estimate_neff`; a float pins it;
            ``None`` skips the downscale entirely — correct only when the counts are already
            molecule/UMI counts rather than reads (thesis p. 86), since then there is no
            oversampling to undo.
        min_total: Testability floor — clonotypes with a combined downscaled count below this
            are classed ``untested`` rather than tested and called non-significant. Below ~6
            the discrete hypergeometric cannot reach a small p at all, so testing them only
            costs multiple-testing burden. Distinct from ``estimate_neff``'s ``min_count``.
        alpha: FDR threshold (Benjamini–Hochberg q) for calling a change significant.
        **neff_kw: Forwarded to :func:`estimate_neff` when ``neff="auto"``.

    Returns:
        One row per clonotype with the ``key`` columns, ``count_a``/``count_b`` (downscaled),
        ``f_a``/``f_b`` (original within-sample frequencies), ``p_value``, ``q_value`` and
        ``dynamics`` — one of ``emergent`` (absent from ``a``), ``expanded``, ``persistent``
        (no evidence of change), ``contracted``, ``vanishing`` (absent from ``b``), or
        ``untested``. The classes partition the frame.

    Raises:
        ValueError: If ``neff="auto"`` and the pair cannot be fit (see :func:`estimate_neff`).
    """
    key = list(key)
    j = _joined(a, b, key)
    f_a, f_b = j["f_a"].to_numpy(), j["f_b"].to_numpy()

    if neff == "auto":
        neff = estimate_neff(a, b, key=key, **neff_kw)
    if neff is None:                       # molecule counts: the observed counts ARE the draw
        c_a, c_b = j["count_a"].to_numpy(), j["count_b"].to_numpy()
    else:
        c_a, c_b = _downscale(f_a, neff), _downscale(f_b, neff)

    # R_a / R_b are the REALIZED downscaled library sizes, not N_eff: rounding makes them
    # near-equal but not equal, which is exactly why the min-likelihood convention is required.
    r_a, r_b = int(c_a.sum()), int(c_b.sum())
    total = c_a + c_b
    tested = total >= min_total

    p = np.ones(total.shape, dtype=float)
    p[tested] = fisher_p(c_a[tested], c_b[tested], r_a - c_a[tested], r_b - c_b[tested],
                         alternative="two-sided-minlike")
    q = np.ones(total.shape, dtype=float)
    if tested.any():
        q[tested] = fdr_bh(p[tested])

    sig = tested & (q < alpha)
    up = (c_b / max(r_b, 1)) > (c_a / max(r_a, 1))
    dynamics = np.where(
        ~tested, "untested",
        np.where(sig & (c_a == 0), "emergent",
                 np.where(sig & (c_b == 0), "vanishing",
                          np.where(sig & up, "expanded",
                                   np.where(sig, "contracted", "persistent")))))

    return (j.select(key)
             .with_columns(pl.Series("count_a", c_a), pl.Series("count_b", c_b),
                           pl.Series("f_a", f_a), pl.Series("f_b", f_b),
                           pl.Series("p_value", p), pl.Series("q_value", q),
                           pl.Series("dynamics", dynamics))
             .sort("q_value"))
