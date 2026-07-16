"""Vectorised 2×2 contingency-test kernels for biomarker association / co-occurrence.

Every kernel takes the four cells of the subject-incidence table as **numpy int arrays**
(one entry per feature / feature-pair) and returns numpy arrays — there is no per-feature
Python loop, so millions of features score in a handful of ``scipy`` calls. The cell
convention throughout is

    ===============  ===========  ===========
                     condition +  condition −
    ===============  ===========  ===========
    feature present  ``a``        ``b``
    feature absent   ``c``        ``d``
    ===============  ===========  ===========

so ``a`` = subjects that are condition-positive **and** carry the feature, ``n = a+b+c+d``,
``n_pos = a+c``, ``n_neg = b+d``, ``present = a+b``. Kernels: Fisher (hypergeometric tail),
Pearson χ² (Yates), Haldane odds ratio, a normal-approximation Bayesian posterior over the
log odds-ratio, a Beta-Binomial Bayes factor, Cochran–Mantel–Haenszel (stratified), a
label-permutation null, and Benjamini-Hochberg FDR.
"""
from __future__ import annotations

import numpy as np
from scipy.special import betaln
from scipy.stats import chi2 as _chi2
from scipy.stats import false_discovery_control, hypergeom, norm

_ALTERNATIVES = ("greater", "less", "two-sided")


def _cells(a, b, c, d):
    """Coerce the four cells to int64 arrays of a common shape."""
    a, b, c, d = (np.asarray(x, dtype=np.int64) for x in (a, b, c, d))
    return np.broadcast_arrays(a, b, c, d)


def fisher_p(a, b, c, d, alternative: str = "greater") -> np.ndarray:
    """Fisher's exact p-value per feature, as the hypergeometric tail (vectorised).

    ``"greater"`` tests enrichment in condition+ (one-tailed, the CMV setting), ``"less"``
    depletion, ``"two-sided"`` either (doubled smaller tail, capped at 1) — identical to the
    convention in :func:`vdjtools.biomarker.fisher_association`.
    """
    if alternative not in _ALTERNATIVES:
        raise ValueError(f"alternative must be one of {_ALTERNATIVES}; got {alternative!r}")
    a, b, c, d = _cells(a, b, c, d)
    n = a + b + c + d
    n_pos = a + c
    present = a + b
    p_greater = hypergeom.sf(a - 1, n, n_pos, present)   # P(X >= a)
    p_less = hypergeom.cdf(a, n, n_pos, present)         # P(X <= a)
    if alternative == "greater":
        return np.asarray(p_greater)
    if alternative == "less":
        return np.asarray(p_less)
    return np.minimum(1.0, 2.0 * np.minimum(p_greater, p_less))


def chi2_p(a, b, c, d, yates: bool = True) -> np.ndarray:
    """Pearson χ² p-value per feature for the 2×2 table (Yates-corrected by default).

    Matches ``scipy.stats.chi2_contingency([[a,b],[c,d]], correction=yates)``. Any table with
    an empty margin yields χ²=0, p=1.
    """
    a, b, c, d = _cells(a, b, c, d)
    n = (a + b + c + d).astype(np.float64)
    r1, r2, c1, c2 = a + b, c + d, a + c, b + d
    denom = (r1 * r2 * c1 * c2).astype(np.float64)
    num_ad_bc = np.abs(a * d - b * c).astype(np.float64)
    if yates:
        num_ad_bc = np.maximum(0.0, num_ad_bc - n / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.where(denom > 0, n * num_ad_bc**2 / denom, 0.0)
    return np.asarray(_chi2.sf(stat, 1))


def odds_ratio(a, b, c, d) -> np.ndarray:
    """Haldane–Anscombe (``+0.5``) odds ratio per feature — never 0/∞."""
    a, b, c, d = _cells(a, b, c, d)
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def direction(a, b, c, d) -> np.ndarray:
    """``"enriched"`` if the feature is over-represented in condition+, else ``"depleted"``."""
    a, b, c, d = _cells(a, b, c, d)
    n_pos, n_neg = a + c, b + d
    return np.where(a * n_neg > b * n_pos, "enriched", "depleted")


def bayes_logodds(a, b, c, d, alpha: float = 0.5) -> dict[str, np.ndarray]:
    """Normal-approximation Bayesian posterior over the log odds-ratio (Woolf).

    With a flat prior on ``log OR`` the posterior is ``Normal(logor, se²)`` where ``logor`` is
    the Haldane-corrected log odds-ratio and ``se² = Σ 1/(cell+alpha)`` (Woolf standard error;
    ``alpha=0.5`` = Haldane). Returns ``logor``, ``se``, a 95% credible interval
    (``ci_lo``/``ci_hi``), and ``p_or_gt1 = P(OR>1) = Φ(logor/se)``.
    """
    a, b, c, d = _cells(a, b, c, d)
    fa, fb, fc, fd = a + alpha, b + alpha, c + alpha, d + alpha
    logor = np.log(fa) + np.log(fd) - np.log(fb) - np.log(fc)
    se = np.sqrt(1.0 / fa + 1.0 / fb + 1.0 / fc + 1.0 / fd)
    return {
        "logor": logor,
        "se": se,
        "ci_lo": logor - 1.959963984540054 * se,
        "ci_hi": logor + 1.959963984540054 * se,
        "p_or_gt1": norm.cdf(logor / se),
    }


def bayes_bf(a, b, c, d, alpha: float = 1.0, beta: float = 1.0) -> np.ndarray:
    """log Bayes factor BF₁₀ for association, Beta-Binomial (analytic, vectorised).

    Compares H₁ (the feature's presence rate differs between conditions: independent
    ``Beta(alpha,beta)`` priors on the two rates) against H₀ (a single shared rate). The
    binomial coefficients cancel, leaving

        log BF₁₀ = betaln(a+α, c+β) + betaln(b+α, d+β) − betaln(α,β) − betaln(a+b+α, c+d+β)

    (positive ⇒ evidence for association). ``alpha=beta=1`` is a uniform prior.
    """
    a, b, c, d = _cells(a, b, c, d)
    return (betaln(a + alpha, c + beta) + betaln(b + alpha, d + beta)
            - betaln(alpha, beta) - betaln(a + b + alpha, c + d + beta))


def cmh(a, b, c, d) -> dict[str, np.ndarray]:
    """Cochran–Mantel–Haenszel stratified test — combined OR + χ² over strata.

    ``a,b,c,d`` are 2-D arrays of shape ``(n_features, n_strata)`` (the same 2×2 table per
    feature, one column per stratum). Returns the Mantel-Haenszel odds ratio ``or_mh``, the
    continuity-corrected ``chi2`` statistic, and its 1-df ``p_value``. Strata with fewer than
    2 subjects (``n_k < 2``) contribute nothing (their variance is undefined).
    """
    a, b, c, d = (np.asarray(x, dtype=np.float64) for x in (a, b, c, d))
    a, b, c, d = np.broadcast_arrays(a, b, c, d)
    if a.ndim == 1:
        a, b, c, d = (x[:, None] for x in (a, b, c, d))
    n = a + b + c + d
    ok = n >= 2
    nz = np.where(n > 0, n, 1.0)                       # avoid /0 in masked-out strata
    r1, c1 = a + b, a + c
    num_or = np.where(ok, a * d / nz, 0.0).sum(axis=1)
    den_or = np.where(ok, b * c / nz, 0.0).sum(axis=1)
    or_mh = np.where(den_or > 0, num_or / np.where(den_or > 0, den_or, 1.0), np.inf)
    exp = np.where(ok, r1 * c1 / nz, 0.0).sum(axis=1)
    var = np.where(ok & (n > 1), r1 * (c + d) * c1 * (b + d) / (nz**2 * (nz - 1.0)), 0.0).sum(axis=1)
    obs = np.where(ok, a, 0.0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.where(var > 0, (np.abs(obs - exp) - 0.5).clip(min=0.0) ** 2 / var, 0.0)
    return {"or_mh": or_mh, "chi2": stat, "p_value": np.asarray(_chi2.sf(stat, 1))}


def permutation_p(present: np.ndarray, labels: np.ndarray, *, n_perm: int = 1000,
                  seed: int = 0, alternative: str = "greater") -> np.ndarray:
    """Label-permutation p-value per feature from the subject-incidence matrix.

    ``present`` is a boolean ``(n_subjects, n_features)`` incidence matrix, ``labels`` a boolean
    ``(n_subjects,)`` condition vector. Permuting ``labels`` holds each feature's incidence
    fixed, so the null is exactly the hypergeometric one — this converges to
    :func:`fisher_p` as ``n_perm → ∞`` and is a robust check for sparse tables. Uses the
    add-one estimator ``(1 + #exceed)/(1 + n_perm)``. Seeded → reproducible.
    """
    if alternative not in ("greater", "less"):
        raise ValueError("permutation supports 'greater' or 'less'")
    present = np.asarray(present, dtype=np.float64)
    labels = np.asarray(labels, dtype=bool)
    a_obs = present[labels].sum(axis=0)                        # per-feature condition+ count
    rng = np.random.default_rng(seed)
    n_sub = present.shape[0]
    n_pos = int(labels.sum())
    exceed = np.zeros(present.shape[1], dtype=np.int64)
    for _ in range(n_perm):
        idx = rng.permutation(n_sub)[:n_pos]
        a_perm = present[idx].sum(axis=0)
        if alternative == "greater":
            exceed += a_perm >= a_obs
        else:
            exceed += a_perm <= a_obs
    return (1.0 + exceed) / (1.0 + n_perm)


def fdr_bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (empty-safe)."""
    p = np.asarray(p, dtype=np.float64)
    return false_discovery_control(p, method="bh") if p.size else p
