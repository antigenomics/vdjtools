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

_ALTERNATIVES = ("greater", "less", "two-sided", "two-sided-minlike")

#: R's ``fisher.test`` relative-error slack when comparing a table's likelihood to the
#: observed one; scipy's ``fisher_exact`` uses the same guard. Without it, floating-point
#: noise drops tables that are exactly as likely as the observed table.
_MINLIKE_RELERR = 1 + 1e-7


def _cells(a, b, c, d):
    """Coerce the four cells to int64 arrays of a common shape."""
    a, b, c, d = (np.asarray(x, dtype=np.int64) for x in (a, b, c, d))
    return np.broadcast_arrays(a, b, c, d)


def _two_sided_minlike(a, n, n_pos, present) -> np.ndarray:
    """Minimum-likelihood two-sided Fisher p: Σ P(X=k) over all k at most as likely as ``a``.

    Vectorised by grouping on the hypergeometric parameters ``(n, n_pos, present)`` — every
    feature in a group shares one support and one pmf, so each group costs a single
    ``hypergeom.pmf`` over its support plus a ``searchsorted``. In the paired-dynamics setting
    the two library sizes are fixed per pair, so a group is just "features with the same
    combined count", and the support is ``0..(a+b)`` because the library sizes dwarf it.
    """
    out = np.zeros(a.shape, dtype=np.float64)
    flat_a, flat_out = a.reshape(-1), out.reshape(-1)
    keys = np.stack([n.reshape(-1), n_pos.reshape(-1), present.reshape(-1)], axis=-1)
    uniq, inv = np.unique(keys, axis=0, return_inverse=True)
    for g, (gn, gk, gd) in enumerate(uniq):
        sel = inv == g
        lo, hi = max(0, gd - (gn - gk)), min(gk, gd)
        support = np.arange(lo, hi + 1)
        pmf = hypergeom.pmf(support, gn, gk, gd)
        order = np.argsort(pmf, kind="stable")
        cum = np.cumsum(pmf[order])
        p_obs = hypergeom.pmf(flat_a[sel], gn, gk, gd)
        # Sum the likelihood of every table no more likely than the observed one.
        idx = np.searchsorted(pmf[order], p_obs * _MINLIKE_RELERR, side="right")
        flat_out[sel] = np.where(idx > 0, cum[np.maximum(idx - 1, 0)], 0.0)
    return np.minimum(out, 1.0)


def fisher_p(a, b, c, d, alternative: str = "greater") -> np.ndarray:
    """Fisher's exact p-value per feature, as the hypergeometric tail (vectorised).

    ``"greater"`` tests enrichment in condition+ (one-tailed, the CMV setting), ``"less"``
    depletion, ``"two-sided"`` either (doubled smaller tail, capped at 1) — identical to the
    convention in :func:`vdjtools.biomarker.fisher_association`.

    ``"two-sided-minlike"`` is the *minimum-likelihood* two-sided convention used by R's
    ``fisher.test`` and :func:`scipy.stats.fisher_exact`: the sum of the likelihood of every
    table at most as likely as the observed one. It differs from ``"two-sided"`` whenever the
    margins are not exactly symmetric — measured on paired-dynamics-shaped tables, exactly
    equal margins agree 300/300, but margins differing by rounding alone disagree 171/300 with
    ``"two-sided"`` up to exactly 2× larger. Use it when comparing against R or scipy, or when
    the margins are only near-equal; ``"two-sided"`` is kept as the default two-sided
    convention because the association/co-occurrence results were computed with it.
    """
    if alternative not in _ALTERNATIVES:
        raise ValueError(f"alternative must be one of {_ALTERNATIVES}; got {alternative!r}")
    a, b, c, d = _cells(a, b, c, d)
    n = a + b + c + d
    n_pos = a + c
    present = a + b
    if alternative == "two-sided-minlike":
        return _two_sided_minlike(a, n, n_pos, present)
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
    # den_or==0 with num_or>0 is a genuinely infinite OR; den_or==0 AND num_or==0 is 0/0 —
    # UNDEFINED, so nan. (Reachable for a near-ubiquitous pair, where b=c=d=0 in every
    # stratum; inf there would rank degenerate pairs top under a `or_mh > x` filter.)
    _ratio = num_or / np.where(den_or > 0, den_or, 1.0)
    or_mh = np.where(den_or > 0, _ratio, np.where(num_or > 0, np.inf, np.nan))
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

    A permutation's counts are ``L @ present`` with ``L`` the 0/1 permuted-label matrix, so a
    chunk of permutations is one BLAS call instead of ``n_perm`` fancy-index copies (~15x). The
    chunk still draws ``rng.permutation(n_sub)`` once per permutation in the original order, so
    p-values are unchanged for a given seed; and summing 0/1 to an integer <= n_sub is exact in
    float64 whatever order BLAS accumulates in.
    """
    if alternative not in ("greater", "less"):
        raise ValueError("permutation supports 'greater' or 'less'")
    present = np.asarray(present, dtype=np.float64)
    labels = np.asarray(labels, dtype=bool)
    a_obs = present[labels].sum(axis=0)                        # per-feature condition+ count
    rng = np.random.default_rng(seed)
    n_sub, n_feat = present.shape
    n_pos = int(labels.sum())
    exceed = np.zeros(n_feat, dtype=np.int64)
    # Bound the (chunk x n_feat) intermediate to ~64 MB: at 100k features a single full
    # (n_perm x n_feat) matmul would allocate ~800 MB.
    chunk = max(1, min(n_perm, 8_000_000 // max(n_feat, 1)))
    for start in range(0, n_perm, chunk):
        k = min(chunk, n_perm - start)
        L = np.zeros((k, n_sub))
        for r in range(k):
            L[r, rng.permutation(n_sub)[:n_pos]] = 1.0
        a_perm = L @ present                                   # [k, n_feat]
        exceed += ((a_perm >= a_obs) if alternative == "greater"
                   else (a_perm <= a_obs)).sum(axis=0)
    return (1.0 + exceed) / (1.0 + n_perm)


def fdr_bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR-adjusted q-values (empty-safe)."""
    p = np.asarray(p, dtype=np.float64)
    return false_discovery_control(p, method="bh") if p.size else p
