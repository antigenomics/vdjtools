"""Oracle tests for the biomarker 2×2 statistical kernels (vs scipy / brute force)."""
import numpy as np
import pytest
from scipy import integrate
from scipy.stats import betabinom, chi2_contingency, fisher_exact

from vdjtools.biomarker import stats

# A spread of 2×2 tables: enriched, depleted, null, sparse, and a zero-margin edge.
TABLES = [
    (8, 1, 2, 9),      # enriched
    (1, 8, 9, 2),      # depleted
    (5, 5, 5, 5),      # null
    (3, 0, 0, 7),      # perfect separation
    (0, 0, 4, 6),      # empty "present" margin
    (12, 7, 40, 60),   # larger, mild
]
A, B, C, D = (np.array(x) for x in zip(*TABLES))


@pytest.mark.parametrize("alt", ["greater", "less"])
def test_fisher_one_tailed_matches_scipy(alt):
    got = stats.fisher_p(A, B, C, D, alternative=alt)
    for i, (a, b, c, d) in enumerate(TABLES):
        _, exp = fisher_exact([[a, b], [c, d]], alternative=alt)
        assert got[i] == pytest.approx(exp, rel=1e-9, abs=1e-12)


def test_fisher_two_sided_is_doubled_min_tail():
    # The module's documented two-sided convention (matching fisher_association) is the
    # doubled smaller one-tail, NOT scipy's "sum of tables <= observed" exact two-sided.
    got = stats.fisher_p(A, B, C, D, alternative="two-sided")
    for i, (a, b, c, d) in enumerate(TABLES):
        _, pg = fisher_exact([[a, b], [c, d]], alternative="greater")
        _, pl = fisher_exact([[a, b], [c, d]], alternative="less")
        assert got[i] == pytest.approx(min(1.0, 2.0 * min(pg, pl)), rel=1e-9, abs=1e-12)


@pytest.mark.parametrize("yates", [True, False])
def test_chi2_matches_scipy(yates):
    got = stats.chi2_p(A, B, C, D, yates=yates)
    for i, (a, b, c, d) in enumerate(TABLES):
        tbl = np.array([[a, b], [c, d]])
        if (tbl.sum(0) == 0).any() or (tbl.sum(1) == 0).any():
            assert got[i] == pytest.approx(1.0)          # degenerate margin -> p=1
            continue
        _, exp, _, _ = chi2_contingency(tbl, correction=yates)
        assert got[i] == pytest.approx(exp, rel=1e-9)


def test_odds_ratio_haldane():
    got = stats.odds_ratio(A, B, C, D)
    exp = [((a + .5) * (d + .5)) / ((b + .5) * (c + .5)) for a, b, c, d in TABLES]
    assert np.allclose(got, exp)
    assert np.all(np.isfinite(got))                       # never 0/inf even with a zero cell


def test_direction():
    d = stats.direction(A, B, C, D)
    assert d[0] == "enriched" and d[1] == "depleted"


def test_bayes_logodds_consistency():
    r = stats.bayes_logodds(A, B, C, D)
    # posterior mean == Haldane log-OR; P(OR>1) monotone with logor; CI brackets the mean
    assert np.allclose(np.exp(r["logor"]), stats.odds_ratio(A, B, C, D))
    assert r["p_or_gt1"][0] > 0.9 and r["p_or_gt1"][1] < 0.1   # enriched vs depleted
    assert np.all(r["ci_lo"] < r["logor"]) and np.all(r["logor"] < r["ci_hi"])
    assert r["p_or_gt1"][2] == pytest.approx(0.5, abs=1e-9)    # symmetric null table


def _bf_brute(a, b, c, d, alpha=1.0, beta=1.0):
    """Independent Bayes factor via BetaBinomial marginals + a numerical H0 integral."""
    n_pos, n_neg = a + c, b + d
    h1 = betabinom.pmf(a, n_pos, alpha, beta) * betabinom.pmf(b, n_neg, alpha, beta)
    from math import comb
    from scipy.stats import beta as beta_dist
    coef = comb(n_pos, a) * comb(n_neg, b)

    def integrand(p):
        return coef * p**(a + b) * (1 - p)**(c + d) * beta_dist.pdf(p, alpha, beta)
    h0, _ = integrate.quad(integrand, 0, 1)
    return np.log(h1 / h0)


def test_bayes_bf_matches_brute_force():
    got = stats.bayes_bf(A, B, C, D)
    for i, (a, b, c, d) in enumerate(TABLES):
        assert got[i] == pytest.approx(_bf_brute(a, b, c, d), rel=1e-7, abs=1e-9)
    # association tables give positive evidence; the null table gives negative (favours H0)
    assert got[0] > 1.0 and got[3] > 1.0
    assert got[2] < 0.0


def test_cmh_single_stratum_reduces_to_mh():
    a, b, c, d = 8, 2, 3, 12
    r = stats.cmh([a], [b], [c], [d])
    assert r["or_mh"][0] == pytest.approx((a * d) / (b * c))   # plain OR, one stratum
    n = a + b + c + d
    e = (a + b) * (a + c) / n
    v = (a + b) * (c + d) * (a + c) * (b + d) / (n**2 * (n - 1))
    chi = (abs(a - e) - 0.5) ** 2 / v
    assert r["chi2"][0] == pytest.approx(chi)


def test_cmh_confounding_two_strata():
    # Same marginal OR in each stratum but strata differ in exposure prevalence -> MH combines
    # the within-stratum effect (OR≈1 here => no real association after stratifying).
    a = np.array([[10, 10]])
    b = np.array([[10, 10]])
    c = np.array([[10, 10]])
    d = np.array([[10, 10]])
    r = stats.cmh(a, b, c, d)
    assert r["or_mh"][0] == pytest.approx(1.0)
    assert r["p_value"][0] > 0.5
    # a genuine within-stratum effect is detected
    a2 = np.array([[18, 16]])
    b2 = np.array([[4, 6]])
    c2 = np.array([[6, 4]])
    d2 = np.array([[16, 18]])
    r2 = stats.cmh(a2, b2, c2, d2)
    assert r2["or_mh"][0] > 3.0 and r2["p_value"][0] < 0.01


def test_permutation_reproducible_and_approaches_fisher():
    rng = np.random.default_rng(0)
    n_sub, n_feat = 60, 5
    labels = np.array([True] * 30 + [False] * 30)
    present = rng.random((n_sub, n_feat)) < 0.4
    # plant a strongly enriched feature 0
    present[:30, 0] = rng.random(30) < 0.8
    p1 = stats.permutation_p(present, labels, n_perm=2000, seed=1)
    p2 = stats.permutation_p(present, labels, n_perm=2000, seed=1)
    assert np.array_equal(p1, p2)                              # seeded -> reproducible
    # permutation p ≈ Fisher greater on the same incidence table
    a = present[labels].sum(0).astype(int)
    pres = present.sum(0).astype(int)
    n_pos = int(labels.sum())
    fis = stats.fisher_p(a, pres - a, n_pos - a, n_sub - n_pos - (pres - a), alternative="greater")
    assert np.allclose(p1, fis, atol=0.03)


def test_fdr_bh_empty_safe():
    assert stats.fdr_bh(np.array([])).size == 0
    q = stats.fdr_bh(np.array([0.001, 0.5, 0.9]))
    assert np.all((q >= 0) & (q <= 1)) and q[0] < q[2]
