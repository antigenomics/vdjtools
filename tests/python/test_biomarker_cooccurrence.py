"""Tests for vdjtools.biomarker.cooccurrence — α-β pairing and same-chain co-specificity."""
import numpy as np
import polars as pl
import pytest

from vdjtools.biomarker import cooccurrence


def _c(sid, v, j, aa):
    return dict(sample_id=sid, v_call=v, j_call=j, junction_aa=aa,
                junction_nt="ACG", duplicate_count=10, frequency=0.0)


def _paired_cohort(seed=1, n=50, n_pair=24):
    """n subjects with both chains; a planted α-β pair co-occurs in n_pair of them."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        sid = f"S{i:02d}"
        if i < n_pair:                                     # the true pair
            rows.append(_c(sid, "TRAV1", "TRAJ1", "CAVRF"))
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSPF"))
        rows.append(_c(sid, "TRAV2", "TRAJ2", "CAVBG"))    # background, both chains everywhere
        rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSBG"))
        if rng.random() < 0.5:                             # independent extras
            rows.append(_c(sid, "TRAV3", "TRAJ3", "CAVXX"))
        if rng.random() < 0.5:
            rows.append(_c(sid, "TRBV3", "TRBJ3", "CASSYY"))
    return pl.DataFrame(rows)


def test_alpha_beta_recovers_planted_pair():
    cohort = _paired_cohort()
    r = cooccurrence(cohort, chain_a="TRA", chain_b="TRB", min_incidence=3,
                     min_cooccurrence=3, evalue=True)
    top = r.sort("p_value").row(0, named=True)
    assert (top["a_junction_aa"], top["b_junction_aa"]) == ("CAVRF", "CASSPF")
    assert top["n_ab"] == 24 and top["theta"] > 1.8 and top["q_value"] < 1e-6
    # an independent pair is not significant
    ind = r.filter((pl.col("a_junction_aa") == "CAVXX") & (pl.col("b_junction_aa") == "CASSYY"))
    assert ind.height == 0 or ind["q_value"][0] > 0.1
    assert {"expected", "e_value"} <= set(r.columns)               # evalue columns present


def test_universe_is_subjects_with_both_chains():
    cohort = _paired_cohort(n=50)
    # add 10 alpha-only subjects — they must not enter the paired universe (n stays 50)
    extra = pl.DataFrame([_c(f"A{i}", "TRAV1", "TRAJ1", "CAVRF") for i in range(10)])
    r = cooccurrence(pl.concat([cohort, extra]), chain_a="TRA", chain_b="TRB",
                     min_incidence=3, min_cooccurrence=3)
    assert r["n"][0] == 50


def test_same_chain_upper_triangle_no_self_pairs():
    # plant two TRB clonotypes that co-occur (co-specific), plus independent background
    rng = np.random.default_rng(11)
    rows = []
    for i in range(40):
        sid = f"S{i:02d}"
        if i < 20:
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSA"))
            rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSB"))       # co-occurs with CASSA
        rows.append(_c(sid, "TRBV9", "TRBJ9", "CASSBG"))
        # Per-subject-unique filler sets repertoire depth INDEPENDENTLY of the planted pair, as
        # in a real repertoire. Without it depth ≡ pair-presence, and conditioning on depth (the
        # default) would legitimately condition the signal away.
        for f in range(rng.integers(20, 60)):
            rows.append(_c(sid, "TRBV8", "TRBJ8", f"CASF{i}_{f}"))
    cohort = pl.DataFrame(rows)
    r = cooccurrence(cohort, chain_a="TRB", chain_b=None, min_incidence=3, min_cooccurrence=3)
    # no self pairs
    assert (r["a_junction_aa"] != r["b_junction_aa"]).all()
    pair = r.filter(pl.col("a_junction_aa").is_in(["CASSA", "CASSB"])
                    & pl.col("b_junction_aa").is_in(["CASSA", "CASSB"]))
    assert pair.height == 1 and pair["n_ab"][0] == 20 and pair["p_value"][0] < 1e-6


def test_empty_result_stable_schema():
    # nothing clears min_cooccurrence -> empty frame, columns still present
    cohort = pl.DataFrame([_c("S0", "TRAV1", "TRAJ1", "CAVA"), _c("S0", "TRBV1", "TRBJ1", "CASSB"),
                           _c("S1", "TRAV2", "TRAJ2", "CAVC"), _c("S1", "TRBV2", "TRBJ2", "CASSD")])
    r = cooccurrence(cohort, chain_a="TRA", chain_b="TRB", min_incidence=2, min_cooccurrence=2)
    assert r.height == 0 and "theta" in r.columns and "a_junction_aa" in r.columns


def test_bad_test_arg():
    with pytest.raises(ValueError, match="test must be"):
        cooccurrence(_paired_cohort(), test="nope")


def _depth_confounded_cohort(seed=7, n=400):
    """Depth as the ONLY link: A and B are drawn INDEPENDENTLY given repertoire size.

    Deep subjects carry more of everything, so A and B co-occur across subjects with no pairing
    whatsoever — the pooled test is fooled; conditioning on depth must not be. Depth is set by
    filler clonotypes so that A/B are a negligible share of it (as in real repertoires), and the
    detection probability is kept away from 0/1 where no lift can be induced.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        sid = f"S{i:03d}"
        depth = max(int(np.exp(rng.normal(np.log(40), 0.9))), 2)  # skewed; CV ~ 1
        p = 1 - np.exp(-depth / 90.0)                             # ~0.36 at median, rises with depth
        # Background present in every subject, so the universe is the whole cohort and A/B
        # absence is observed (a single-feature `candidates` list would collapse the universe
        # to the subjects that HAVE the feature, making n_a == n_ab == n and theta == 1).
        rows.append(_c(sid, "TRAV2", "TRAJ2", "CAVBG"))
        rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSBG"))
        if rng.random() < p:
            rows.append(_c(sid, "TRAV1", "TRAJ1", "CAVRF"))
        if rng.random() < p:
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSPF"))
        for f in range(depth):                    # per-subject-unique filler → sets depth only
            rows.append(_c(sid, "TRAV9", "TRAJ9", f"CAVF{i}_{f}"))
            rows.append(_c(sid, "TRBV9", "TRBJ9", f"CASF{i}_{f}"))
    return pl.DataFrame(rows)


def test_depth_strata_zero_reproduces_the_pooled_test_exactly():
    """depth_strata=0 is the uncorrected pooled test — the oracle the CMH branch is checked against."""
    df = _depth_confounded_cohort()
    pooled = cooccurrence(df, min_cooccurrence=2, depth_strata=0)
    assert "or_mh" not in pooled.columns and "chi2" not in pooled.columns
    strat = cooccurrence(df, min_cooccurrence=2, depth_strata=10)
    shared = ["a_junction_aa", "b_junction_aa", "n", "n_a", "n_b", "n_ab", "theta"]
    # The pooled lift/counts are untouched by stratification; only the p-value changes.
    assert pooled.sort(shared[:2]).select(shared).equals(strat.sort(shared[:2]).select(shared))
    assert "or_mh" in strat.columns


def _ab(df):
    return df.filter((pl.col("a_junction_aa") == "CAVRF") & (pl.col("b_junction_aa") == "CASSPF"))


def test_cmh_over_depth_strata_kills_a_depth_induced_false_positive():
    """The headline: depth alone must not produce a significant pair once conditioned on.

    A and B are drawn independently given depth — there is no pairing to find. The pooled test
    is fooled (that is the 45-49% false-positive regime measured on simulated null pairs);
    conditioning on depth must not be.
    """
    df = _depth_confounded_cohort()
    kw = dict(min_incidence=2, min_cooccurrence=2)
    pooled = _ab(cooccurrence(df, depth_strata=0, **kw))
    strat = _ab(cooccurrence(df, depth_strata=10, **kw))
    assert pooled.height == 1 and strat.height == 1

    # Depth inflates the pooled lift above 1 even though A ⫫ B given depth ...
    assert pooled["theta"][0] > 1.15, f"expected depth-induced lift, got {pooled['theta'][0]}"
    assert pooled["p_value"][0] < 0.05, "pooled test should be fooled by depth here"
    # ... and conditioning on depth removes both the significance and the lift.
    assert strat["p_value"][0] > 0.05, (
        f"CMH still significant (p={strat['p_value'][0]:.3g}, or_mh={strat['or_mh'][0]:.3f}) "
        "— depth conditioning failed")
    assert strat["or_mh"][0] < pooled["odds_ratio"][0]


def test_depth_conditioning_is_skipped_when_depth_barely_varies():
    """Conditioning is not free — below theta_depth ~1.2 there is nothing to correct.

    Guards the shallow/uniform-depth cohort: stratifying there sheds power (and on a shallow
    repertoire a single clonotype is a large share of depth, making it a mediator of the very
    pair under test), so the pooled test is kept.
    """
    df = _paired_cohort()                                  # CV(depth) ~ 0.30 -> theta_depth ~1.09
    strat = cooccurrence(df, min_cooccurrence=2, depth_strata=10)
    assert "or_mh" not in strat.columns, "should have fallen back to the pooled test"
    pooled = cooccurrence(df, min_cooccurrence=2, depth_strata=0)
    assert strat.equals(pooled)


def test_depth_strata_retains_a_genuine_pair_in_a_deep_cohort():
    """Calibration must not cost the real signal: a planted pair survives depth conditioning."""
    rng = np.random.default_rng(3)
    rows = []
    for i in range(400):                       # deep, skewed cohort -> conditioning is active
        sid = f"S{i:03d}"
        depth = max(int(np.exp(rng.normal(np.log(40), 0.9))), 2)
        rows.append(_c(sid, "TRAV2", "TRAJ2", "CAVBG"))
        rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSBG"))
        if i < 160:                            # a REAL pair: both chains of one clone, together
            rows.append(_c(sid, "TRAV1", "TRAJ1", "CAVRF"))
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSPF"))
        for f in range(depth):
            rows.append(_c(sid, "TRAV9", "TRAJ9", f"CAVF{i}_{f}"))
            rows.append(_c(sid, "TRBV9", "TRBJ9", f"CASF{i}_{f}"))
    hit = _ab(cooccurrence(pl.DataFrame(rows), min_cooccurrence=2, depth_strata=10))
    assert hit.height == 1
    assert hit["p_value"][0] < 1e-6, f"planted pair lost by stratification (p={hit['p_value'][0]})"
    assert hit["or_mh"][0] > 1


def _mutually_exclusive_cohort(seed=5, n=400, n_both=10):
    """A and B strongly ANTI-correlated, on a skewed-depth cohort so CMH engages.

    They must still co-occur in ``n_both`` subjects — far BELOW the ~n_a*n_b/n expected — so the
    pair clears ``min_cooccurrence`` and actually appears in the output. (A strictly alternating
    fixture gives n_ab=0, the pair is dropped, and the test silently asserts nothing.)
    """
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        sid = f"S{i:03d}"
        depth = max(int(np.exp(rng.normal(np.log(40), 0.9))), 2)
        rows.append(_c(sid, "TRAV2", "TRAJ2", "CAVBG"))
        rows.append(_c(sid, "TRBV2", "TRBJ2", "CASSBG"))
        if i < n_both:                              # the few overlaps: n_ab = n_both
            rows.append(_c(sid, "TRAV1", "TRAJ1", "CAVRF"))
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSPF"))
        elif i % 2:                                 # otherwise A xor B -> n_a ~ n_b ~ n/2
            rows.append(_c(sid, "TRAV1", "TRAJ1", "CAVRF"))
        else:
            rows.append(_c(sid, "TRBV1", "TRBJ1", "CASSPF"))
        for f in range(depth):
            rows.append(_c(sid, "TRAV9", "TRAJ9", f"CAVF{i}_{f}"))
            rows.append(_c(sid, "TRBV9", "TRBJ9", f"CASF{i}_{f}"))
    return pl.DataFrame(rows)


def test_cmh_honours_alternative_and_never_ranks_anti_correlated_pairs_as_cooccurring():
    """A mutually exclusive pair must NOT be a top co-occurrence hit under the CMH default.

    stats.cmh returns a two-sided chi2 p; handing it back unmodified made `alternative` inert, so
    an ANTI-correlated pair (or_mh -> 0) came back as the rank-1 "co-occurrence" hit at p=3e-78.
    Callers filter on q_value alone, so direction must be honoured inside the test.
    """
    df = _mutually_exclusive_cohort()
    r = cooccurrence(df, min_incidence=2, min_cooccurrence=1, depth_strata=10)
    assert "or_mh" in r.columns, "CMH must be engaged for this fixture (skewed depth)"
    hit = _ab(r)
    assert hit.height == 1, "the anti-correlated pair must be present to be asserted on"
    assert hit["theta"][0] < 0.5 and hit["or_mh"][0] < 1, "fixture should be anti-correlated"
    assert hit["p_value"][0] > 0.5, (
        f"anti-correlated pair significant under alternative='greater' "
        f"(p={hit['p_value'][0]:.3g}, or_mh={hit['or_mh'][0]:.3g})")
    # ... and it must not outrank genuine co-occurrence.
    top = r.sort("p_value").row(0, named=True)
    assert not (top["a_junction_aa"] == "CAVRF" and top["b_junction_aa"] == "CASSPF"), \
        "anti-correlated pair is the rank-1 co-occurrence hit"

    # 'less' finds it; 'greater' does not — the direction must actually do something.
    rl = _ab(cooccurrence(df, min_incidence=2, min_cooccurrence=1, depth_strata=10,
                          alternative="less"))
    rg = _ab(cooccurrence(df, min_incidence=2, min_cooccurrence=1, depth_strata=10,
                          alternative="greater"))
    if rl.height and rg.height:
        assert rl["p_value"][0] < rg["p_value"][0], "alternative= is inert under CMH"


def test_cooccurrence_is_reproducible_across_runs():
    """max_features cuts through an incidence tie band, so candidate order must be deterministic."""
    rng = np.random.default_rng(2)
    rows = []
    for i in range(60):                              # many features sharing one incidence value
        sid = f"S{i:02d}"
        for f in range(40):
            if rng.random() < 0.25:
                rows.append(_c(sid, "TRAV1", "TRAJ1", f"CAV{f:02d}"))
                rows.append(_c(sid, "TRBV1", "TRBJ1", f"CAS{f:02d}"))
    df = pl.DataFrame(rows)
    kw = dict(min_incidence=2, min_cooccurrence=2, max_features=15, depth_strata=0)
    runs = [cooccurrence(df, **kw) for _ in range(6)]
    assert len({r.height for r in runs}) == 1, f"tested counts vary: {[r.height for r in runs]}"
    for r in runs[1:]:
        assert r.equals(runs[0]), "identical inputs produced a different pair set"


def test_alternative_is_validated():
    with pytest.raises(ValueError, match="alternative must be"):
        cooccurrence(_paired_cohort(), alternative="nope")


def test_bincount_is_exact_vs_int64_matmul():
    """The BLAS route must be bitwise-identical to the int64 matmul it replaced.

    Entries are dot products of 0/1 over n_universe subjects, so they are integers well under
    2**53 and float64 carries them exactly. Pin that, since the whole reason to route through
    float64 is speed and a silent rounding would be a wrong biomarker count.
    """
    from vdjtools.biomarker.cooccurrence import _bincount

    rng = np.random.default_rng(0)
    for n_sub, n_feat, p in ((572, 200, 0.11), (50, 33, 0.5), (7, 3, 1.0)):
        m_a = rng.random((n_sub, n_feat)) < p
        m_b = rng.random((n_sub, n_feat)) < p
        got = _bincount(m_a, m_b)
        want = m_a.T.astype(np.int64) @ m_b.astype(np.int64)
        assert np.array_equal(got, want)
        assert got.dtype == want.dtype
        assert got.max() <= n_sub                      # the bound the exactness argument rests on
