"""Cross-validation of vdjtools basic analytics against mirpy oracle values.

Self-contained: it hard-codes mirpy's asserted oracle constants (each cited to the
mirpy test that asserts it) and checks that OUR functions reproduce them on the same
inputs. mirpy is deliberately NOT imported — it must remain a non-dependency of the
shipped test suite. Where our legacy-faithful definition intentionally differs from
mirpy, we assert OUR value and document why.

All CDR3 inputs are junction/anchored strings so the CDR3 convention lines up: mirpy's
``junction_aa`` equals our canonical ``cdr3_aa`` (both include the conserved anchors).
"""
import math

import numpy as np
import polars as pl

from vdjtools.io import schema as S
from vdjtools import stats, features as F, overlap as O


def _frame(cdr3, counts, v="TRBV12-3*01", j="TRBJ2-2*01"):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: [v] * n, S.J_CALL: [j] * n, S.CDR3_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


# --- diversity — mir tests/test_diversity.py -------------------------------

def test_diversity_matches_mirpy_summarize_counts():
    # mir::test_diversity_summary_core_metrics: summarize_counts([10,5,2,1,1,1]) ->
    #   abundance=20, diversity=6, singletons=3, doubletons=1, chao1=7.5.
    c = np.array([10, 5, 2, 1, 1, 1])
    assert int(c.sum()) == 20                              # mirpy abundance
    assert stats.observed_richness(c) == 6                 # mirpy diversity
    assert stats.chao1(c) == 7.5                           # mirpy chao1 (same formula)
    # mirpy hill q1 == our shannon_wiener (exp of Shannon entropy); mirpy hill q2 ==
    # our inverse_simpson. Constants computed exactly in Python from p=c/20.
    assert math.isclose(stats.inverse_simpson(c), 3.0303030303030303, rel_tol=1e-12)
    assert math.isclose(stats.shannon_wiener(c), 3.9462490923720126, rel_tol=1e-12)


def test_hill_matches_mirpy_reference_values():
    # mir::test_hill_curve_reference_values: hill_curve([8,4,2,1], q=[0,1,2]) ->
    #   q0=4, q2 = 1/Σp² = 225/85 ≈ 2.6470588 (p=[8,4,2,1]/15).
    c = np.array([8, 4, 2, 1])
    assert stats.observed_richness(c) == 4                 # mirpy hill q0
    assert math.isclose(stats.inverse_simpson(c), 225 / 85, rel_tol=1e-12)   # hill q2


# --- overlap — mir tests/test_overlap.py -----------------------------------

def test_overlap_identical_matches_mirpy():
    # mir::test_identical_reps_all_match: rep=[CASSF,CASSY,CASSW] vs itself ->
    #   f_similarity=1.0, f2_similarity=1.0. With 3 shared clonotypes R is defined
    #   (legacy guard n>2) and equals 1.0.
    rep = _frame(["CASSF", "CASSY", "CASSW"], [1, 2, 3])
    m = O.overlap_metrics(rep, rep, key=(S.CDR3_AA,))
    assert m["d12"] == 3
    assert math.isclose(m["F"], 1.0, rel_tol=1e-12)
    assert math.isclose(m["F2"], 1.0, rel_tol=1e-12)
    assert math.isclose(m["R"], 1.0, rel_tol=1e-9)


def test_overlap_partial_matches_mirpy_F_and_F2():
    # mir::test_partial_overlap: rep1=[CASSF(1),CASSY(2)] vs rep2=[CASSF(1),CASSW(2)]
    #   shared = {CASSF}; mirpy F=1/3, F2=1/3.
    r1 = _frame(["CASSF", "CASSY"], [1, 2])
    r2 = _frame(["CASSF", "CASSW"], [1, 2])
    m = O.overlap_metrics(r1, r2, key=(S.CDR3_AA,))
    assert m["d12"] == 1
    assert math.isclose(m["F"], 1 / 3, rel_tol=1e-12)
    assert math.isclose(m["F2"], 1 / 3, rel_tol=1e-12)
    # Only 1 shared clonotype -> R undefined. Legacy returns None here; mirpy returns
    # NaN (its correlation needs >=2 overlapping, ours needs >=3). Both "undefined".
    assert m["R"] is None


# --- INTENTIONAL differences: assert OURS (legacy), note where mirpy differs ---

def test_overlap_D_is_legacy_definition_not_mirpy():
    # ours D = d12/(d1*d2) (legacy OverlapEvaluator: div12/div1/div2);
    # mirpy d_similarity = n12/sqrt(n1*n2). For the partial case ours = 1/(2*2) = 0.25
    # (mir::test_partial_overlap asserts mirpy's 0.5). We assert OUR legacy value.
    r1 = _frame(["CASSF", "CASSY"], [1, 2])
    r2 = _frame(["CASSF", "CASSW"], [1, 2])
    assert O.overlap_metrics(r1, r2, key=(S.CDR3_AA,))["D"] == 0.25


def test_overlap_F2_and_R_on_disjoint_are_legacy_values():
    # ours F2 on disjoint = 0.0 (empty sum, legacy); mir::test_disjoint_reps_zero_overlap
    # asserts mirpy's NaN. ours R undefined -> None; mirpy NaN. Assert OUR values.
    r1 = _frame(["CASSF"], [1])
    r2 = _frame(["CASSY"], [1])
    m = O.overlap_metrics(r1, r2, key=(S.CDR3_AA,))
    assert m["d12"] == 0
    assert m["F2"] == 0.0
    assert m["R"] is None


# --- kmer — mir tests/test_tokens.py ---------------------------------------

def test_kmer_matches_mirpy_tokenize():
    # mirpy tokenize builds sliding, overlapping k-mers
    # (mir::test_tokens: tokenize("CASSL",3) == [CAS,ASS,SSL]). For CASSLAP, k=4 the
    # tokens are {CASS, ASSL, SSLA, SLAP}.
    prof = F.kmer_profile(_frame(["CASSLAP"], [1]), k=4, weight="unique", by_locus=False)
    assert set(prof["kmer"].to_list()) == {"CASS", "ASSL", "SSLA", "SLAP"}


# --- usage allele-stripping — mir tests/test_alleles.py --------------------

def test_usage_strips_allele_like_mirpy():
    # mirpy strip_allele("TRBV12-3*01") == "TRBV12-3" (default gene-level rollup).
    stripped = pl.DataFrame({"g": ["TRBV12-3*01"]}).select(S.strip_allele(pl.col("g")))["g"][0]
    assert stripped == "TRBV12-3"
    u = stats.segment_usage(_frame(["CASSF"], [1], v="TRBV12-3*01"), "v")
    assert u["v_call"].to_list() == ["TRBV12-3"]
