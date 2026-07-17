"""Tests for vdjtools.features.kmer — CDR3 k-mer and V+kmer+C profiles."""
import polars as pl

from vdjtools.io import schema as S
from vdjtools import features as F


def _frame():
    df = pl.DataFrame({
        S.V_CALL: ["TRBV1*01", "TRBV1*01"],
        S.J_CALL: ["TRBJ1", "TRBJ1"],
        S.C_CALL: ["IGHM", None],
        S.JUNCTION_AA: ["CASSL", "CAS"],
        S.COUNT: [3, 2],
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_kmer_profile_reads():
    prof = F.kmer_profile(_frame(), k=2, weight="reads")
    got = dict(zip(prof["kmer"].to_list(), prof["weight"].to_list()))
    # CASSL -> CA,AS,SS,SL (each x3); CAS -> CA,AS (each x2)
    assert got == {"CA": 5, "AS": 5, "SS": 3, "SL": 3}


def test_kmer_profile_unique_short_cdr3_skipped():
    df = pl.DataFrame({
        S.V_CALL: ["TRBV1"], S.J_CALL: ["TRBJ1"], S.C_CALL: [None],
        S.JUNCTION_AA: ["AC"], S.COUNT: [1],
    })
    df = S.add_locus(S.normalize(df, recompute_freq=True))
    prof = F.kmer_profile(df, k=3, weight="unique")     # len 2 < k -> no kmers
    assert prof.height == 0


def test_v_kmer_c_profile_null_c_retained():
    prof = F.v_kmer_c_profile(_frame(), k=2, weight="unique")
    # the null-c clonotype (CAS) contributes CA/AS rows with c_call null
    null_c = prof.filter(pl.col("c_call").is_null())
    assert set(null_c["kmer"].to_list()) == {"CA", "AS"}
    assert prof["v_call"].unique().to_list() == ["TRBV1"]   # allele stripped


# --- anchor-aware core k-mers ------------------------------------------------------------------

def test_flank_matches_seqtree_core_kmers():
    """`flank` must agree with seqtree.seeds.core_kmers exactly — it is the same definition, and
    seqtree's measured background (central 4-mer 0.080% of control vs CASS 56.5%) only applies if
    the cores are identical. Two edge cases this pins: `str.len_chars()` is UInt32 so
    `len - 2*flank` UNDERFLOWS on a short junction, and polars' str.slice with a bad length
    silently returns a germline tail instead of nothing.
    """
    import random
    core_kmers = __import__("seqtree.seeds", fromlist=["core_kmers"]).core_kmers
    from vdjtools.features.kmer import _explode_kmers
    random.seed(0)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    seqs = sorted({"CASSVGLYSTDTQYF", "CASSVGGFGDTQYF", "CASSFGH", "CAF", "C", "CA",
                   *("C" + "".join(random.choice(aa) for _ in range(n)) + "F" for n in range(22))})
    df = pl.DataFrame({S.JUNCTION_AA: seqs})
    for k in range(2, 8):
        for flank in (0, 1, 3, 4, 6):
            got = {d[S.JUNCTION_AA]: set(d["kmer"]) for d in
                   _explode_kmers(df, k, flank).group_by(S.JUNCTION_AA)
                   .agg(pl.col("kmer").unique()).to_dicts()}
            for s in seqs:
                assert set(core_kmers(s, k, flank=flank)) == got.get(s, set()), \
                    f"k={k} flank={flank} {s!r}"


def test_flank_drops_the_anchors():
    df = pl.DataFrame({S.JUNCTION_AA: ["CASSVGLYSTDTQYF"]})
    from vdjtools.features.kmer import _explode_kmers
    assert set(_explode_kmers(df, 4, 4)["kmer"]) == {"VGLY", "GLYS", "LYST", "YSTD"}  # CASS/TQYF gone
    assert "CASS" in set(_explode_kmers(df, 4, 0)["kmer"])                            # flank=0 keeps them


def test_kmer_cohort_feeds_association():
    """kmer_cohort -> association(key=('v_call','kmer'), match='exact') is the V+kmer test."""
    from vdjtools.biomarker import association
    from vdjtools.features.kmer import kmer_cohort
    n = 8
    df = pl.DataFrame({
        "sample_id": [f"s{i}" for i in range(n)],
        S.V_CALL: ["TRBV9*01"] * n,
        S.J_CALL: ["TRBJ2-3"] * n,
        # the first 4 carry the VGLY core, the last 4 do not
        S.JUNCTION_AA: ["CASSVGLYSTDTQYF"] * 4 + ["CASSQQQQSTDTQYF"] * 4,
        S.COUNT: [1] * n,
    })
    km = kmer_cohort(df, k=4, flank=4)
    assert set(km.columns) == {"sample_id", S.V_CALL, "kmer"}
    design = pl.DataFrame({"sample_id": [f"s{i}" for i in range(n)],
                           "case": [1] * 4 + [0] * 4})
    res = association(km, design, pheno_col="case", key=(S.V_CALL, "kmer"),
                      match="exact", min_incidence=1)
    vgly = res.filter(pl.col("kmer") == "VGLY")
    assert vgly.height == 1
    assert vgly["n_pos_present"][0] == 4 and vgly["n_neg_present"][0] == 0


def test_association_rejects_fuzzy_without_junction_aa():
    from vdjtools.biomarker import association
    df = pl.DataFrame({"sample_id": ["a"], S.V_CALL: ["TRBV9"], "kmer": ["VGLY"]})
    design = pl.DataFrame({"sample_id": ["a"], "case": [1]})
    import pytest
    with pytest.raises(ValueError, match="searches on"):
        association(df, design, pheno_col="case", key=(S.V_CALL, "kmer"), match="fuzzy")
    with pytest.raises(ValueError, match="absent from the cohort"):
        association(df, design, pheno_col="case", key=(S.V_CALL, "nope"), match="exact")
