"""Tests for vdjtools.features.kmer — CDR3 k-mer and V+kmer+C profiles."""
import polars as pl

from vdjtools.io import schema as S
from vdjtools import features as F


def _frame():
    df = pl.DataFrame({
        S.V_CALL: ["TRBV1*01", "TRBV1*01"],
        S.J_CALL: ["TRBJ1", "TRBJ1"],
        S.C_CALL: ["IGHM", None],
        S.CDR3_AA: ["CASSL", "CAS"],
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
        S.CDR3_AA: ["AC"], S.COUNT: [1],
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
