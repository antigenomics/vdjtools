"""Tests for vdjtools.stats.usage — segment and V-J usage profiles."""
import polars as pl

from vdjtools.io import schema as S
from vdjtools import stats


def _frame():
    df = pl.DataFrame({
        S.V_CALL: ["TRBV12-3*01", "TRBV12-3*02", "TRBV20-1*01", "TRBV20-1*01"],
        S.D_CALL: [None, None, "TRBD1", None],
        S.J_CALL: ["TRBJ1-1*01", "TRBJ2-1*01", "TRBJ1-1*01", "TRBJ2-1*01"],
        S.C_CALL: [None, None, None, None],
        S.JUNCTION_AA: ["CASSA", "CASSB", "CASSC", "CASSD"],
        S.COUNT: [10, 20, 30, 40],
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_segment_usage_reads_strips_allele():
    u = stats.segment_usage(_frame(), "v", weight="reads")
    got = dict(zip(u["v_call"].to_list(), u["weight"].to_list()))
    # alleles collapsed: TRBV12-3 = 10+20 = 30 reads; TRBV20-1 = 30+40 = 70
    assert got == {"TRBV12-3": 30, "TRBV20-1": 70}


def test_segment_usage_unique_and_keep_allele():
    u = stats.segment_usage(_frame(), "v", weight="unique", keep_allele=True)
    got = dict(zip(u["v_call"].to_list(), u["weight"].to_list()))
    assert got == {"TRBV12-3*01": 1, "TRBV12-3*02": 1, "TRBV20-1*01": 2}


def test_segment_usage_empty_when_all_null():
    u = stats.segment_usage(_frame(), "c", weight="reads")
    assert u.height == 0


def test_vj_usage_pairs():
    vj = stats.vj_usage(_frame(), weight="unique")
    pairs = {(r["v_call"], r["j_call"]): r["weight"] for r in vj.iter_rows(named=True)}
    assert pairs[("TRBV12-3", "TRBJ1-1")] == 1
    assert pairs[("TRBV20-1", "TRBJ1-1")] == 1
    assert pairs[("TRBV20-1", "TRBJ2-1")] == 1
    assert "locus" in vj.columns


def test_strip_allele_never_drops_a_gene_from_a_tie():
    """A comma-ambiguous call keeps every distinct gene; a same-gene allele tie collapses to one.

    Regression: the old greedy `\\*.*$` matched from the first `*` to end and silently dropped
    IGHV3-23D from "IGHV3-23*01,IGHV3-23D*01", so it reported zero usage cohort-wide.
    """
    from vdjtools.io.schema import strip_allele

    df = pl.DataFrame({"g": ["TRBV12-3*01", "IGHV3-23*01,IGHV3-23D*01",
                             "IGHV1-2*02,IGHV1-2*04", "TRBV20-1", None]})
    out = df.select(strip_allele(pl.col("g")).alias("g"))["g"].to_list()
    assert out[0] == "TRBV12-3"
    assert out[1] == "IGHV3-23,IGHV3-23D"          # both genes kept
    assert out[2] == "IGHV1-2"                      # same-gene allele tie -> single gene
    assert out[3] == "TRBV20-1"
    assert out[4] is None


def test_segment_usage_counts_the_second_gene_of_a_tie():
    """IGHV3-23D must not vanish from segment usage when it appears only in ties."""
    from vdjtools.stats import segment_usage

    df = pl.DataFrame({
        "v_call": ["IGHV3-23*01,IGHV3-23D*01"] * 3 + ["IGHV1-2*02"] * 2,
        "j_call": ["IGHJ4*02"] * 5, "junction_aa": [f"CAS{i}F" for i in range(5)],
        "junction_nt": ["ACG"] * 5, "duplicate_count": [1] * 5, "frequency": [0.2] * 5,
    })
    u = segment_usage(df, "v", weight="unique")
    genes = set(u["v_call"].to_list())
    assert any("IGHV3-23D" in g for g in genes), f"IGHV3-23D dropped from usage: {genes}"
