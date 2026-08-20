"""Tests for vdjtools.preprocess.filter — functional / frequency / segment / by-sample."""
import math

import polars as pl
import pytest

from vdjtools.io import schema as S
from vdjtools import preprocess as pp


def _sample(cdr3, counts, v=None, j=None, d=None):
    n = len(cdr3)
    df = pl.DataFrame({
        S.V_CALL: v or ["TRBV1"] * n, S.D_CALL: d or ["TRBD1"] * n,
        S.J_CALL: j or ["TRBJ1"] * n, S.JUNCTION_AA: cdr3, S.COUNT: counts,
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_filter_functional_drops_stop_and_oof():
    # CA*SF has a stop codon; CASaL carries a lowercase nucleotide (legacy out-of-frame
    # marker [atgc#~_?]); CASSL is clean coding.
    a = _sample(["CASSL", "CA*SF", "CASaL"], [10, 20, 30])
    coding = pp.filter_functional(a, keep="coding")
    assert coding[S.JUNCTION_AA].to_list() == ["CASSL"]
    assert math.isclose(coding[S.FREQ].sum(), 1.0, rel_tol=1e-12)
    noncoding = pp.filter_functional(a, keep="noncoding")
    assert set(noncoding[S.JUNCTION_AA].to_list()) == {"CA*SF", "CASaL"}


def test_filter_frequency_min_freq():
    a = _sample(["A", "B", "C"], [1, 9, 90])                 # freqs .01 .09 .90
    out = pp.filter_frequency(a, min_freq=0.05)
    assert set(out[S.JUNCTION_AA].to_list()) == {"B", "C"}       # drops the .01 clone


def test_filter_frequency_top_quantile():
    # counts 60,30,10 -> freqs .6,.3,.1. top_quantile=0.25 keeps only clones whose
    # cumulative read mass stays <= 0.25 -> just the top (0.6 > 0.25 excludes rest,
    # but the top itself: cumulative 0.6 > 0.25 -> even it fails? legacy keeps while
    # cumulative <= q; here 0.6 > 0.25 so nothing passes). Use q=0.65 -> keep top only.
    a = _sample(["T", "M", "S"], [60, 30, 10])
    q = pp.filter_frequency(a, top_quantile=0.65)
    assert q[S.JUNCTION_AA].to_list() == ["T"]                   # cum .6<=.65, +.3 -> .9>.65
    q2 = pp.filter_frequency(a, top_quantile=0.95)
    assert q2[S.JUNCTION_AA].to_list() == ["T", "M"]             # .6,.9<=.95 ; +.1 excluded


def test_filter_segment_keep_and_remove():
    a = _sample(["A", "B", "C"], [1, 1, 1],
                v=["TRBV12-3*01", "TRBV20-1*01", "TRBV12-4*01"])
    keep = pp.filter_segment(a, v=["TRBV12"])               # prefix, allele-insensitive
    assert set(keep[S.V_CALL].to_list()) == {"TRBV12-3*01", "TRBV12-4*01"}
    remove = pp.filter_segment(a, v=["TRBV12"], keep=False)
    assert remove[S.V_CALL].to_list() == ["TRBV20-1*01"]


def test_filter_segment_vj_conjunction():
    a = _sample(["A", "B"], [1, 1], v=["TRBV1", "TRBV1"], j=["TRBJ1-1", "TRBJ2-1"])
    out = pp.filter_segment(a, v=["TRBV1"], j=["TRBJ2"])    # V AND J must match
    assert out[S.J_CALL].to_list() == ["TRBJ2-1"]


def test_filter_by_sample_keep_and_remove():
    a = _sample(["CASSL", "CASSF", "CASSX"], [1, 2, 3])
    other = _sample(["CASSL", "CASSF"], [5, 5])
    keep = pp.filter_by_sample(a, other, keep=True)
    assert set(keep[S.JUNCTION_AA].to_list()) == {"CASSL", "CASSF"}
    remove = pp.filter_by_sample(a, other, keep=False)
    assert remove[S.JUNCTION_AA].to_list() == ["CASSX"]


def test_filter_by_sample_key_includes_vj():
    a = _sample(["CASSL"], [1], v=["TRBV1"], j=["TRBJ1"])
    other = _sample(["CASSL"], [1], v=["TRBV9"], j=["TRBJ1"])   # same aa, different V
    assert pp.filter_by_sample(a, other, keep=True).height == 0  # default key has V/J
    assert pp.filter_by_sample(a, other, keep=True, key=(S.JUNCTION_AA,)).height == 1


def test_filter_segment_matches_a_gene_named_only_in_a_tie():
    """A comma-ambiguity tie must match a query for its non-first gene."""
    from vdjtools.preprocess.filter import filter_segment

    df = pl.DataFrame({
        "v_call": ["IGHV3-23*01,IGHV3-23D*01", "IGHV1-2*02"],
        "d_call": [None, None], "j_call": ["IGHJ4*02", "IGHJ4*02"],
        "junction_aa": ["CASF", "CASG"], "junction_nt": ["ACG", "ACG"],
        "duplicate_count": [1, 1], "frequency": [0.5, 0.5],
    })
    kept = filter_segment(df, v=["IGHV3-23D"])
    assert kept.height == 1 and kept["v_call"][0].startswith("IGHV3-23*01,")
    # keep=False removes it (the legacy --negative path must also see the tie)
    dropped = filter_segment(df, v=["IGHV3-23D"], keep=False)
    assert dropped.height == 1 and dropped["v_call"][0] == "IGHV1-2*02"


# ------------------------------------------------- the three axes, kept apart (2026-08-20)


def _f(**kw):
    base = dict(junction_aa=["CASSIRSSYEQYF", "CASS*YEQYF", "CASS_YEQYF"],
                duplicate_count=[100, 50, 50], frequency=[0.5, 0.25, 0.25])
    base.update(kw)
    return pl.DataFrame(base)


class TestFilterProductive:
    def test_derives_from_junction_aa_when_no_airr_columns(self):
        from vdjtools.preprocess import filter_productive, productive_mask
        _, src = productive_mask(_f())
        assert src == "junction_aa"
        assert filter_productive(_f())["junction_aa"].to_list() == ["CASSIRSSYEQYF"]

    def test_the_files_productive_column_wins_over_the_junction(self):
        """The whole point of reading the annotation: the file may disagree, and it is right."""
        from vdjtools.preprocess import filter_productive, productive_mask
        d = _f(productive=[False, True, True])
        _, src = productive_mask(d)
        assert src == "productive"
        assert filter_productive(d)["junction_aa"].to_list() == ["CASS*YEQYF", "CASS_YEQYF"]

    def test_components_used_when_the_composite_is_absent(self):
        from vdjtools.preprocess import filter_productive, productive_mask
        d = _f(stop_codon=[False, True, False], vj_in_frame=[True, True, False])
        _, src = productive_mask(d)
        assert src == "stop_codon+vj_in_frame"
        assert filter_productive(d)["junction_aa"].to_list() == ["CASSIRSSYEQYF"]

    def test_recompute_frequencies_is_a_switch(self):
        from vdjtools.preprocess import filter_productive
        assert filter_productive(_f())["frequency"].to_list() == [1.0]
        assert filter_productive(_f(), recompute_frequencies=False)["frequency"].to_list() == [0.5]

    def test_nonproductive_is_the_complement(self):
        from vdjtools.preprocess import filter_productive
        assert filter_productive(_f(), keep="nonproductive")["junction_aa"].to_list() == [
            "CASS*YEQYF", "CASS_YEQYF"]

    def test_legacy_alias_warns_but_works(self):
        from vdjtools.preprocess import filter_functional
        with pytest.warns(DeprecationWarning, match="filter_productive"):
            assert filter_functional(_f(), keep="coding").height == 1


class TestFilterLength:
    def test_bounds_are_inclusive(self):
        from vdjtools.preprocess import filter_length
        d = pl.DataFrame({"junction_aa": ["C" * 4, "C" * 5, "C" * 60, "C" * 61],
                          "duplicate_count": [1] * 4, "frequency": [0.25] * 4})
        assert filter_length(d)["junction_aa"].str.len_chars().to_list() == [5, 60]

    def test_outside_is_the_complement_and_takes_the_nulls(self):
        from vdjtools.preprocess import filter_length
        d = pl.DataFrame({"junction_aa": ["C" * 4, "C" * 10, None],
                          "duplicate_count": [1] * 3, "frequency": [0.3] * 3})
        assert filter_length(d, keep="outside").height == 2

    def test_min_above_max_raises(self):
        from vdjtools.preprocess import filter_length
        with pytest.raises(ValueError, match="exceeds"):
            filter_length(_f(), min_len=30, max_len=10)


class TestFilterFunctionalGenes:
    def test_an_unresolvable_gene_is_kept_not_dropped(self):
        """A name we do not recognise is a vocabulary gap, never evidence of a pseudogene."""
        from vdjtools.preprocess import filter_functional_genes
        d = pl.DataFrame({"v_call": ["TRBV20-1", "NOT-A-REAL-GENE"], "j_call": ["TRBJ2-2"] * 2,
                          "junction_aa": ["CASSIRSSYEQYF"] * 2,
                          "duplicate_count": [10, 10], "frequency": [0.5, 0.5]})
        assert filter_functional_genes(d, locus="TRB").height == 2
