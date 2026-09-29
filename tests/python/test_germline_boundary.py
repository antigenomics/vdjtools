"""`germline_boundary`: the V/J boundary a germline alignment supports, in nucleotides.

`Scenario.v_end` is the argmax *history*'s boundary and is the wrong statistic to read a boundary
off (#182). These pin the replacement, and the two properties that make it safe: it uses the
model's own germline rather than arda's, and it does not need the DP to have succeeded.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import germline_boundary, infer_nt_batch, load_bundled

# (junction, V, J, v_end nt, j_start nt). The first two carry external nucleotide truth from
# `isalgo/airr_control`: row 1's V truth is 14 -- one of the 19.6 % where a synonymous codon inside
# the templated run breaks the nucleotide match while leaving the amino-acid match intact, which no
# amino-acid evidence can see -- and every other value here is exact.
CASES = [
    ("CATSDPGTGHQPQHF", "TRBV24-1", "TRBJ1-5", 15, 28),
    ("CASSQVGTGVYEQYF", "TRBV3-1", "TRBJ2-7", 16, 29),
]


@pytest.fixture(scope="module")
def trb():
    return load_bundled("TRB")


@pytest.mark.parametrize("cdr3, v, j, v_end, j_start", CASES)
def test_the_boundary_is_where_the_germline_stops(trb, cdr3, v, j, v_end, j_start):
    got = germline_boundary(trb, [cdr3], [v], [j])
    assert (got["v_end"][0], got["j_start"][0]) == (v_end, j_start)


def test_one_row_per_input_row_in_input_order(trb):
    cdr3s = [c[0] for c in CASES] + ["", "CASSLAPGATNEKLFF"]
    got = germline_boundary(trb, cdr3s, [c[1] for c in CASES] + ["TRBV3-1", "TRBV5-1"],
                            [c[2] for c in CASES] + ["TRBJ2-7", "TRBJ1-4"])
    assert got.height == len(cdr3s)
    assert got["v_end"].to_list()[:2] == [c[3] for c in CASES]
    assert got["v_end"][2] is None and got["j_start"][2] is None       # nothing to align


def test_a_gene_name_an_allele_and_an_ambiguous_call_agree(trb):
    cdr3 = CASES[0][0]
    one = germline_boundary(trb, [cdr3], ["TRBV24-1"], ["TRBJ1-5"])
    two = germline_boundary(trb, [cdr3], ["TRBV24-1*01"], ["TRBJ1-5*01"])
    three = germline_boundary(trb, [cdr3], ["TRBV24-1,TRBV7-9"], ["TRBJ1-5,TRBJ2-7"])
    assert one.equals(two) and one.equals(three)


def test_no_call_declines_that_side_rather_than_guessing(trb):
    got = germline_boundary(trb, [CASES[0][0]], None, [CASES[0][2]])
    assert got["v_end"][0] is None
    assert got["j_start"][0] == CASES[0][4]


def test_a_call_the_model_does_not_carry_declines(trb):
    # `TRBV21-1` is a pseudogene the olga model has no germline for, and `*03` is an allele it
    # does not carry. A boundary from some other allele's germline would be a guess.
    got = germline_boundary(trb, [CASES[0][0]] * 2, ["TRBV21-1", "TRBV7-9*03"], ["TRBJ1-5"] * 2)
    assert got["v_end"].to_list() == [None, None]


@pytest.mark.parametrize("bad", ["v", "j"])
def test_a_call_column_of_the_wrong_length_raises(trb, bad):
    kw = {"v": ["TRBV24-1"], "j": ["TRBJ1-5"]}
    kw[bad] = kw[bad] * 2
    with pytest.raises(ValueError, match="same length"):
        germline_boundary(trb, [CASES[0][0]], **kw)


def test_the_answer_never_exceeds_the_model_s_own_germline(trb):
    """And so cannot have come from `cut_segment`, whose palindromic nt are not germline."""
    v_len = dict(zip(trb.genomic["genes_v"]["v_allele"], trb.genomic["genes_v"]["cdr3_segment"]))
    j_len = dict(zip(trb.genomic["genes_j"]["j_allele"], trb.genomic["genes_j"]["cdr3_segment"]))
    cdr3s = [c[0] for c in CASES]
    got = germline_boundary(trb, cdr3s, ["TRBV24-1*01", "TRBV3-1*01"],
                            ["TRBJ1-5*01", "TRBJ2-7*01"])
    for i, (cdr3, va, ja) in enumerate(zip(cdr3s, ["TRBV24-1*01", "TRBV3-1*01"],
                                           ["TRBJ1-5*01", "TRBJ2-7*01"])):
        assert got["v_end"][i] <= len(v_len[va])
        assert 3 * len(cdr3) - got["j_start"][i] <= len(j_len[ja])


def test_infer_nt_batch_carries_the_boundary_beside_the_history(trb):
    cdr3s = [c[0] for c in CASES]
    df = infer_nt_batch(trb, cdr3s, [c[1] for c in CASES], [c[2] for c in CASES])
    assert {"v_end_germline", "j_start_germline"} <= set(df.columns)
    alone = germline_boundary(trb, cdr3s, [c[1] for c in CASES], [c[2] for c in CASES])
    assert df["v_end_germline"].to_list() == alone["v_end"].to_list()
    assert df["j_start_germline"].to_list() == alone["j_start"].to_list()


def test_a_row_the_dp_declined_keeps_its_germline_boundary(trb):
    """The germline explains what it explains whether or not a recombination was found."""
    df = infer_nt_batch(trb, [CASES[0][0], "CASSXQVGTGVYEQYF"], ["TRBV24-1", "TRBV3-1"],
                        ["TRBJ1-5", "TRBJ2-7"])
    assert df["cdr3_nt"][1] is None and df["v_end"][1] is None
    assert df["v_end_germline"][1] is not None
    assert df["j_start_germline"][1] is not None


def test_the_boundary_does_not_depend_on_the_history_being_the_argmax(trb):
    """`n_best` changes which nt sequence wins; the germline boundary cannot move with it."""
    cdr3s = [c[0] for c in CASES]
    a = infer_nt_batch(trb, cdr3s, [c[1] for c in CASES], [c[2] for c in CASES], n_best=1)
    b = infer_nt_batch(trb, cdr3s, [c[1] for c in CASES], [c[2] for c in CASES], n_best=8)
    for col in ("v_end_germline", "j_start_germline"):
        assert a[col].to_list() == b[col].to_list()


def test_a_non_functional_allele_has_no_germline_and_declines(trb):
    empty = [r["v_allele"] for r in trb.genomic["genes_v"].iter_rows(named=True)
             if not r["cdr3_segment"]]
    assert empty, "the olga TRB model ships pseudogenes with no CDR3-region germline"
    got = germline_boundary(trb, [CASES[0][0]], [empty[0]], [CASES[0][2]])
    assert got["v_end"][0] is None


def test_the_frame_is_typed_even_when_every_row_declines(trb):
    got = germline_boundary(trb, ["", ""], None, None)
    assert got.schema == {"v_end": pl.Int64, "j_start": pl.Int64}
