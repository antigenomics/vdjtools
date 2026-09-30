"""The junction pipeline: `(junction_aa, V, J)` -> repair -> nucleotides -> D.

What is pinned here is the **contract and the coordinate algebra**, not the accuracy: accuracy is
measured against real nucleotide rearrangements from `isalgo/airr_control` by
`appendix/bench_junction_pipeline.py`, and recorded in the CHANGELOG and in
`docs/junction_pipeline.md`.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import JUNCTION_COLUMNS, annotate_junctions, load_bundled
from vdjtools.model.junction import _aa_span, _reachable

TRB = (["CASSLAPGATNEKLFF", "CASSPGQGAYEQYF", "CAWSVSDLAKNIQYF"],
       ["TRBV5-1*01", "TRBV5-1*01", "TRBV30*01"],
       ["TRBJ1-4*01", "TRBJ2-7*01", "TRBJ2-4*01"])


def test_one_row_out_per_row_in_in_input_order():
    """The contract a consumer joins on. A row nothing can explain is present with nulls."""
    cdr3 = [*TRB[0], "CAVTDDKIIF", "GARBAGE", ""]
    v = [*TRB[1], "TRAV12-2*01", "TRBV9*01", "TRBV9*01"]
    j = [*TRB[2], "TRAJ30*01", "TRBJ1-1*01", "TRBJ1-1*01"]
    out = annotate_junctions(cdr3, v, j, species="human")
    assert out.height == len(cdr3)
    assert out["cdr3_aa"].to_list() == cdr3, "input order, row for row"
    assert set(JUNCTION_COLUMNS) >= set(out.columns)
    assert out["d_call"][3] is None, "a VJ locus has no D"


def test_ragged_input_raises_rather_than_truncating():
    with pytest.raises(ValueError, match="ragged"):
        annotate_junctions(["CASSLAPGATNEKLFF"], [], [])
    with pytest.raises(ValueError, match="ragged"):
        annotate_junctions(TRB[0], TRB[1], TRB[2], species=["human"])
    assert annotate_junctions([], [], []).height == 0


def test_the_nucleotide_junction_is_exactly_three_times_the_repaired_one():
    """Stage 2's output has to be in the same coordinate frame as stage 1's, or stage 3 is nonsense."""
    out = annotate_junctions(*TRB, species="human").filter(pl.col("cdr3_nt").is_not_null())
    assert out.height
    for nt, aa in zip(out["cdr3_nt"], out["cdr3_repaired"]):
        assert len(nt) == 3 * len(aa), (nt, aa)


def test_the_d_sits_inside_the_boundaries_stage_one_placed():
    """Stage 3 is given the interior, not left to re-derive it from an inferred sequence.

    Re-deriving it matched germline against the CALLED allele while the nucleotides came from the
    model's representative allele of the same gene, so the prefix broke at the first synonymous
    difference and a spurious D won inside the V -- `CATSIRFTDTQYF` placed a TRBD2 at nucleotide 6.
    """
    out = annotate_junctions(*TRB, species="human").filter(pl.col("d_call").is_not_null())
    assert out.height
    for r in out.iter_rows(named=True):
        assert r["d_start_nt"] > r["v_end_nt"], r
        assert r["d_end_nt"] <= r["j_start_nt"], r


@pytest.mark.parametrize("start_nt, end_nt, want", [
    (4, 6, (1, 1)),        # exactly one codon
    (5, 7, (1, 2)),        # straddles a boundary: two residues, one codon's worth of nucleotides
    (1, 3, (0, 0)),
    (13, 19, (4, 6)),
    (-1, -1, (None, None)),
])
def test_a_nucleotide_span_folds_onto_the_residues_whose_codons_it_touches(start_nt, end_nt, want):
    """The point of the nucleotide detour: a 1-residue D can still report a 2-residue span.

    It is part of both flanking codons, and those codons are exactly the evidence a translated
    junction throws away.
    """
    assert _aa_span(start_nt, end_nt) == want


def test_flooring_the_model_creates_answers_and_changes_none():
    """`_reachable`, and the invariant that makes it safe.

    36 of 66 human TRB V alleles come out of the EM fit at p = 0, and 34 have an all-zero deletion
    profile, so conditioning on one returned no scenario at all. Flooring is only allowed because
    this stage pins the V and the J, which makes their usage a constant factor per row: measured over
    every allele the fit DID see, not one answer moves.
    """
    from vdjtools.model import infer_nt

    m = load_bundled("TRB", source="arda", organism="human")
    f = _reachable(m)
    t = m.tables["v_choice"]
    col = t.columns[0]
    fitted = [r[col] for r in t.filter(pl.col("p") > 0).to_dicts()]
    zeroed = [r[col] for r in t.filter(pl.col("p") <= 0).to_dicts()]
    assert len(zeroed) > 10, "the premise: the shipped fit really does zero many alleles"

    moved = 0
    for allele in fitted:
        a = infer_nt(m, "CASSTQENTEAFF", v=allele, j="TRBJ1-1*01")
        b = infer_nt(f, "CASSTQENTEAFF", v=allele, j="TRBJ1-1*01")
        if a is not None and b is not None and (a.cdr3_nt, a.d_call) != (b.cdr3_nt, b.d_call):
            moved += 1
    assert moved == 0, f"{moved} fitted alleles changed answer -- flooring is not neutral"
    assert sum(infer_nt(f, "CASSTQENTEAFF", v=x, j="TRBJ1-1*01") is not None
               for x in zeroed[:8]) >= 6, "zero-usage alleles must become answerable"


def test_the_model_names_the_d_and_the_alignment_places_it():
    """One gene estimator, one placer, and the bounds they agree on.

    Letting the aligner choose the gene as well was measured and dropped: on 4,000 real human TRB
    rearrangements the model's posterior is right on 74.30 % of ALL rows against the gated
    alignment's 47.93 %, it is no worse where the alignment is confident (86.32 % against 85.96 %
    on those 2,230 rows), and a hybrid of the two scores below the model alone. What the alignment
    is for is saying WHERE, ungated, so a row has coordinates to draw -- 99.70 % of rows do.
    """
    out = annotate_junctions(*TRB, species="human")
    placed = 0
    for r in out.iter_rows(named=True):
        if not r["d_call"]:
            assert r["d_start_nt"] is None and r["d_posterior"] is None
            continue
        assert "*" not in r["d_call"], "the D is reported as a GENE; the allele is not identifiable"
        assert 0.0 < r["d_posterior"] <= 1.0
        if r["d_start_nt"] is None:
            continue
        placed += 1
        # The D sits inside the junction, the way round it was written, and its residues are the
        # ones its codons touch -- one definition of where it is, not two that can disagree.
        assert 1 <= r["d_start_nt"] <= r["d_end_nt"] <= len(r["cdr3_nt"])
        assert r["d_start_aa"] == (r["d_start_nt"] - 1) // 3
        assert r["d_end_aa"] == (r["d_end_nt"] - 1) // 3
        assert r["cdr3_nt"][r["v_end_nt"]:r["d_start_nt"] - 1] == (r["np1"] or "")
        assert r["cdr3_nt"][r["d_end_nt"]:r["j_start_nt"]] == (r["np2"] or "")
    assert placed, "no D was placed at all"


def test_the_retired_d_columns_are_gone():
    """`d_best` existed only to reconcile two gene estimators, and there is one now."""
    for gone in ("d_best", "d_best_source", "d_posterior_call", "d_entropy", "d_support",
                 "d2_call", "np3"):
        assert gone not in JUNCTION_COLUMNS
    import vdjtools.model as M
    assert not hasattr(M, "posterior_d"), "the dominated posterior must not ship"


def test_auto_is_a_chain_and_answers_at_least_as_many_rows_as_either_set():
    """`model_source="auto"` runs OLGA's fit then arda's on what it left, and that is measured.

    Neither set dominates: OLGA's is better calibrated and cheaper but declines rows outright (143
    of 4,000 on human TRA) and has no mouse; arda's answers everything and is thinner -- 36 of its
    66 human TRB V alleles sit at probability zero. On 4,000 real human rearrangements the chain
    wins every column: TRB nucleotide-exact 14.40 / 17.32 / 17.32 % and D gene 72.58 / 74.08 /
    74.35 % for arda / OLGA / the chain.

    What is pinned here is the structural property that makes that possible -- the chain never
    answers FEWER rows than a set pinned by name -- not the accuracy numbers, which need the truth
    set the benchmark uses.
    """
    auto = annotate_junctions(*TRB, species="human")                      # the default
    answered = lambda f: f["cdr3_nt"].is_not_null().sum()                 # noqa: E731
    for pinned in ("olga", "arda"):
        one = annotate_junctions(*TRB, species="human", model_source=pinned)
        assert answered(auto) >= answered(one), (
            f"the chain answered {answered(auto)} rows, {pinned} alone answered {answered(one)}")
    assert answered(auto) > 0


def test_a_pinned_set_that_does_not_exist_is_not_silently_ignored():
    """Mouse has no OLGA fit, so `auto` must reach arda's -- and a bad name must answer nothing
    rather than quietly falling back to whatever does exist."""
    mouse = (["CASSLAPGATNEKLFF"], ["TRBV13-1*01"], ["TRBJ1-4*01"])
    assert annotate_junctions(*mouse, species="mouse")["cdr3_nt"][0] is not None
    assert annotate_junctions(*mouse, species="mouse", model_source="nosuchset")["cdr3_nt"][0] is None
