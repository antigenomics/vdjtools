"""The junction pipeline: `(junction_aa, V, J)` -> repair -> nucleotides -> D.

What is pinned here is the **contract and the coordinate algebra**, not the accuracy: accuracy is
measured against real nucleotide rearrangements from `isalgo/airr_control` by
`appendix/bench_junction_pipeline.py`.
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


def test_the_nucleotides_always_translate_to_the_junction_asked_about():
    """#186/#187: a germline flank must never overwrite a residue; aa <> nt always agree.

    3,000 junctions drawn from the model: 20 of them had a residue overwritten before the fix.
    """
    from vdjtools.model.generate import generate
    from vdjtools.model.reference import translate
    g = generate(load_bundled("TRB"), 3000, productive_only=True).filter(
        pl.col("junction_aa").str.len_chars() > 5)
    out = annotate_junctions(g["junction_aa"].to_list(), g["v_call"].to_list(),
                             g["j_call"].to_list()).filter(pl.col("cdr3_nt").is_not_null())
    assert out.height > 2000
    bad = [(n, a) for n, a in zip(out["cdr3_nt"], out["cdr3_repaired"]) if translate(n) != a]
    assert not bad, bad[:3]


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


def test_auto_routes_to_isalgo_and_matches_compatibility_aliases(monkeypatch, toy_model):
    """Supported contexts use one isalgo family; aliases cannot restore an old fallback chain."""
    from vdjtools.model import bundled, junction, native, viterbi

    sources = []
    def load(locus, source, **kw):
        sources.append(source)
        return toy_model

    batches = []
    def infer(model, sequences, **kw):
        batches.append(len(sequences))
        return pl.DataFrame({"cdr3_nt": [None] * len(sequences)}, schema={"cdr3_nt": pl.String})

    monkeypatch.setattr(bundled, "load_bundled", load)
    monkeypatch.setattr(junction, "_reachable", lambda model: model)
    monkeypatch.setattr(viterbi, "infer_nt_batch", infer)
    monkeypatch.setattr(native, "best_aa_scenarios_batch", lambda *a, **kw: pl.DataFrame())
    results = [annotate_junctions(*TRB, species="human", model_source=source, threads=1)
               for source in ("auto", "isalgo", "learned", "arda")]
    assert all(result.equals(results[0]) for result in results[1:])
    assert sources == ["isalgo", "isalgo", "learned", "arda"]
    assert batches == [len(TRB[0])] * 4


def test_a_pinned_set_that_does_not_exist_is_not_silently_ignored():
    """A missing organism/model combination is a hole; an unknown source name is an error."""
    mouse = (["CASSLAPGATNEKLFF"], ["TRBV13-1*01"], ["TRBJ1-4*01"])
    assert annotate_junctions(*mouse, species="mouse")["cdr3_nt"][0] is not None
    assert annotate_junctions(*mouse, species="mouse", model_source="olga")["cdr3_nt"][0] is None
    with pytest.raises(ValueError, match="model_source"):
        annotate_junctions(*mouse, species="mouse", model_source="nosuchset")


def test_a_blank_call_is_proposed_and_says_so_in_the_frame():
    """`proposed` has to REACH the caller: it is what lets a consumer delete its segment proposer.

    A submission may leave one side out -- 3,130 of VDJdb's 192,726 curation keys do -- and arda
    proposes that side from the junction, the locus coming from the side that IS named. Saying which
    side was never curated is a different fact from the `allele` flag, which means the submission
    named a different allele of the same gene.
    """
    out = annotate_junctions(["CASSLAPGATNEKLFF"] * 3, ["", "TRBV9*01", "TRBV9*01"],
                             ["TRBJ2-1*01", "", "TRBJ1-4*01"], species="human")
    assert "proposed" in JUNCTION_COLUMNS and "proposed" in out.columns
    assert out["proposed"].to_list() == ["V", "J", ""]
    # and a proposed side is a complete row, not a stub
    for r in out.head(2).iter_rows(named=True):
        assert r["good"] and r["cdr3_nt"] and r["v_end_nt"] >= 0 and r["j_start_nt"] >= 0


def test_an_organism_with_no_fitted_model_still_answers():
    """Rhesus has no bundled fit in any set, and used to come back with nothing at all.

    A FITTED model exists for human (seven loci) and mouse (TRA/TRB); every other organism arda
    ships germline for got no nucleotides, no D and no bounds. On VDJdb's own corpus that was 1,457
    rhesus keys answered ZERO times. ``"auto"``'s last rung is a germline scaffold built by
    ``from_arda(locus, organism=)`` -- the templated flanks come from that organism's own germline
    and only the N region is a default, which is exactly the part a fit would have improved. That
    lifts rhesus TRB to 1,379 of 1,383 keys and rhesus TRA to 73 of 74.
    """
    got = annotate_junctions(["CASSLAPGATNEKLFF", "CAGGQEQFF"],
                             ["TRBV7-9*01", "TRBV9*01"], ["TRBJ1-4*01", "TRBJ2-1*01"],
                             species="MacacaMulatta")
    assert got.height == 2
    assert got["cdr3_nt"].null_count() == 0, "a scaffold organism came back with no nucleotides"
    assert got["d_call"].null_count() == 0, "a scaffold organism came back with no D"
    assert (got["locus"] == "TRB").all()


def test_a_junction_naming_neither_side_gets_a_locus_proposed():
    """Both calls blank is a curation gap, not a refusal -- 522 of VDJdb's `chunks` records.

    arda proposes the locus as well as the two sides (every locus competes, both ends must agree),
    so the row reaches this pipeline like any other and comes back with nucleotides and a D.
    """
    got = annotate_junctions(["CASSSANYGYTF", "CAAHYGNKLVF"], ["", ""], ["", ""], species="human")
    assert list(got["locus"]) == ["TRB", "TRA"]
    assert got["cdr3_nt"].null_count() == 0
    assert all("V" in p and "J" in p for p in got["proposed"]), got["proposed"].to_list()
    assert got["d_call"][0] is not None and got["d_call"][1] is None   # TRB has a D, TRA has none


def test_d_posterior_reduction_has_one_value_across_repeated_hash_grouping():
    """Consumers may initialise multithreaded Polars before importing vdjtools."""
    import os
    import subprocess
    import sys

    probe = r"""
import polars as pl
assert pl.thread_pool_size() == 4
from vdjtools.model.junction import _d_posteriors
import sys
scenarios = pl.read_parquet(sys.argv[1])
values = {_d_posteriors(scenarios).filter(pl.col("row") == 1228)["_p"][0]
          for _ in range(128)}
assert values == {0.8072042852066905}, values
"""
    from pathlib import Path

    fixture = Path(__file__).parents[1] / "data" / "junction_d_scenarios.parquet"
    result = subprocess.run([sys.executable, "-c", probe, str(fixture)], timeout=30,
                            capture_output=True, text=True,
                            env={**os.environ, "POLARS_MAX_THREADS": "4"})
    assert result.returncode == 0, result.stdout + result.stderr
