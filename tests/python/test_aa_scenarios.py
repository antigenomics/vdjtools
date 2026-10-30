"""The amino-acid scenario API: reachable, resolvable, batched, and honest about an empty result.

Filed as antigenomics/vdjtools#179. ``native.best_aa_scenarios`` was the right primitive for "all
plausible ways this junction could have been made" and could not be used as one: it was not exported
at :mod:`vdjtools.model`, it raised on the gene-level V/J calls that every real repertoire carries,
and the per-row Python loop was most of the cost of running it over a clonotype table.
"""
from __future__ import annotations

from vdjtools.cores import available_cores
import polars as pl
import pytest

from vdjtools.model import (
    best_aa_scenarios,
    best_aa_scenarios_batch,
    gene_to_allele,
    load_bundled,
)

#: Real human TRB junctions with the V/J calls as an aligner reports them -- GENE level, no allele.
CASES = [("CASSQDLNTEAFF", "TRBV4-3", "TRBJ1-1"),
         ("CASSLGQGAYEQYF", "TRBV5-1", "TRBJ2-7"),
         ("CASSPRTGELFF", "TRBV7-9", "TRBJ2-2")]


@pytest.fixture(scope="module")
def trb():
    return load_bundled("TRB")


def test_the_scenario_api_is_reachable_from_vdjtools_model():
    """Both entry points must be exported beside ``infer_nt``, not reached through ``model.native``."""
    import vdjtools.model as M

    for name in ("best_aa_scenarios", "best_aa_scenarios_batch", "gene_to_allele"):
        assert name in M.__all__ and hasattr(M, name), name


def test_a_gene_level_call_is_resolved_rather_than_refused(trb):
    """Every V/J call in a real repertoire is gene-level; refusing them all made this unusable."""
    for aa, v, j in CASES:
        got = best_aa_scenarios(trb, aa, v, j, 4)
        assert got, (aa, v, j)
        assert got[0][1].startswith(v + "*") and got[0][3].startswith(j + "*")
        assert all(a[0] >= b[0] for a, b in zip(got, got[1:])), "not ordered by weight"


def test_single_representative_needs_no_allele_selection(trb):
    aa, v, j = CASES[0]
    aliases = gene_to_allele(trb)
    expected = best_aa_scenarios(trb, aa, aliases[v], aliases[j], resolve_genes=False)
    assert expected
    assert best_aa_scenarios(trb, aa, v, j, resolve_genes=False) == expected


def test_refusing_resolution_rejects_a_multiallele_gene(toy_model):
    with pytest.raises(KeyError, match="gene name"):
        best_aa_scenarios(toy_model, "CAF", "TOYV1", "TOYJ1", resolve_genes=False)


def test_an_unknown_call_raises_and_names_it_rather_than_returning_empty(trb):
    """The filed complaint: an empty list with no reason. A naming mismatch must raise, not decline."""
    with pytest.raises(KeyError, match="TRBV999"):
        best_aa_scenarios(trb, "CASSKSQDEQYF", "TRBV999", "TRBJ2-7")


def test_an_empty_result_means_the_dp_explains_nothing_not_a_naming_problem(trb):
    """A residue outside the genetic code has no codons, so there is genuinely no scenario.

    The other cause -- a pinned V whose germline cannot reach the junction -- used to fire on 731 of
    25,000 real TRB clonotypes, 698 of them ``TRBV4-3``, and every one of those was the mis-anchored
    germline of :func:`vdjtools.model.io.repair_anchors`, not a property of the sequence. It is
    deliberately not pinned here: a case that only exists while a germline is broken is a test of the
    breakage.
    """
    assert best_aa_scenarios(trb, "CASSXQDLNTEAFF", "TRBV4-3", "TRBJ1-1") == []
    assert best_aa_scenarios(trb, "CASSQDLNTEAFF", "TRBV4-3", "TRBJ1-1")


def test_the_batch_is_the_per_row_loop_exactly(trb):
    """Same numbers, one row per scenario, and ``row`` indexes the input."""
    aa = [c[0] for c in CASES] * 40
    v = [c[1] for c in CASES] * 40
    j = [c[2] for c in CASES] * 40
    got = best_aa_scenarios_batch(trb, aa, v, j, k=4)
    assert got.columns == ["row", "rank", "w", "v_call", "len_v", "j_call", "len_j",
                           "d_call", "idx5", "idx3", "pos", "d2_call", "idx5_2", "idx3_2", "pos2"]
    for i in (0, 1, 2, 77):
        ref = best_aa_scenarios(trb, aa[i], v[i], j[i], 4)
        rows = got.filter(pl.col("row") == i).sort("rank")
        assert rows.height == len(ref)
        for r, (w, va, lv, ja, lj, da, i5, i3, pos, d2, k5, k3, pos2) in zip(rows.iter_rows(named=True), ref):
            assert r["w"] == w and r["v_call"] == va and r["len_v"] == lv
            assert r["j_call"] == ja and r["len_j"] == lj and r["d_call"] == da
            assert (r["idx5"], r["idx3"], r["pos"]) == (i5, i3, pos)
            assert (r["d2_call"], r["idx5_2"], r["idx3_2"], r["pos2"]) == (d2, k5, k3, pos2)
    assert got.schema['d2_call'] == pl.String


def test_the_batch_result_does_not_depend_on_the_thread_count(trb):
    """The parallelism is disjoint, so it must be identical -- not merely close -- at any width."""
    aa = [c[0] for c in CASES] * 40
    v = [c[1] for c in CASES] * 40
    j = [c[2] for c in CASES] * 40
    one = best_aa_scenarios_batch(trb, aa, v, j, k=4, threads=1)
    many = best_aa_scenarios_batch(trb, aa, v, j, k=4, threads=min(8, available_cores()))
    assert one.equals(many)


def test_a_declined_query_contributes_no_rows_rather_than_a_marginalized_one(trb):
    """An absent ``row`` is how a decline shows up; a silently V/J-agnostic row would be the trap."""
    aa = ["CASSQDLNTEAFF", "CASSXQDLNTEAFF", "CASSPRTGELFF"]
    got = best_aa_scenarios_batch(trb, aa, ["TRBV4-3", "TRBV4-3", "TRBV7-9"],
                                  ["TRBJ1-1", "TRBJ1-1", "TRBJ2-2"], k=2)
    assert set(got["row"].to_list()) == {0, 2}


def test_a_length_mismatch_is_rejected_rather_than_zipped_short(trb):
    with pytest.raises(ValueError, match="same length"):
        best_aa_scenarios_batch(trb, ["CASSPRTGELFF", "CASSPRTGELFF"], ["TRBV7-9"], None)


def test_the_resolver_has_one_definition(trb):
    """``sc.pgen`` had its own copy; #179's point was that every caller reimplements it."""
    import inspect

    from vdjtools.sc import pgen as scp

    assert "gene_to_allele" in inspect.getsource(scp.paired_pgen.__wrapped__
                                                 if hasattr(scp.paired_pgen, "__wrapped__")
                                                 else scp.paired_pgen)
    assert not hasattr(scp, "_gene_to_allele")
    assert gene_to_allele(trb)["TRBV4-3"] == "TRBV4-3*01"
