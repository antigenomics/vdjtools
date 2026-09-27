"""A V CDR3-region germline starts at Cys104 and a J one ends at Phe/Trp118 -- for every allele.

``TRBV4-3*02`` did neither, and nothing said so. OLGA's ``human_T_beta`` records its anchor as 267
against its own 287-nt germline while ``*01`` is 284 nt at the same 267, so OLGA's cut segment for
``*02`` begins one Leu codon (``CTC``) before Cys104. OLGA never returns 0 for a *gene*-level TRBV4-3
query -- it marginalises over the gene's alleles and ``*01`` carries the mass -- but
``collapse_alleles`` ranked the representative by germline length first, ``*02``'s broken cut is
24 nt against ``*01``'s correct 21, so the collapsed gene inherited the broken one relabelled ``*01``
and ``pgen_aa`` was **exactly 0** for every real TRBV4-3 junction, with no error raised.

Two independent defences, because either alone leaves the other's failure reachable: the germline is
re-anchored on load (so an allele-level query is right too), and the representative ranking refuses an
out-of-frame candidate (so a germline nobody can repair is never promoted).
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import collapse_alleles, load_bundled, native, translate
from vdjtools.model.io import _anchored, repair_anchors

LOCI = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")

#: J alleles that do not end at Phe/Trp and must be left alone rather than re-anchored onto a
#: coincidental codon. Two different reasons:
#:
#: * ``TRBJ2-7*02`` / ``IGHJ6*02`` -- ORF alleles that genuinely template another residue
#:   (``SYEQYV`` rather than the conserved ``SYEQYF``). Nothing to fix.
#: * ``TRAJ35*01`` -- **unresolved, and not resolvable from our reference**: arda's recorded anchor
#:   (25) points at a Cys codon, so the germline reads ``CGSGTQVIVLP`` and neither the stored slice
#:   nor arda's own says where Phe118 is. Settling it needs IMGT arbitration, not a search, and a
#:   search is exactly what re-anchors onto a coincidence. It is the survivor of the 11-allele family
#:   in ``from_olga(derive_orf=True)``; the repair fixes the other 9 here and 7 of those now match
#:   arda exactly. Sizing unchanged: 0 of VDJdb's 30,937 human TRA records use any of them.
GENUINE_EXCEPTIONS = {"TRBJ2-7*02", "IGHJ6*02", "TRAJ35*01"}


@pytest.mark.parametrize("source", ["olga", "learned"])
@pytest.mark.parametrize("locus", LOCI)
def test_every_v_germline_starts_at_the_conserved_cysteine(source, locus):
    """Out of frame after the repair means no offset works -- a missing Cys, not a shifted index."""
    g = load_bundled(locus, source, collapse=False).genomic["genes_v"]
    bad = [a for a, c in zip(g["v_allele"], g["cut_segment"]) if c and translate(c[:3]) != "C"]
    # A pseudogene with no Cys at all is a hole, not a shifted anchor; what must never appear is an
    # allele whose germline simply starts a codon early.
    for a in bad:
        cut = g.filter(pl.col("v_allele") == a)["cut_segment"][0]
        full = g.filter(pl.col("v_allele") == a)["full_germline"][0]
        assert not any(translate(full[i:i + 3]) == "C" and full[i:] == cut[len(cut) - len(full[i:]):]
                       for i in range(0, max(len(full) - 3, 0), 3)), f"{a} is repairable and was not"


@pytest.mark.parametrize("source", ["olga", "learned"])
@pytest.mark.parametrize("locus", LOCI)
def test_every_j_germline_ends_at_the_conserved_phe_or_trp(source, locus):
    g = load_bundled(locus, source, collapse=False).genomic["genes_j"]
    bad = {a for a, c in zip(g["j_allele"], g["cut_segment"])
           if c and translate(c[-3:]) not in ("F", "W")}
    assert bad <= GENUINE_EXCEPTIONS | {a for a in bad if a.endswith("P*01") or "/OR" in a}, bad


def test_trbv4_3_is_anchored_and_scores_real_junctions():
    """The filed case, at allele resolution and through the collapsed gene."""
    for source in ("olga", "learned"):
        raw = load_bundled("TRB", source, collapse=False).genomic["genes_v"]
        for allele in ("TRBV4-3*01", "TRBV4-3*02"):
            cut = raw.filter(pl.col("v_allele") == allele)["cut_segment"][0]
            assert translate(cut[:3]) == "C", f"{source} {allele} is not Cys-anchored: {cut[:9]}"
        m = load_bundled("TRB", source)
        p = native.pgen_aa(m, "CASSQDLNTEAFF", "TRBV4-3*01", "TRBJ1-1*01")
        assert p > 0, f"{source}: a real TRBV4-3 junction still scores 0"


def test_the_repair_never_shortens_a_germline_it_cannot_confirm():
    """A J search finds an earlier in-frame Phe by coincidence; mouse ``TRAJ19*01`` measured that.

    Its arda germline is 30 nt and does not end at Phe/Trp at all, and a two-codon search cut it to
    24. So J is repaired only by undoing the known wrong-side slice, never by searching.
    """
    g = load_bundled("TRA", "arda", organism="mouse", collapse=False).genomic["genes_j"]
    row = g.filter(pl.col("j_allele") == "TRAJ19*01")
    assert row.height == 1 and len(row["cdr3_segment"][0]) == 30


def test_the_representative_refuses_an_out_of_frame_germline_even_when_it_is_longer():
    """The second defence, exercised on a frame built to fail: longer must not beat anchored."""
    genomic = {"genes_v": pl.DataFrame({
        "v_allele": ["TESTV1*01", "TESTV1*02"],
        "gene": ["TESTV1", "TESTV1"],
        "full_germline": ["", ""],
        "cdr3_segment": ["TGCGCCAGCAGC", "CTCTGCGCCAGCAGC"],
        "cut_segment": ["TGCGCCAGCAGC", "CTCTGCGCCAGCAGC"],   # *02 is longer AND out of frame
        "anchor": [-1, -1],
        "functional": [True, True],
    })}
    from vdjtools.model.collapse import _representative
    cut = dict(zip(genomic["genes_v"]["v_allele"], genomic["genes_v"]["cut_segment"]))
    assert _representative(list(cut), cut, {}, {}, "v") == "TESTV1*01"
    # ... and with no anchored candidate the gate is inert, so length still leads (TRBV23/OR9-2).
    empty = {"X*01": "", "X*02": "CTCTGCGCC"}
    assert _representative(list(empty), empty, {}, {}, "v") == "X*02"


def test_repair_anchors_reports_what_it_moved():
    """A repair that edits a germline silently is the shape of the bug it fixes."""
    frame = {"genes_v": pl.DataFrame({
        "v_allele": ["TESTV1*01"],
        "gene": ["TESTV1"],
        "full_germline": ["AAA" * 10 + "CTC" + "TGCGCCAGCAGCCAAGA"],
        "cdr3_segment": ["CTCTGCGCCAGCAGCCAAGA"],
        "cut_segment": ["CTCTGCGCCAGCAGCCAAGA"],
        "anchor": [30],
        "functional": [True],
    })}
    fixed, rep = repair_anchors(frame, report=True)
    assert rep == ["TESTV1*01: anchor 30 -> 33, 20 -> 17 nt"]
    assert _anchored("v", fixed["genes_v"]["cut_segment"][0])


def test_collapsing_a_bundled_model_leaves_no_out_of_frame_v():
    """The end-to-end property: whatever the artifact holds, the model you get is anchored."""
    for locus in LOCI:
        g = collapse_alleles(load_bundled(locus, "olga", collapse=False)).genomic["genes_v"]
        bad = [a for a, c in zip(g["v_allele"], g["cut_segment"]) if c and translate(c[:3]) != "C"]
        assert bad == [], f"{locus}: {bad}"
