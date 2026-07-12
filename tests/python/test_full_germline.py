"""arda full-length V/J germline helper (Phase 13 / P1c prerequisite)."""
import pytest

pytest.importorskip("arda")

from vdjtools.model import reference as ref  # noqa: E402


def test_load_full_vj_germline_has_vj_alleles():
    full = ref.load_full_vj_germline("human")
    assert {seg for seg, _ in full} == {"V", "J"}
    assert ("V", "TRBV20-1*01") in full
    assert ("J", "TRBJ2-1*01") in full
    # full V-REGION is far longer than the CDR3-region germline (frameworks present).
    assert len(full[("V", "TRBV20-1*01")]) > 250


def test_arda_full_germline_anchor_consistent_with_cdr3_region():
    fg = ref.arda_full_germline("TRB")
    gl = ref.load_germline("TRB")
    cdr3 = {(r["segment"], r["allele"]): r["sequence"] for r in gl.iter_rows(named=True)}
    assert fg  # non-empty
    for (seg, allele), (full, anchor) in fg.items():
        c = cdr3[(seg, allele)]
        if seg == "V":
            assert full[anchor:] == c          # framework 5' of Cys104 is full[:anchor]
        else:
            assert full[:anchor + 3] == c      # framework 3' of [FW]118 is full[anchor+3:]


def test_arda_germline_stitches_full_contig():
    """full[:anchor_V] + CDR3 + full[anchor_J+3:] assembles a contig with correct flanks."""
    fg = ref.arda_full_germline("TRB")
    fv, av = fg[("V", "TRBV20-1*01")]
    fj, aj = fg[("J", "TRBJ2-1*01")]
    cdr3 = "TGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC"  # a real TRB junction (Cys→Phe)
    contig = fv[:av] + cdr3 + fj[aj + 3:]
    assert cdr3 in contig
    assert contig.startswith(fv[:av]) and contig.endswith(fj[aj + 3:])
    assert len(contig) == av + len(cdr3) + (len(fj) - aj - 3)
    # frameworks flank the junction — no CDR3-region germline duplication at the seam.
    assert contig[:av] == fv[:av]
