"""Building a model from a custom V(D)J germline library (``from_germline`` and friends).

The point of these tests is that a user's own reference works end to end — scaffold, validate,
generate, score — and that the mistakes which otherwise produce a model that *builds cleanly and
scores wrongly* are caught at build time. The anchor-frame checks matter most: a CDR3 anchor off
by one codon shifts every deletion profile by a constant and is invisible downstream.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import native
from vdjtools.model.generate import generate
from vdjtools.model.io import from_arda, from_germline
from vdjtools.model.reference import read_germline_fasta, validate_germline


def _errors(issues: pl.DataFrame) -> list[str]:
    return issues.filter(pl.col("severity") == "error")["check"].to_list()


def _checks(issues: pl.DataFrame, severity: str) -> list[str]:
    return issues.filter(pl.col("severity") == severity)["check"].to_list()


def test_toy_germline_is_clean(toy_germline):
    assert validate_germline(toy_germline).is_empty()


def test_from_germline_builds_a_vj_model(toy_model):
    assert toy_model.chain_type == "VJ"
    assert toy_model.locus == "TOY"
    assert toy_model.genomic["genes_v"].height == 3
    assert toy_model.genomic["genes_j"].height == 2
    assert "genes_d" not in toy_model.genomic
    assert "vj_ins" in toy_model.tables and "d_gene" not in toy_model.tables
    toy_model.validate()


def test_a_d_allele_makes_it_vdj(toy_model_vdj):
    assert toy_model_vdj.chain_type == "VDJ"
    assert toy_model_vdj.genomic["genes_d"].height == 1
    for ev in ("d_gene", "d_del", "n_d", "vd_ins", "dj_ins", "vd_dinucl", "dj_dinucl"):
        assert ev in toy_model_vdj.tables
    toy_model_vdj.validate()


def test_custom_model_generates_and_scores(toy_model):
    """The whole point: a hand-made reference produces sequences its own Pgen can score."""
    gen = generate(toy_model, 50, seed=0)
    assert gen.height == 50
    assert set(gen["v_call"]) <= set(toy_model.genomic["genes_v"]["v_allele"])
    for seq in gen["junction_nt"].to_list()[:10]:
        assert native.pgen_nt(toy_model, seq) > 0


@pytest.mark.parametrize("drop", ["allele", "segment", "sequence"])
def test_missing_required_column_is_an_error(toy_germline, drop):
    issues = validate_germline(toy_germline.drop(drop))
    assert "germline_columns" in _errors(issues)


def test_no_j_alleles_is_an_error(toy_germline):
    issues = validate_germline(toy_germline.filter(pl.col("segment") != "J"))
    assert "germline_segment_missing" in _errors(issues)


def test_duplicate_allele_is_an_error(toy_germline):
    issues = validate_germline(pl.concat([toy_germline, toy_germline.head(1)]))
    assert "germline_duplicate_allele" in _errors(issues)


def test_bad_segment_letter_is_an_error(toy_germline):
    bad = toy_germline.with_columns(
        segment=pl.when(pl.col("segment") == "V").then(pl.lit("X")).otherwise(pl.col("segment")))
    assert "germline_segment" in _errors(validate_germline(bad))


def test_v_not_starting_on_cys_warns(toy_germline):
    """The commonest custom-library mistake: the anchor is one codon off."""
    shifted = toy_germline.with_columns(
        sequence=pl.when(pl.col("allele") == "TOYV1*01")
        .then(pl.lit("GCCAGCTGT")).otherwise(pl.col("sequence")))
    assert "germline_anchor_frame" in _checks(validate_germline(shifted), "warn")


def test_j_not_ending_on_phe_or_trp_warns(toy_germline):
    shifted = toy_germline.with_columns(
        sequence=pl.when(pl.col("allele") == "TOYJ2*01")
        .then(pl.lit("AACGAGCAGTTTAA")).otherwise(pl.col("sequence")))
    assert "germline_anchor_frame" in _checks(validate_germline(shifted), "warn")


def test_ambiguous_bases_are_dropped_and_reported(toy_germline):
    """IUPAC germline cannot be encoded by the native DP — dropped, but never silently."""
    amb = toy_germline.with_columns(
        sequence=pl.when(pl.col("allele") == "TOYV2*01")
        .then(pl.lit("TGTGCNTCC")).otherwise(pl.col("sequence")))
    assert "germline_ambiguous" in _checks(validate_germline(amb), "warn")
    m = from_germline(amb, locus="TOY", ins_max=3)
    assert "TOYV2*01" not in m.genomic["genes_v"]["v_allele"].to_list()


def test_gene_level_allele_name_warns(toy_germline):
    """The model is allele-keyed; a bare gene name is the documented 2.38x-too-high trap."""
    bare = toy_germline.with_columns(
        allele=pl.when(pl.col("allele") == "TOYV2*01")
        .then(pl.lit("TOYV2")).otherwise(pl.col("allele")))
    assert "germline_gene_level_name" in _checks(validate_germline(bare), "warn")


def test_strict_raises_on_an_error(toy_germline):
    with pytest.raises(ValueError, match="error"):
        from_germline(toy_germline.filter(pl.col("segment") != "J"), locus="TOY")


def test_optional_columns_default(toy_germline):
    """A library with only the three required columns still builds."""
    minimal = toy_germline.select(["allele", "segment", "sequence"])
    m = from_germline(minimal, locus="TOY", ins_max=3)
    assert m.genomic["genes_v"].height == 3
    assert m.genomic["genes_v"]["gene"].to_list()[0] == "TOYV1"


def test_read_germline_fasta(tmp_path):
    v = tmp_path / "v.fasta"
    j = tmp_path / "j.fasta"
    v.write_text(">TOYV1*01\nTGTGCCAGC\n>TOYV2*01\nTGTGCTTCC\n")
    j.write_text(">TRGJ|TOYJ1*01\nAACTATGGCTATACCTTT\n")
    gl = read_germline_fasta(v, j)
    assert set(gl["allele"]) == {"TOYV1*01", "TOYV2*01", "TOYJ1*01"}
    assert set(gl.filter(pl.col("segment") == "V")["allele"]) == {"TOYV1*01", "TOYV2*01"}
    assert validate_germline(gl).is_empty()
    m = from_germline(gl, locus="TOY", ins_max=3)
    assert m.chain_type == "VJ" and m.genomic["genes_v"].height == 2


def test_read_germline_fasta_with_anchors_slices_full_length(tmp_path):
    """With an anchor CSV the FASTAs are full-length and get sliced to the CDR3 region."""
    v = tmp_path / "v.fasta"
    j = tmp_path / "j.fasta"
    anchors = tmp_path / "anchors.csv"
    v.write_text(">TOYV1*01\nAAAAAATGTGCCAGC\n")            # 6 nt of framework, then the Cys codon
    j.write_text(">TOYJ1*01\nAACTATGGCTATACCTTTGGGGGG\n")    # Phe codon at 15, then framework
    anchors.write_text("gene,anchor_index,function\nTOYV1*01,6,F\nTOYJ1*01,15,F\n")
    gl = read_germline_fasta(v, j, anchors=anchors)
    assert gl.filter(pl.col("allele") == "TOYV1*01")["sequence"][0] == "TGTGCCAGC"
    assert gl.filter(pl.col("allele") == "TOYJ1*01")["sequence"][0] == "AACTATGGCTATACCTTT"
    assert validate_germline(gl).is_empty()


def test_from_arda_still_works():
    """The refactor must not change from_arda: same gene set, same chain type, still valid."""
    m = from_arda("TRG")
    assert m.chain_type == "VJ"
    assert m.genomic["genes_v"].height > 5
    assert m.manifest.source == "arda:TRG"
    m.validate()
