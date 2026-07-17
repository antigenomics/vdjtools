"""Conformance tests for the legacy-format converters (vdjtools.io.convert).

Oracles are worked examples extracted from the legacy Groovy parsers + their sample fixtures
(tests/python/fixtures/legacy/*.txt.gz, carried over from the legacy-1.x branch).
"""
from pathlib import Path

import polars as pl
import pytest

from vdjtools.io import convert
from vdjtools.io.schema import COLUMNS, LOCUS

FIX = Path(__file__).parent / "fixtures" / "legacy"

# format -> (reader, fixture, [expected clonotypes present in the output]).
# Each expected row is a verified worked example (input row -> canonical values).
CASES = {
    "mixcr": (convert.read_mixcr, "mixcr.txt.gz", [
        dict(v="TRBV13", j="TRBJ2-4", c="TRBC1", aa="CASSLGENIQYF",
             nt="TGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC", count=16988),
    ]),
    "mixcr3": (convert.read_mixcr, "mixcr.3.txt.gz", None),
    "migec": (convert.read_migec, "migec.txt.gz", None),
    # MiTCR / tcR dotted dialect. `D.gene` here is the ambiguous "TRBD1, TRBD2" -> first call.
    "mitcr": (convert.read_mitcr, "mitcr.txt.gz", [
        dict(v="TRBV13", j="TRBJ2-4", d="TRBD1", aa="CASSLGENIQYF",
             nt="TGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC", count=4354),
        dict(v="TRBV29-1", j="TRBJ2-1", d="TRBD1", aa="CSVEIWDSSYNEQFF",
             nt="TGCAGCGTTGAAATTTGGGATAGCTCCTACAATGAGCAGTTCTTC", count=2312),
        dict(v="TRBV7-6", j="TRBJ1-4", d="TRBD2", aa="CASSLAPGATNEKLFF", count=1435),
    ]),
    "rtcr": (convert.read_rtcr, "rtcr.txt.gz", [
        dict(v="TRBV19", j="TRBJ2-7", aa="CARMGQLSYEQYF",
             nt="TGTGCCAGGATGGGACAACTTTCCTACGAGCAGTACTTC", count=11),
        dict(v="TRBV7-9", j="TRBJ2-5", aa="CATLAQDQETQYF",
             nt="TGTGCCACGTTAGCTCAGGATCAAGAGACCCAGTACTTC", count=7),
    ]),
    "imgt": (convert.read_imgt, "imgthighvquest.txt.gz", [
        dict(v="IGHV2-26", d="IGHD1-20", j="IGHJ1", aa="CA**L_TFQHW",
             nt="TGTGCATAATAACTGGAACCTTCCAGCACTGG"),
    ]),
    "immunoseqv2": (convert.read_immunoseq, "immunoseqv2.txt.gz", [
        dict(v="TRBV29-1", j="TRBJ2-6", aa="CSVEDSGANVLTF",
             nt="TGCAGCGTCGAAGACTCTGGGGCCAACGTCCTGACTTTC"),
    ]),
    "immunoseq": (convert.read_immunoseq, "immunoseq.txt.gz", [
        dict(v="TRBV29-1", j="TRBJ2-6", aa="CSVEDSGANVLTF"),  # same underlying data as v2
    ]),
    # NB: two Vidjil clones share this junction, so the collapsed count (201582) exceeds the
    # per-clone reads (189991) — presence + nt is the oracle here, not the pre-collapse count.
    "vidjil": (convert.read_vidjil, "vidjil.txt.gz", [
        dict(v="IGHV3-9", j="IGHJ6", aa="CAPGGMDVW", nt="TGTGCACCCGGAGGTATGGACGTCTGG"),
    ]),
    # TRUST4 report: CDR3 = junction (anchors included); C is the constant/isotype hit.
    "trust4": (convert.read_trust4, "trust4.txt.gz", [
        dict(v="TRBV13", j="TRBJ2-4", c="TRBC2", aa="CASSLGENIQYF",
             nt="TGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC", count=152),
        dict(v="TRBV19", j="TRBJ2-7", c="TRBC1", aa="CARMGQLSYEQYF", count=88),
    ]),
}


_SNIFF = {
    "mixcr.txt.gz": "mixcr", "mixcr.3.txt.gz": "mixcr", "migec.txt.gz": "migec",
    "mitcr.txt.gz": "mitcr",
    "rtcr.txt.gz": "rtcr", "imgthighvquest.txt.gz": "imgt", "vidjil.txt.gz": "vidjil",
    "immunoseq.txt.gz": "immunoseq", "immunoseqv2.txt.gz": "immunoseq",
    "trust4.txt.gz": "trust4",
}


@pytest.mark.parametrize("fixture,fmt", list(_SNIFF.items()))
def test_sniff_and_autoread(fixture, fmt):
    """io.read(fmt='auto') detects each legacy format and returns the canonical frame."""
    from vdjtools import io as vio
    assert vio.sniff_format(FIX / fixture) == fmt
    df = vio.read(FIX / fixture)  # auto-detect + dispatch
    assert df.columns == [*COLUMNS, LOCUS]
    assert df.height > 0


@pytest.mark.parametrize("name", list(CASES))
def test_converter_conformance(name):
    reader, fixture, oracles = CASES[name]
    df = reader(FIX / fixture)

    # structural: canonical schema, non-empty, no null keys, clean nt.
    assert df.columns == [*COLUMNS, LOCUS], df.columns
    assert df.height > 0
    assert df[["v_call", "j_call", "junction_aa", "junction_nt"]].null_count().sum_horizontal()[0] == 0
    assert df["junction_nt"].str.contains(r"^[ACGT]+$").all()
    assert (df["duplicate_count"] > 0).all()
    assert abs(df["frequency"].sum() - 1.0) < 1e-9

    for want in oracles or []:
        hit = df.filter(pl.col("junction_aa") == want["aa"])
        assert hit.height >= 1, f"{name}: {want['aa']} not found; sample={df.head(3).to_dicts()}"
        row = hit.row(0, named=True)
        assert row["v_call"] == want["v"], f"{name} v_call: {row['v_call']} != {want['v']}"
        assert row["j_call"] == want["j"], f"{name} j_call: {row['j_call']} != {want['j']}"
        if "d" in want:
            assert row["d_call"] == want["d"], f"{name} d_call: {row['d_call']} != {want['d']}"
        if "c" in want:
            assert row["c_call"] == want["c"], f"{name} c_call: {row['c_call']} != {want['c']}"
        if "nt" in want:
            assert row["junction_nt"] == want["nt"], f"{name} junction_nt mismatch"
        if "count" in want:
            assert row["duplicate_count"] == want["count"], \
                f"{name} count: {row['duplicate_count']} != {want['count']}"


def test_reader_rejects_wrong_columns(tmp_path):
    """Each TSV reader fails loudly on a table whose columns don't match its format."""
    p = tmp_path / "wrong.tsv"
    p.write_text("colA\tcolB\n1\t2\n")
    for reader in (convert.read_mixcr, convert.read_migec, convert.read_mitcr,
                   convert.read_rtcr, convert.read_imgt, convert.read_immunoseq):
        with pytest.raises(ValueError, match="not a|not an"):
            reader(p)


def test_mitcr_not_confused_with_migec(tmp_path):
    """MiGEC and MiTCR both carry a CDR3-nt column; the dialects differ only by dots vs spaces,
    so each sniff must reject the other's header rather than mis-dispatch."""
    from vdjtools import io as vio
    migec_hdr = ("Count\tCDR3 nucleotide sequence\tCDR3 amino acid sequence\t"
                 "V segments\tJ segments\tD segments\n1\tTGT\tC\tTRBV9\tTRBJ2-3\tTRBD1\n")
    p = tmp_path / "m.tsv"
    p.write_text(migec_hdr)
    assert vio.sniff_format(p) == "migec"
    with pytest.raises(ValueError, match="not a MiTCR"):
        convert.read_mitcr(p)


def test_read_vidjil_skips_segless_clones(tmp_path):
    """A Vidjil clone with no ``seg`` (or no ``seg.junction``) is skipped."""
    import json
    p = tmp_path / "x.vidjil"
    p.write_text(json.dumps({"clones": [{"sequence": "ACGT"}, {"seg": {}}]}))
    assert convert.read_vidjil(p).height == 0


def test_convert_helper_edges():
    """Translation / normalisation / count helpers on their edge inputs."""
    assert convert._finalize([]).columns == [*COLUMNS, LOCUS]
    assert convert._finalize([]).height == 0
    assert convert.translate("") == ""
    assert convert.to_unified_cdr3aa(None) is None
    assert convert._to_int("null", "x") == 0          # nothing numeric → 0
    assert convert._to_int("0", "50") == 50           # skip non-positive → fallback


def test_mixcr_reads_c_gene_and_bcr_isotype():
    """MiXcr C-hit column → ``c_call``: TCR constant (v1/2) and BCR isotype (v3/4)."""
    tcr = convert.read_mixcr(FIX / "mixcr.txt.gz")
    assert tcr["c_call"].drop_nulls().is_in(["TRBC1", "TRBC2"]).all()
    bcr = convert.read_mixcr(FIX / "mixcr.3.txt.gz")            # allCHitsWithScore = isotype
    assert bcr.filter(pl.col("junction_aa") == "CARKKLEGRFDYW")["c_call"][0] == "IGHG4"
    assert bcr["c_call"].drop_nulls().str.starts_with("IGH").all()


def test_trust4_skips_partial_and_missing_segment_clones():
    """TRUST4 ``partial`` / non-ACGT CDR3s and gene-less (``*``) clones are dropped."""
    df = convert.read_trust4(FIX / "trust4.txt.gz")
    assert df.height == 2                                        # 2 of 4 rows survive
    assert set(df["c_call"]) == {"TRBC1", "TRBC2"}
    assert df["junction_nt"].str.contains(r"^[ACGT]+$").all()


def test_read_arda_maps_airr_and_nulls_empty_calls(tmp_path):
    """arda's AIRR output reads via read_arda; its literal ``\"\"`` empty gene call → null."""
    p = tmp_path / "x.airr.tsv"
    p.write_text(
        "sequence_id\tlocus\tv_call\td_call\td2_call\tj_call\tc_call\tjunction\tjunction_aa\tduplicate_count\n"
        'r1\tIGK\tIGKV1-33*01\t""\t\tIGKJ2*02\tIGKC\tTGTCAACAGTATGACAATCTCCCGTACACTTTT\tCQQYDNLPYTF\t40\n'
    )
    from vdjtools import io as vio
    assert vio.sniff_format(p) == "arda"                        # d2_call routes to read_arda
    df = vio.read(p)                                            # auto-detect + dispatch
    assert df.columns == [*COLUMNS, LOCUS]
    r = df.row(0, named=True)
    assert r["d_call"] is None                                  # arda's ``""`` empty-D nulled
    assert (r["c_call"], r["junction_aa"], r["duplicate_count"], r["locus"]) == \
        ("IGKC", "CQQYDNLPYTF", 40, "IGK")


def test_immunoseq_count_falls_back_when_templates_zero(tmp_path):
    """immunoSEQ v1 count uses ``reads`` when ``templates`` is 'null' OR '0' (regression)."""
    hdr = ("rearrangement\tamino_acid\tframe_type\tv_index\tcdr3_length\t"
           "v_gene\tj_gene\ttemplates\treads\n")
    nt = "TGTGCCAGCAGCTTAGGGGAAAACATTCAGTACTTC"
    for templ in ("null", "0"):
        p = tmp_path / f"is_{templ}.tsv"
        p.write_text(hdr + f"{nt}\tCASSLGENIQYF\tIn\t0\t15\tTCRBV13-01\tTCRBJ02-04\t{templ}\t50\n")
        df = convert.read_immunoseq(p)
        assert df.height == 1, f"templates={templ!r} dropped the clonotype"
        assert df["duplicate_count"][0] == 50
