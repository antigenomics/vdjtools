"""Tests for vdjtools.io — schema coercion, native/AIRR readers, auto-detect, batch."""
import gzip

import polars as pl
import pytest

from vdjtools import io as vio
from vdjtools.io import schema as S

NATIVE_HEADER = "count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\tVEnd\tDStart\tDEnd\tJStart\n"


def _write_native(path, rows, gzipped=True):
    """Write a tiny native vdjtools table; rows are (count, cdr3nt, cdr3aa, v, d, j)."""
    opener = gzip.open if gzipped else open
    with opener(path, "wt") as f:
        f.write(NATIVE_HEADER)
        for count, nt, aa, v, d, j in rows:
            f.write(f"{count}\t0\t{nt}\t{aa}\t{v}\t{d}\t{j}\t3\t-1\t-1\t5\n")


def test_schema_normalize_adds_missing_and_types():
    df = pl.DataFrame({S.V_CALL: ["TRBV1"], S.JUNCTION_AA: ["CASSL"], S.COUNT: ["7"]})
    out = S.normalize(df, recompute_freq=True)
    assert set(S.COLUMNS) <= set(out.columns)
    assert out[S.COUNT].dtype == pl.Int64
    assert out[S.FREQ].dtype == pl.Float64
    assert out[S.C_CALL].to_list() == [None]        # missing -> null
    assert out[S.FREQ].to_list() == [1.0]           # single clonotype -> freq 1


def test_locus_helpers():
    assert S.locus_of("TRBV12-3*01") == "TRB"
    assert S.locus_of(None) is None
    df = S.add_locus(pl.DataFrame({S.V_CALL: ["TRBV1", "IGHV3", None]}))
    assert df[S.LOCUS].to_list() == ["TRB", "IGH", None]


def test_read_vdjtools_native(tmp_path):
    p = tmp_path / "s.tsv.gz"
    _write_native(p, [
        (10, "TGTGCC", "CASTV", "TRBV12-4, TRBV12-3", ".", "TRBJ1-1"),
        (30, "TGTGCA", "CASSF", "TRBV20-1", "TRBD1", "TRBJ2-1"),
    ])
    df = vio.read_vdjtools(p)
    assert df.height == 2
    # ambiguous V reduced to the first call; locus derived; freq recomputed
    assert df[S.V_CALL].to_list() == ["TRBV12-4", "TRBV20-1"]
    assert df[S.LOCUS].to_list() == ["TRB", "TRB"]
    assert df[S.D_CALL].to_list() == [None, "TRBD1"]      # '.' -> null
    assert df[S.C_CALL].to_list() == [None, None]         # no C in native
    assert df[S.FREQ].to_list() == [0.25, 0.75]
    assert df[S.COUNT].dtype == pl.Int64


def test_read_vdjtools_commented_header(tmp_path):
    # Some legacy exports comment out the header line (``#count freq ...``); the
    # reader must still map the first column to ``count``.
    p = tmp_path / "h.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write("#" + NATIVE_HEADER)
        f.write("5\t0\tTGTGCC\tCASTV\tTRBV7-2\tTRBD1\tTRBJ1-1\t3\t-1\t-1\t5\n")
    df = vio.read_vdjtools(p)
    assert df.height == 1
    assert df[S.COUNT].to_list() == [5]
    assert df[S.JUNCTION_AA].to_list() == ["CASTV"]
    assert df[S.LOCUS].to_list() == ["TRB"]


def test_read_airr_collapses_per_read(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text(
        "v_call\tj_call\tjunction_aa\tc_call\n"
        "TRBV12-3\tTRBJ1-1\tCCASSLF\tIGHM\n"
        "TRBV12-3\tTRBJ1-1\tCCASSLF\tIGHM\n"
        "TRBV20-1\tTRBJ2-1\tCCASSFF\t\n"
    )
    df = vio.read_airr(p)
    assert df.height == 2                                  # two identical rows collapsed
    row = df.filter(pl.col(S.JUNCTION_AA) == "CCASSLF").row(0, named=True)
    assert row[S.COUNT] == 2
    assert abs(row[S.FREQ] - 2 / 3) < 1e-9
    assert df.filter(pl.col(S.JUNCTION_AA) == "CCASSFF")[S.C_CALL].to_list() == [None]


def test_read_airr_prefers_junction_over_imgt_cdr3(tmp_path):
    # Both junction_aa (anchors included) and IMGT cdr3_aa (anchors excluded, 2 aa
    # shorter) are present; canonical cdr3_aa must be the JUNCTION.
    p = tmp_path / "j.tsv"
    p.write_text(
        "v_call\tj_call\tcdr3_aa\tjunction_aa\tcdr3\tjunction\n"
        "TRBV12-3\tTRBJ1-1\tASSLR\tCASSLRF\tGCTAGT\tTGTGCTAGTTTT\n"
    )
    df = vio.read_airr(p)
    assert df[S.JUNCTION_AA].to_list() == ["CASSLRF"]          # junction, not IMGT ASSLR
    assert df[S.JUNCTION_AA].str.len_chars().to_list() == [7]  # 2 longer than ASSLR
    assert df[S.JUNCTION_NT].to_list() == ["TGTGCTAGTTTT"]     # junction nt, not cdr3 nt


def test_read_airr_collapse_key_ignores_d_call(tmp_path):
    # Same (v, j, cdr3nt) with d_call TRBD1 in one row and null in another must
    # collapse to ONE clonotype (D is not part of clonotype identity).
    p = tmp_path / "d.tsv"
    p.write_text(
        "v_call\td_call\tj_call\tjunction\tjunction_aa\n"
        "TRBV12-3\tTRBD1\tTRBJ1-1\tTGTGCTAGT\tCASS\n"
        "TRBV12-3\t\tTRBJ1-1\tTGTGCTAGT\tCASS\n"
    )
    df = vio.read_airr(p)
    assert df.height == 1
    assert df[S.COUNT].to_list() == [2]                    # counts summed
    assert df[S.D_CALL].to_list() == ["TRBD1"]             # representative non-null D


def test_read_vdjtools_rejects_foreign_header(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("foo\tbar\n1\t2\n")
    with pytest.raises(ValueError):
        vio.read_vdjtools(p)
