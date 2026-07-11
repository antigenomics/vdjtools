"""Tests for Parquet clonotype ingestion (vdjtools.io.read_parquet + auto-detect).

Parquet is the at-scale input format: these pin the canonical round-trip (a frame
written by polars reads back byte-identical), AIRR-name mapping, extension-based
format sniffing, and the count/CDR3 edge cases.
"""
import polars as pl
import pytest

from vdjtools import io as vio
from vdjtools.io import schema as S


def _canon():
    df = pl.DataFrame({
        S.V_CALL: ["TRBV5-1", "TRBV7-9"], S.J_CALL: ["TRBJ2-1", "TRBJ2-7"],
        S.CDR3_AA: ["CASSL", "CASST"], S.CDR3_NT: ["TGTGCC", "TGTAGC"],
        S.COUNT: [10, 5],
    })
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_parquet_round_trip_is_identical(tmp_path):
    """A canonical frame written with write_parquet reads back byte-identical —
    including cdr3_nt, which is NOT an AIRR alias and must be preserved, not dropped."""
    canon = _canon()
    p = tmp_path / "s.parquet"
    canon.write_parquet(p)
    back = vio.read_parquet(p)
    assert back.equals(canon)
    assert back.columns == S.COLUMNS + [S.LOCUS]


def test_sniff_and_read_detect_parquet_by_extension(tmp_path):
    canon = _canon()
    for suffix in (".parquet", ".pq"):
        p = tmp_path / f"s{suffix}"
        canon.write_parquet(p)
        assert vio.sniff_format(p) == "parquet"
        assert vio.read(p).equals(canon)          # auto-detect, no header open


def test_parquet_maps_airr_names(tmp_path):
    """A parquet whose columns use AIRR source names maps to canonical, preferring
    junction_aa (anchors included) and duplicate_count."""
    p = tmp_path / "airr.parquet"
    pl.DataFrame({
        "v_call": ["TRBV5-1"], "j_call": ["TRBJ2-1"],
        "junction_aa": ["CASSL"], "junction": ["TGTGCC"],
        "cdr3_aa": ["ASS"],                        # IMGT (anchors excluded) — ignored
        "duplicate_count": [7],
    }).write_parquet(p)
    df = vio.read_parquet(p)
    assert df[S.CDR3_AA][0] == "CASSL"             # junction preferred over cdr3_aa
    assert df[S.COUNT][0] == 7


def test_parquet_first_call_and_default_count(tmp_path):
    """Ambiguous calls reduce to the first; a missing count defaults to 1."""
    p = tmp_path / "x.parquet"
    pl.DataFrame({
        "v_call": ["TRBV12-4, TRBV12-3"], "j_call": ["TRBJ2-1"],
        "junction_aa": ["CASSL"],
    }).write_parquet(p)
    df = vio.read_parquet(p)
    assert df[S.V_CALL][0] == "TRBV12-4"
    assert df[S.COUNT][0] == 1
    assert df[S.FREQ][0] == 1.0


def test_parquet_without_cdr3_aa_raises(tmp_path):
    p = tmp_path / "bad.parquet"
    pl.DataFrame({"v_call": ["TRBV5-1"], "duplicate_count": [3]}).write_parquet(p)
    with pytest.raises(ValueError, match="CDR3 aa"):
        vio.read_parquet(p)
