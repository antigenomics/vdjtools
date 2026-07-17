"""Parquet round-trip must preserve comma-ambiguous gene calls.

read_airr and scan_cohort keep a tie like "IGHV3-23*01,IGHV3-23D*01" whole; parquet is the
at-scale storage of a canonical frame, so writing and reading it back must not silently drop the
second gene (which model.infer.call_alleles exists to keep).
"""
from __future__ import annotations

import polars as pl

from vdjtools import io as vio


def test_parquet_preserves_ambiguous_v_call(tmp_path):
    src = pl.DataFrame({
        "v_call": ["IGHV3-23*01,IGHV3-23D*01", "IGHV1-2*02"],
        "j_call": ["IGHJ4*02", "IGHJ4*02"],
        "junction_nt": ["TGTACG", "TGTAAA"],
        "junction_aa": ["CASSF", "CASSG"],
        "duplicate_count": [10, 5],
    })
    p = tmp_path / "s.parquet"
    src.write_parquet(p)
    df = vio.read(p, fmt="parquet")
    calls = set(df["v_call"].to_list())
    assert "IGHV3-23*01,IGHV3-23D*01" in calls, "the ambiguous tie was collapsed to one gene"


def test_airr_parquet_roundtrip_is_stable(tmp_path):
    """AIRR -> read_airr -> write_parquet -> read_parquet keeps the same v_call values."""
    tsv = tmp_path / "s.tsv"
    tsv.write_text("sequence_id\tv_call\tj_call\tjunction_aa\tduplicate_count\n"
                   "A\tIGHV3-23*01,IGHV3-23D*01\tIGHJ4*02\tCASSF\t10\n")
    a = vio.read(tsv, fmt="airr")
    p = tmp_path / "s.parquet"
    a.write_parquet(p)
    b = vio.read(p, fmt="parquet")
    assert set(a["v_call"].to_list()) == set(b["v_call"].to_list())
