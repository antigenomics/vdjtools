"""Float-formatted counts must not become singletons.

pandas writes any integer column that ever held a NaN as float ("5000.0"), and the AIRR/vdjtools
readers parse the TSV as all-Utf8. A naive Utf8->Int64 cast yields null for "5000.0", and the
fill_null(1) / recompute_frequency that follow then silently invert the clonal hierarchy.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools import io as vio


def _write(tmp_path, header, rows):
    p = tmp_path / "s.tsv"
    p.write_text(header + "\n" + "\n".join(rows) + "\n")
    return p


def test_airr_float_count_is_preserved(tmp_path):
    p = _write(tmp_path, "sequence_id\tv_call\tj_call\tjunction_aa\tduplicate_count",
               ["A\tTRBV19*01\tTRBJ2-7*01\tCASSIRSSYEQYF\t5000.0",
                "B\tTRBV20-1*01\tTRBJ2-1*01\tCASSLGETQYF\t3.0"])
    df = vio.read(p, fmt="airr").sort("junction_aa")
    counts = dict(zip(df["junction_aa"], df["duplicate_count"]))
    assert counts["CASSIRSSYEQYF"] == 5000, "a 5000-read clone became a singleton"
    assert counts["CASSLGETQYF"] == 3
    # frequency must reflect the real hierarchy, not 0.5/0.5
    freqs = dict(zip(df["junction_aa"], df["frequency"]))
    assert freqs["CASSIRSSYEQYF"] == pytest.approx(5000 / 5003)


def test_vdjtools_float_count_is_preserved(tmp_path):
    p = _write(tmp_path, "count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj",
               ["5000.0\t0\tTGT\tCASSIRSSYEQYF\tTRBV19\t.\tTRBJ2-7",
                "3.0\t0\tTGT\tCASSLGETQYF\tTRBV20-1\t.\tTRBJ2-1"])
    df = vio.read(p, fmt="vdjtools").sort("junction_aa")
    counts = dict(zip(df["junction_aa"], df["duplicate_count"]))
    assert counts["CASSIRSSYEQYF"] == 5000
    assert counts["CASSLGETQYF"] == 3


def test_integer_count_still_works(tmp_path):
    p = _write(tmp_path, "sequence_id\tv_call\tj_call\tjunction_aa\tduplicate_count",
               ["A\tTRBV19*01\tTRBJ2-7*01\tCASSIRSSYEQYF\t42"])
    df = vio.read(p, fmt="airr")
    assert df["duplicate_count"][0] == 42


def test_unparseable_count_falls_back_to_one(tmp_path):
    p = _write(tmp_path, "sequence_id\tv_call\tj_call\tjunction_aa\tduplicate_count",
               ["A\tTRBV19*01\tTRBJ2-7*01\tCASSIRSSYEQYF\t"])
    df = vio.read(p, fmt="airr")
    assert df["duplicate_count"][0] == 1
