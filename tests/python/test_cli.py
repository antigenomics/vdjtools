"""Phase 8 — the ``vdjtools`` CLI (typer). Smoke-tests each command via CliRunner.

Self-contained: uses the bundled OLGA models shipped in the wheel, so no OLGA install, network,
or HuggingFace fetch is needed.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest
from typer.testing import CliRunner

from vdjtools.cli import app

runner = CliRunner()


def _airr_sample(df: pl.DataFrame, seed: int) -> pl.DataFrame:
    rng = np.random.default_rng(seed)
    return df.with_columns(
        pl.Series("duplicate_count", rng.integers(1, 20, df.height)),
        pl.col("junction_aa").alias("junction_aa"),
        pl.col("junction_nt").alias("junction"),
    ).select("v_call", "d_call", "j_call", "junction_aa", "junction", "duplicate_count", "productive")


@pytest.fixture(scope="module")
def gen():
    from vdjtools.model import load_bundled
    from vdjtools.model.generate import generate

    return generate(load_bundled("TRB", "olga"), 150, seed=3, productive_only=True)


def test_models_lists_builtins():
    r = runner.invoke(app, ["models"])
    assert r.exit_code == 0, r.stdout
    assert "olga" in r.stdout and "TRB" in r.stdout


def test_generate(tmp_path):
    out = tmp_path / "gen.tsv"
    r = runner.invoke(app, ["generate", "-m", "TRB", "-n", "8", "--seed", "1", "--productive", "-o", str(out)])
    assert r.exit_code == 0, r.stdout
    df = pl.read_csv(out, separator="\t")
    assert df.height == 8
    assert {"junction_aa", "junction_nt", "v_call", "j_call"} <= set(df.columns)


def test_pgen_exact_and_hamming1(tmp_path):
    seqs = tmp_path / "seqs.tsv"
    pl.DataFrame({"cdr3_aa": ["CASSLAPGATNEKLFF", "CAWSVAPDRGGYTF"]}).write_csv(seqs, separator="\t")
    out = tmp_path / "pgen.tsv"
    r = runner.invoke(app, ["pgen", str(seqs), "-m", "TRB", "-o", str(out)])
    assert r.exit_code == 0, r.stdout
    exact = pl.read_csv(out, separator="\t")
    assert "pgen" in exact.columns and (exact["pgen"] > 0).all()

    out1 = tmp_path / "pgen1.tsv"
    r1 = runner.invoke(app, ["pgen", str(seqs), "-m", "TRB", "--mismatches", "1", "-o", str(out1)])
    assert r1.exit_code == 0, r1.stdout
    ball = pl.read_csv(out1, separator="\t")
    # the Hamming-1 ball contains the exact sequence, so its Pgen is strictly larger
    assert (ball["pgen"] >= exact["pgen"]).all() and (ball["pgen"] > exact["pgen"]).any()


def test_pgen_bad_model():
    r = runner.invoke(app, ["pgen", "x.tsv", "-m", "TRX"])
    assert r.exit_code != 0  # unknown locus is a clean error, not a crash


def test_diversity(tmp_path, gen):
    a = tmp_path / "A.tsv"
    _airr_sample(gen[:100], 1).write_csv(a, separator="\t")
    r = runner.invoke(app, ["diversity", str(a)])
    assert r.exit_code == 0, r.stdout
    assert "observed_diversity" in r.stdout and "shannon_wiener" in r.stdout


def test_overlap_and_usage(tmp_path, gen):
    a, b = tmp_path / "A.tsv", tmp_path / "B.tsv"
    _airr_sample(gen[:100], 1).write_csv(a, separator="\t")
    _airr_sample(gen[50:150], 2).write_csv(b, separator="\t")  # 50 shared clonotypes

    ro = runner.invoke(app, ["overlap", str(a), str(b)])
    assert ro.exit_code == 0, ro.stdout
    assert "sample_a\tsample_b" in ro.stdout and "\tD\t" in ro.stdout.replace("  ", "")

    ru = runner.invoke(app, ["segment-usage", str(a), "--segment", "v"])
    assert ru.exit_code == 0, ru.stdout
    assert "v_call" in ru.stdout and "TRB" in ru.stdout

    rs = runner.invoke(app, ["spectratype", str(a)])
    assert rs.exit_code == 0, rs.stdout
    assert "length" in rs.stdout


def test_overlap_needs_two_samples(tmp_path, gen):
    a = tmp_path / "A.tsv"
    _airr_sample(gen[:50], 1).write_csv(a, separator="\t")
    r = runner.invoke(app, ["overlap", str(a)])
    assert r.exit_code != 0  # single sample → clean error
