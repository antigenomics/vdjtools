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


# --- data commands: convert / downsample / filter / pool ---------------------------------

def test_convert_to_parquet_and_stdout(tmp_path, gen):
    src = tmp_path / "s.tsv"
    _write_airr(src, gen, seed=61)

    pq = tmp_path / "out.parquet"                      # format-aware -o: .parquet -> Parquet
    r = runner.invoke(app, ["convert", str(src), "-o", str(pq)])
    assert r.exit_code == 0, r.stdout
    df = pl.read_parquet(pq)
    assert {"v_call", "j_call", "junction_aa", "duplicate_count", "frequency"} <= set(df.columns)

    r2 = runner.invoke(app, ["convert", str(src)])     # no -o -> canonical TSV on stdout
    assert r2.exit_code == 0, r2.stdout
    assert "junction_aa" in r2.stdout.splitlines()[0]


def test_downsample_reduces_reads(tmp_path, gen):
    src, out = tmp_path / "s.tsv", tmp_path / "ds.tsv"
    _write_airr(src, gen, seed=62)                     # ~31k reads over 150 clonotypes
    r = runner.invoke(app, ["downsample", str(src), "500", "--seed", "1", "-o", str(out)])
    assert r.exit_code == 0, r.stdout
    assert pl.read_csv(out, separator="\t")["duplicate_count"].sum() == 500


def test_filter_coding_and_guard(tmp_path, gen):
    src, out = tmp_path / "s.tsv", tmp_path / "f.tsv"
    _write_airr(src, gen, seed=63)
    r = runner.invoke(app, ["filter", str(src), "--coding", "-o", str(out)])
    assert r.exit_code == 0, r.stdout
    df = pl.read_csv(out, separator="\t")
    assert df.height > 0 and not df["junction_aa"].str.contains(r"\*").any()  # no stop codons
    # --coding and --noncoding together is a clean error, not both-filters-applied
    assert runner.invoke(app, ["filter", str(src), "--coding", "--noncoding"]).exit_code != 0


def test_pool_flat_and_join(tmp_path, gen):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    _write_airr(a, gen[:100], seed=64)
    _write_airr(b, gen[50:150], seed=65)               # 50 clonotypes shared with a
    flat = tmp_path / "pool.tsv"
    r = runner.invoke(app, ["pool", str(a), str(b), "-o", str(flat)])
    assert r.exit_code == 0, r.stdout
    df = pl.read_csv(flat, separator="\t")
    assert df.height <= 200 and "incidence" in df.columns   # shared clonotypes collapse

    joint = tmp_path / "joint.tsv"
    rj = runner.invoke(app, ["pool", str(a), str(b), "--join", "--min-samples", "2", "-o", str(joint)])
    assert rj.exit_code == 0, rj.stdout
    assert pl.read_csv(joint, separator="\t").height >= 1   # the shared ones survive the 2-sample join


def test_pool_needs_two_samples(tmp_path, gen):
    a = tmp_path / "a.tsv"
    _write_airr(a, gen, seed=66)
    assert runner.invoke(app, ["pool", str(a)]).exit_code != 0


# --- Phase 14: dynamics + the two enrichment nulls ---------------------------------------

def _write_airr(path, df, seed, boost=None, scale=1):
    """Write an AIRR sample; `boost` multiplies the first N clonotypes' counts."""
    rng = np.random.default_rng(seed)
    c = rng.integers(20, 400, df.height) * scale
    if boost:
        c[:boost] = c[:boost] * 50
    out = df.with_columns(pl.Series("duplicate_count", c)).select(
        "v_call", "j_call", "junction_aa", pl.col("junction_nt").alias("junction"),
        "duplicate_count")
    out.write_csv(path, separator="\t")


def test_dynamics_classifies_a_pair(tmp_path, gen):
    # --neff is pinned deliberately: this fixture is 150 clonotypes, which cannot support a
    # mean-variance fit (estimate_neff correctly raises on it). CLI tests exercise the CLI;
    # the estimator is tested against planted truths in test_dynamics_paired.py.
    pre, post = tmp_path / "pre.tsv", tmp_path / "post.tsv"
    _write_airr(pre, gen, seed=11)
    _write_airr(post, gen, seed=12, boost=10)
    r = runner.invoke(app, ["dynamics", str(pre), str(post), "--min-total", "4",
                            "--neff", "50000"])
    assert r.exit_code == 0, r.stdout
    head, *rows = r.stdout.strip().splitlines()
    assert "dynamics" in head and "q_value" in head
    assert rows, "no clonotypes returned"
    # the 50x-boosted clonotypes must come back as changed, not as `persistent`
    assert "expanded" in r.stdout


def test_dynamics_umi_skips_the_downscale(tmp_path, gen):
    # --umi is the thesis's p.86 case: counts are already molecule counts, so there is no
    # oversampling to undo and estimating N_eff would downscale twice.
    pre, post = tmp_path / "a.tsv", tmp_path / "b.tsv"
    _write_airr(pre, gen, seed=21)
    _write_airr(post, gen, seed=22)
    r = runner.invoke(app, ["dynamics", str(pre), str(post), "--umi", "--min-total", "4"])
    assert r.exit_code == 0, r.stdout


def test_dynamics_errors_on_an_unfittable_pair(tmp_path, gen):
    # A pair too shallow to fit N_eff must fail loudly with a usable message, not return a
    # guessed number. 150 clonotypes is genuinely unfittable -- and a silent fallback (say,
    # "use the read depth") would be wrong by the whole oversampling factor with no error.
    p = tmp_path / "t.tsv"
    _write_airr(p, gen, seed=51)
    r = runner.invoke(app, ["dynamics", str(p), str(p)])
    assert r.exit_code != 0
    assert "bins" in r.output or "shallow" in r.output


def test_tcrnet_command(tmp_path, gen):
    pytest.importorskip("vdjmatch")
    p = tmp_path / "s.tsv"
    _write_airr(p, gen, seed=31)
    r = runner.invoke(app, ["tcrnet", str(p), "--locus", "TRB"])
    assert r.exit_code == 0, r.stdout
    assert "p_enrichment" in r.stdout and "q_value" in r.stdout


def test_alice_command(tmp_path, gen):
    pytest.importorskip("vdjmatch")
    p = tmp_path / "s.tsv"
    _write_airr(p, gen, seed=41)
    r = runner.invoke(app, ["alice", str(p), "--locus", "TRB"])
    assert r.exit_code == 0, r.stdout
    assert "pgen_ball" in r.stdout and "q_value" in r.stdout
