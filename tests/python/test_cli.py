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


@pytest.mark.parametrize('command', ['pgen', 'generate'])
@pytest.mark.parametrize('organism', ['human', 'mouse'])
def test_model_commands_forward_organism(tmp_path, monkeypatch, command, organism):
    import importlib
    import vdjtools.model as models
    from vdjtools.model import native
    marker, calls = object(), []
    def load(locus, *, source, organism):
        calls.append((locus, source, organism))
        return marker
    monkeypatch.setattr(models, 'load_bundled', load)
    def sample(model, n, **kwargs):
        assert model is marker and kwargs['engine'] == 'native'
        return pl.DataFrame({'junction_aa': ['CASSF'] * n})
    monkeypatch.setattr(importlib.import_module('vdjtools.model.generate'), 'generate', sample)
    def score(model, sequences, *args, **kwargs):
        assert model is marker
        return [.25] * len(sequences)
    monkeypatch.setattr(native, 'pgen_aa_batch', score)
    source = tmp_path / 'queries.tsv'
    source.write_text('junction_aa\nCASSF\n')
    arguments = ([command, str(source)] if command == 'pgen' else [command, '-n', '2'])
    arguments += ['--model', 'TRA', '-o', str(tmp_path / 'out.tsv')]
    if organism == 'mouse':
        arguments += ['--organism', 'mouse']
    result = runner.invoke(app, arguments)
    assert result.exit_code == 0, result.output
    assert calls == [('TRA', 'learned', organism)]


@pytest.mark.parametrize('command', ['pgen', 'generate'])
@pytest.mark.parametrize('organism,message', [('mouse', 'no bundled mouse TRG model'),
                                            ('rabbit', '--organism must be human or mouse')])
def test_unsupported_mouse_model_is_a_friendly_cli_error(tmp_path, monkeypatch, command, organism, message):
    import vdjtools.model as models
    def missing(locus, *, source, organism):
        assert locus == 'TRG' and organism == 'mouse'
        raise FileNotFoundError('no bundled mouse TRG model')
    monkeypatch.setattr(models, 'load_bundled', missing)
    source = tmp_path / 'queries.tsv'
    source.write_text('junction_aa\nCASSF\n')
    arguments = [command, str(source)] if command == 'pgen' else [command]
    result = runner.invoke(app, arguments + ['--model', 'TRG', '--organism', organism])
    assert result.exit_code == 1
    assert 'error: ' + message in result.output
    assert isinstance(result.exception, SystemExit)


@pytest.mark.parametrize("source,engine,custom,expected_engine", [
    ("learned", None, False, "native"),
    ("learned", "reference", False, "reference"),
    ("learned", "native", False, "native"),
    ("olga", None, False, "reference"),
    ("arda", None, False, "reference"),
    ("olga", "native", False, "native"),
    ("learned", None, True, "reference"),
    ("learned", "native", True, "native"),
])
def test_generate_cli_engine_policy(tmp_path, toy_model, source, engine, custom, expected_engine):
    from vdjtools.model import load_bundled
    from vdjtools.model.generate import generate
    from vdjtools.model.io import load_model, save_model

    if custom:
        path = tmp_path / "model"
        save_model(toy_model, path)
        model = load_model(path)
        arguments = ["--model-path", str(path)]
    else:
        model = load_bundled("TRG", source)
        arguments = ["--model", "TRG"]
        if source != "learned":
            arguments += ["--source", source]
    if engine is not None:
        arguments += ["--engine", engine]
    out = tmp_path / "draws.parquet"
    result = runner.invoke(app, ["generate", *arguments, "-n", "64", "--seed", "7", "-o", str(out)])
    assert result.exit_code == 0, result.output
    expected = generate(model, 64, seed=7, engine=expected_engine, threads=1)
    assert pl.read_parquet(out).equals(expected)


@pytest.mark.parametrize('organism,locus',
                         [('human', locus) for locus in ('TRA', 'TRB', 'TRG', 'TRD', 'IGH', 'IGK', 'IGL')]
                         + [('mouse', locus) for locus in ('TRA', 'TRB', 'IGH', 'IGK', 'IGL')])
def test_generate_native_cli_matches_library(tmp_path, organism, locus):
    from vdjtools.model import load_bundled
    from vdjtools.model.generate import generate

    out = tmp_path / "native.parquet"
    result = runner.invoke(app, ["generate", "-m", locus, "--organism", organism,
                                 "-n", "64", "--seed", "7",
                                 "--engine", "native", "--threads", "1", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert pl.read_parquet(out).equals(generate(load_bundled(locus, organism=organism), 64, seed=7,
                                               engine="native", threads=1))


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


def test_correct_vj_writes_usage_and_corrected_tables(tmp_path, gen):
    """`correct-vj` is the CLI face of the V/J-usage batch correction (Vlasova et al. 2026).

    Both label forms are exercised, because the TSV branch matches on the FILE STEM and that is
    the part a caller gets wrong silently -- a metadata table keyed on anything else simply finds
    no batch for any sample.
    """
    paths = []
    for i in range(4):
        p = tmp_path / f"s{i}.tsv"
        _write_airr(p, gen, seed=70 + i)
        paths.append(str(p))

    usage, outdir = tmp_path / "usage.tsv", tmp_path / "corrected"
    r = runner.invoke(app, ["correct-vj", *paths, "-b", "A,A,B,B",
                            "--transform", "sigmoid",
                            "--usage-out", str(usage), "--outdir", str(outdir)])
    assert r.exit_code == 0, r.stdout
    assert pl.read_csv(usage, separator="\t").height > 0
    assert sorted(p.name for p in outdir.glob("*.tsv")) == [f"s{i}.corrected.tsv" for i in range(4)]

    meta = tmp_path / "meta.tsv"
    pl.DataFrame({"sample_id": [f"s{i}" for i in range(4)],
                  "batch": ["A", "A", "B", "B"]}).write_csv(meta, separator="\t")
    r = runner.invoke(app, ["correct-vj", *paths, "-b", str(meta), "--usage-out", str(tmp_path / "u2.tsv")])
    assert r.exit_code == 0, r.stdout

    # one sample is not a between-batch quantity, and a short label list is a caller error
    assert runner.invoke(app, ["correct-vj", paths[0], "-b", "A"]).exit_code != 0
    assert runner.invoke(app, ["correct-vj", *paths, "-b", "A,B"]).exit_code != 0


@pytest.fixture(scope="module")
def tiny_corpus(tmp_path_factory):
    """One small artifact, built once: the CLI tests are about the command, not the fit."""
    from vdjtools.signature.corpus import synthesize

    art, _ = synthesize("memory", loci=("TRG",), n_samples=24, size=80, seed=4, n_components=4)
    return art.save(tmp_path_factory.mktemp("corpus") / "tiny")


def _one_sample(path):
    path.write_text("junction_aa\tv_call\tj_call\tduplicate_count\n"
                    + "".join(f"CASS{'ACDEFGHIKLMNPQRSTVWY'[i % 20] * 3}YEQYF\tTRBV20-1\t"
                              f"TRBJ2-2\t{i + 1}\n" for i in range(60)))
    return path


def test_signature_requires_a_corpus(tmp_path):
    """There is no default corpus, deliberately.

    A silently chosen rotation makes two matrices look comparable when they are not, which is the
    mixed-units failure the rewrite exists to end.
    """
    src = _one_sample(tmp_path / "s1.tsv")
    assert CliRunner().invoke(app, ["signature", str(src)]).exit_code != 0


def test_describe_names_exactly_the_columns_this_invocation_emits(tmp_path, tiny_corpus):
    """The one output whose entire job is to be right about the width.

    In 3.18.0 ``--describe`` printed all 688 columns while the command emitted 161 -- it named
    columns you would not get. So it is resolved against --corpus, --components AND --columns, and
    compared here against the real emitted header rather than against a layout constant.
    """
    src = _one_sample(tmp_path / "s1.tsv")
    runner = CliRunner()

    for extra in ([], ["--components", "2"]):
        desc, emit = tmp_path / "d.tsv", tmp_path / "e.tsv"
        r1 = runner.invoke(app, ["signature", "--corpus", str(tiny_corpus), "--describe",
                                 "-o", str(desc), *extra])
        assert r1.exit_code == 0, r1.output
        r2 = runner.invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src),
                                 "-o", str(emit), *extra])
        assert r2.exit_code == 0, r2.output
        described = pl.read_csv(desc, separator="\t")["column"].to_list()
        emitted = [c for c in pl.read_csv(emit, separator="\t").columns if c != "sample_id"]
        assert described == emitted, f"--describe disagrees with the output for {extra}"

    # and it carries the two things a caller needs per column
    desc = pl.read_csv(tmp_path / "d.tsv", separator="\t")
    assert set(desc["kind"].unique()) <= {"rotated", "channel"}
    assert set(desc["support"].unique()) <= set(__import__(
        "vdjtools.signature.layout", fromlist=["SUPPORTS"]).SUPPORTS)


def test_describe_honours_a_column_subset(tmp_path, tiny_corpus):
    want = tmp_path / "cols.txt"
    want.write_text("vsig:pc:TRG:PC01\nvsig:cov:TRG:cstar\n")
    out = tmp_path / "d2.tsv"
    r = CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), "--describe",
                                 "--columns", str(want), "-o", str(out)])
    assert r.exit_code == 0, r.output
    assert pl.read_csv(out, separator="\t")["column"].to_list() == [
        "vsig:pc:TRG:PC01", "vsig:cov:TRG:cstar"]


def test_the_coverage_target_is_a_runtime_choice_not_a_baked_constant(tmp_path, tiny_corpus):
    """The target used to live in a scaling artifact, which is how one reference ended up
    standardising tissue samples to a level measured on blood."""
    src = _one_sample(tmp_path / "s1.tsv")

    def run(tag, *extra):
        out = tmp_path / f"o_{tag}.tsv"
        r = CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src),
                                     "-o", str(out), *extra])
        assert r.exit_code == 0, r.output
        return pl.read_csv(out, separator="\t")

    shared, own = run("min"), run("own", "--cstar-target", "own")
    fixed = run("fixed", "--cstar-target", "0.5")
    for frame in (shared, own, fixed):
        # coverage is emitted whatever the target, because it is a property of the sample
        assert frame["vsig:cov:TRB:cstar"][0] is not None
    assert CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src),
                                    "--cstar-target", "not-a-number"]).exit_code != 0


def test_whatever_the_bound_clamps_is_reported(tmp_path, tiny_corpus):
    """A bound that edits a number the caller cannot see is the failure winsor_frac prevents."""
    src = _one_sample(tmp_path / "s1.tsv")

    def frac(mode):
        out = tmp_path / f"w_{mode}.tsv"
        r = CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src),
                                     "--winsorize", mode, "-o", str(out)])
        assert r.exit_code == 0, r.output
        return pl.read_csv(out, separator="\t")["vsig:qc:-:winsor_frac"][0]

    assert frac("none") == 0.0
    assert frac("features") >= 0.0
    assert CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src),
                                    "--winsorize", "nonsense"]).exit_code != 0


@pytest.mark.parametrize("flag", [
    ["--preset", "classify"], ["--tier", "core"], ["--pgen-n-max", "200"],
    ["--cstar", "0.2"], ["--threads", "4"], ["--channels"], ["--standardize", "none"],
    ["--scale", "blood"], ["--clip", "8"], ["--squash", "soft"],
])
def test_every_removed_flag_is_rejected(tmp_path, tiny_corpus, flag):
    src = _one_sample(tmp_path / "s1.tsv")
    r = CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src), *flag])
    assert r.exit_code != 0, f"{flag} was accepted"


def test_the_presets_command_is_gone():
    from typer.main import get_command

    assert "presets" not in get_command(app).commands
    assert "corpus" in get_command(app).commands


def test_the_preflight_line_names_the_corpus_and_its_hash(tmp_path, tiny_corpus):
    """A saved matrix's provenance should be in the log even without the manifest to hand."""
    src = _one_sample(tmp_path / "s1.tsv")
    r = CliRunner().invoke(app, ["signature", "--corpus", str(tiny_corpus), str(src),
                                 "-o", str(tmp_path / "o.tsv")])
    assert r.exit_code == 0, r.output
    assert "corpus memory" in r.output and "winsorize=features" in r.output and "k=" in r.output


def test_diversity_process_backend(tmp_path, gen):
    from vdjtools.cores import available_cores
    if available_cores() < 2:
        pytest.skip('requires two allocated cores')
    a, b = tmp_path / 'a.tsv', tmp_path / 'b.tsv'
    _airr_sample(gen[:100], 1).write_csv(a, separator="\t")
    _airr_sample(gen[100:200], 2).write_csv(b, separator="\t")
    args = ['diversity', str(a), str(b), '--threads', '2']
    threaded = runner.invoke(app, args)
    processes = runner.invoke(app, args + ['--backend', 'process'])
    assert threaded.exit_code == processes.exit_code == 0, processes.output
    assert threaded.stdout == processes.stdout
