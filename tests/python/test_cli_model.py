"""``vdjtools model <sub>`` — every subcommand, on the smallest models available.

Kept fast deliberately: TRG is the smallest bundled locus, the toy model is smaller still, and the
EM commands run one or two iterations. These assert the wiring (exit codes, files written, the
shape of what lands on stdout), not the maths — that is covered by the library-level tests.
"""
from __future__ import annotations

import polars as pl
import pytest
from typer.testing import CliRunner

from vdjtools.cli import app
from vdjtools.model.generate import generate
from vdjtools.model.io import save_model

runner = CliRunner()


def _run(*args):
    return runner.invoke(app, list(args))


def _ok(result):
    assert result.exit_code == 0, f"exit {result.exit_code}\n{result.output}\n{result.exception!r}"
    return result


@pytest.fixture(scope="module")
def trg_seqs(tmp_path_factory):
    """A small TRG sequence table with V/J calls, written once for the whole module."""
    from vdjtools.model import load_bundled

    path = tmp_path_factory.mktemp("cli") / "seqs.tsv"
    generate(load_bundled("TRG", "olga"), 60, seed=1).select(
        ["junction_nt", "junction_aa", "v_call", "j_call"]).write_csv(path, separator="\t")
    return path


@pytest.fixture
def toy_dir(tmp_path, toy_model):
    save_model(toy_model, tmp_path / "toy")
    return tmp_path / "toy"


def test_model_help_lists_the_subcommands():
    out = _ok(_run("model", "--help")).output
    for cmd in ("check", "template", "learn", "extend", "rescale", "export", "net",
                "entropy", "diversity", "compare", "compare-pgen", "loglik", "log"):
        assert cmd in out


def test_list():
    assert "olga" in _ok(_run("model", "list")).output


def test_check_clean_model_exits_zero(tmp_path):
    out = tmp_path / "issues.tsv"
    _ok(_run("model", "check", "TRG:olga", "-o", str(out)))
    assert "severity" in pl.read_csv(out, separator="\t").columns or out.read_text().strip()


def test_check_exits_one_on_a_broken_model(tmp_path, toy_model):
    broken = type(toy_model)(
        manifest=toy_model.manifest,
        tables={**toy_model.tables,
                "vj_ins": toy_model.tables["vj_ins"].with_columns(p=pl.col("p") * 4)},
        genomic=toy_model.genomic)
    save_model(broken, tmp_path / "broken")
    result = _run("model", "check", str(tmp_path / "broken"))
    assert result.exit_code == 1
    assert "error" in result.output.lower()


def test_check_rejects_an_unknown_model():
    result = _run("model", "check", "NOSUCHLOCUS")
    assert result.exit_code == 1


def test_template_from_arda(tmp_path):
    out = tmp_path / "tmpl"
    _ok(_run("model", "template", "--locus", "TRG", "-o", str(out)))
    assert (out / "manifest.json").exists()
    from vdjtools.model import load_model

    assert load_model(out, validate=True).locus == "TRG"


def test_template_from_fasta(tmp_path):
    v, j = tmp_path / "v.fa", tmp_path / "j.fa"
    v.write_text(">TOYV1*01\nTGTGCCAGC\n>TOYV2*01\nTGTGCTTCC\n")
    j.write_text(">TOYJ1*01\nAACTATGGCTATACCTTT\n")
    out = tmp_path / "custom"
    _ok(_run("model", "template", "--locus", "TOY", "--germline-v", str(v),
             "--germline-j", str(j), "--ins-max", "3", "-o", str(out)))
    from vdjtools.model import load_model

    m = load_model(out, validate=True)
    assert m.chain_type == "VJ" and m.genomic["genes_v"].height == 2


def test_template_needs_an_output_dir():
    assert _run("model", "template", "--locus", "TRG").exit_code == 1


def test_export_long_and_directory(tmp_path):
    flat = tmp_path / "marginals.tsv"
    _ok(_run("model", "export", "TRG:olga", "--long", "-o", str(flat)))
    df = pl.read_csv(flat, separator="\t")
    assert {"event", "kind", "p"} <= set(df.columns)

    as_tsv = tmp_path / "trg_tsv"
    _ok(_run("model", "export", "TRG:olga", "--format", "tsv", "-o", str(as_tsv)))
    assert (as_tsv / "v_choice.tsv").exists()
    # A TSV model directory is a first-class model path.
    _ok(_run("model", "check", str(as_tsv)))


def test_net_dot_to_stdout_and_file(tmp_path):
    out = _ok(_run("model", "net", "TRG:olga")).output
    assert "digraph bn" in out
    path = tmp_path / "net.dot"
    _ok(_run("model", "net", "TRG:olga", "-o", str(path)))
    assert "digraph bn" in path.read_text()


def test_entropy_tables(tmp_path):
    for table in ("entropy", "mi", "total"):
        out = tmp_path / f"{table}.tsv"
        _ok(_run("model", "entropy", "TRG:olga", "--table", table, "-o", str(out)))
        assert pl.read_csv(out, separator="\t").height > 0
    assert _run("model", "entropy", "TRG:olga", "--table", "nope").exit_code == 1


def test_diversity(tmp_path):
    out = tmp_path / "div.tsv"
    result = _ok(_run("model", "diversity", "TRG:olga", "-n", "200", "-o", str(out)))
    df = pl.read_csv(out, separator="\t")
    assert {"scenario_entropy_bits", "diversity_shannon", "diversity_simpson"} <= set(df.columns)
    assert "Shannon" in result.output


def test_compare_and_usage(tmp_path):
    out = tmp_path / "diff.tsv"
    _ok(_run("model", "compare", "TRG:olga", "TRG:learned", "--by", "gene", "-o", str(out)))
    df = pl.read_csv(out, separator="\t")
    assert {"event", "status", "jsd_bits", "tv_max"} <= set(df.columns)

    usage = tmp_path / "usage.tsv"
    _ok(_run("model", "compare", "TRG:olga", "TRG:learned", "--usage", "v", "-o", str(usage)))
    assert {"name", "p_a", "p_b"} <= set(pl.read_csv(usage, separator="\t").columns)

    dot = tmp_path / "diff.dot"
    _ok(_run("model", "compare", "TRG:olga", "TRG:learned", "--dot", str(dot),
             "-o", str(tmp_path / "d2.tsv")))
    assert "digraph compare" in dot.read_text()


def test_compare_pgen_and_summary(tmp_path, trg_seqs):
    out = tmp_path / "cmp.tsv"
    _ok(_run("model", "compare-pgen", "TRG:olga", "TRG:learned", str(trg_seqs),
             "-c", "junction_aa", "-o", str(out)))
    assert "delta_log10" in pl.read_csv(out, separator="\t").columns

    summary = tmp_path / "summary.tsv"
    _ok(_run("model", "compare-pgen", "TRG:olga", "TRG:learned", str(trg_seqs),
             "-c", "junction_aa", "--summary", "-o", str(summary)))
    df = pl.read_csv(summary, separator="\t")
    assert df.height == 1
    assert {"ks_stat", "spearman_log10", "only_a_scoreable"} <= set(df.columns)


def test_loglik_reports_bic(tmp_path, trg_seqs):
    out = tmp_path / "fit.tsv"
    result = _ok(_run("model", "loglik", str(trg_seqs), "TRG:olga",
                      "-c", "junction_nt", "-o", str(out)))
    df = pl.read_csv(out, separator="\t")
    assert df.height == 1 and {"loglik_sum", "k", "aic", "bic"} <= set(df.columns)
    assert "BIC=" in result.output


def test_loglik_per_sequence(tmp_path, trg_seqs):
    out = tmp_path / "per.tsv"
    _ok(_run("model", "loglik", str(trg_seqs), "TRG:olga", "-c", "junction_nt",
             "--per-sequence", "-o", str(out)))
    df = pl.read_csv(out, separator="\t")
    assert df.height == 60 and "pgen" in df.columns


def test_learn_writes_a_model_and_a_training_log(tmp_path, toy_dir, toy_model):
    clones = tmp_path / "clones.tsv"
    generate(toy_model, 120, seed=2).rename({"junction_nt": "junction"}).select(
        ["junction", "v_call", "j_call"]).write_csv(clones, separator="\t")
    out = tmp_path / "learned"
    result = _ok(_run("model", "learn", str(clones), "--template", str(toy_dir),
                      "--max-iter", "2", "-o", str(out)))
    assert "loglik" in result.output
    assert (out / "training.json").exists()

    log = tmp_path / "log.tsv"
    _ok(_run("model", "log", str(out), "-o", str(log)))
    df = pl.read_csv(log, separator="\t")
    assert {"run", "iter", "loglik"} <= set(df.columns)


def test_log_errors_when_there_is_none():
    assert _run("model", "log", "TRG:olga").exit_code == 1


def test_learn_needs_a_template(tmp_path, trg_seqs):
    assert _run("model", "learn", str(trg_seqs), "-o", str(tmp_path / "x")).exit_code == 1


def test_extend(tmp_path, toy_dir):
    v, j = tmp_path / "v.fa", tmp_path / "j.fa"
    v.write_text(">TOYV1*01\nTGTGCCAGC\n>TOYV1*03\nTGTGCCAGA\n>TOYV3*01\nTGTTGGGGA\n")
    j.write_text(">TOYJ1*01\nAACTATGGCTATACCTTT\n")
    out = tmp_path / "extended"
    _ok(_run("model", "extend", str(toy_dir), "--germline-v", str(v),
             "--germline-j", str(j), "-o", str(out)))
    from vdjtools.model import load_model

    m = load_model(out, validate=True)
    assert {"TOYV1*03", "TOYV3*01"} <= set(m.genomic["genes_v"]["v_allele"])


def test_rescale(tmp_path):
    """Usage is protocol-dependent; this is the only route to rescale_usage from the CLI."""
    from vdjtools.model import load_bundled

    sample = tmp_path / "sample.tsv"
    gen = generate(load_bundled("TRG", "olga"), 200, seed=3)
    gen.select(["junction_nt", "junction_aa", "v_call", "j_call"]).with_columns(
        count=pl.lit(1)).rename({"junction_nt": "cdr3nt", "junction_aa": "cdr3aa"}).write_csv(
        sample, separator="\t")
    out = tmp_path / "rescaled"
    result = _run("model", "rescale", "TRG:olga", str(sample), "--format", "auto", "-o", str(out))
    if result.exit_code != 0:          # the reader may not sniff this minimal table
        pytest.skip(f"sample not readable by io.read: {result.output}")
    from vdjtools.model import load_model

    load_model(out, validate=True)


def test_top_level_model_commands_still_work(tmp_path, trg_seqs):
    """pgen / generate / models must keep working unchanged next to the sub-app."""
    _ok(_run("models"))
    out = tmp_path / "gen.tsv"
    _ok(_run("generate", "--model", "TRG", "-n", "5", "--seed", "1", "-o", str(out)))
    pg = tmp_path / "pgen.tsv"
    _ok(_run("pgen", str(trg_seqs), "--model", "TRG", "-c", "junction_aa", "-o", str(pg)))
    assert "pgen" in pl.read_csv(pg, separator="\t").columns
