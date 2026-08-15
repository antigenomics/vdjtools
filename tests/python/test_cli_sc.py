"""``vdjtools sc <sub>`` — every subcommand, through the real argv.

Deliberately invoked through the CLI runner rather than `--help`: typer short-circuits
`--help` before argument parsing, so a help smoke test cannot catch a broken option or a
renamed flag. These assert the wiring (exit codes, what lands on stdout, files written),
not the analysis -- that is covered by the library-level tests.
"""
from __future__ import annotations

import polars as pl
import pytest
from typer.testing import CliRunner

from vdjtools import sc
from vdjtools.cli import app
from vdjtools.sc.read import SC_COLUMNS

runner = CliRunner()


def _run(*args):
    return runner.invoke(app, list(args))


def _ok(result):
    assert result.exit_code == 0, f"exit {result.exit_code}\n{result.output}\n{result.exception!r}"
    return result


@pytest.fixture(scope="module")
def airr_tsv(tmp_path_factory):
    """Two complete alpha/beta cells, written as an AIRR Rearrangement TSV."""
    cells = pl.DataFrame({
        "cell_id": ["c1", "c1", "c2", "c2"],
        "sequence_id": ["a", "b", "c", "d"],
        "locus": ["TRA", "TRB", "TRA", "TRB"],
        "v_call": ["TRAV1-2*01", "TRBV20-1*01", "TRAV12-1*01", "TRBV9*01"],
        "d_call": [None, "TRBD1*01", None, None],
        "j_call": ["TRAJ33*01", "TRBJ2-7*01", "TRAJ8*01", "TRBJ2-3*01"],
        "c_call": ["TRAC", "TRBC2", "TRAC", "TRBC2"],
        "junction_aa": ["CAVRDSNYQLIW", "CASSLGQAYEQYF", "CAVNTGGFKTIF", "CASSVDTQYF"],
        "junction_nt": ["TGTGCC", "TGTGCCAGC", "TGTGCTGTG", "TGTGCCAGCAGC"],
        "duplicate_count": [120, 340, 88, 150],
        "umi_count": [6, 14, 4, 9],
        "clone_id": ["x", "x", "y", "y"],
        "productive": [True, True, True, True],
    }).select(SC_COLUMNS)
    path = tmp_path_factory.mktemp("sc") / "airr_rearrangement.tsv"
    return sc.write_airr(cells, path)


def test_convert_emits_the_canonical_frame(airr_tsv):
    out = _ok(_run("sc", "convert", str(airr_tsv))).output
    assert out.splitlines()[0].split("\t") == SC_COLUMNS


def test_convert_airr_flag_switches_to_the_airr_spelling(airr_tsv):
    header = _ok(_run("sc", "convert", str(airr_tsv), "--airr")).output.splitlines()[0]
    assert "junction" in header.split("\t") and "junction_nt" not in header.split("\t")
    assert "consensus_count" in header.split("\t")


def test_qc_reports_the_multiplicity_quadrants(airr_tsv):
    result = _ok(_run("sc", "qc", str(airr_tsv)))
    assert "n_light\tn_heavy\tcell_count" in result.output
    assert "2 cells, 4 contigs" in result.output


def test_pair_emits_one_row_per_receptor(airr_tsv):
    lines = _ok(_run("sc", "pair", str(airr_tsv))).output.strip().splitlines()
    assert lines[0].startswith("cell_id\tpair_id\talpha_v_call")
    assert len(lines) == 3          # header + two paired cells


def test_pair_can_add_mispairing_flags(airr_tsv):
    header = _ok(_run("sc", "pair", str(airr_tsv), "--flag-mispairing")).output.splitlines()[0]
    assert "mispairing_flag" in header and "mispairing_reason" in header


def test_pair_rejects_an_unknown_locus_pair(airr_tsv):
    result = _run("sc", "pair", str(airr_tsv), "--locus-pair", "TRA_TRA")
    assert result.exit_code != 0


def test_pgen_adds_the_three_probability_columns(airr_tsv):
    header = _ok(_run("sc", "pgen", str(airr_tsv))).output.splitlines()[0].split("\t")
    assert {"pgen_alpha", "pgen_beta", "pgen_paired"} <= set(header)


def test_export_writes_each_plain_text_target(airr_tsv, tmp_path):
    for target, name in [("airr", "a.tsv"), ("screpertoire", "s.tsv"),
                         ("screpertoire-10x", "s.csv"), ("airr-cell", "cells.yaml")]:
        out = tmp_path / name
        _ok(_run("sc", "export", str(airr_tsv), "--to", target, "-o", str(out)))
        assert out.exists() and out.stat().st_size > 0


def test_export_scirpy_writes_an_h5ad(airr_tsv, tmp_path):
    pytest.importorskip("scirpy")
    out = tmp_path / "vdj.h5ad"
    _ok(_run("sc", "export", str(airr_tsv), "--to", "scirpy", "-o", str(out)))
    import anndata as ad
    assert "airr" in ad.read_h5ad(out).obsm


def test_export_dandelion_writes_an_h5ddl(airr_tsv, tmp_path):
    pytest.importorskip("dandelion")
    out = tmp_path / "vdj.h5ddl"
    _ok(_run("sc", "export", str(airr_tsv), "--to", "dandelion", "-o", str(out)))
    # Readable back without dandelion.
    assert sc.read_h5ddl(out).height == 4


def test_export_rejects_an_unknown_target(airr_tsv, tmp_path):
    result = _run("sc", "export", str(airr_tsv), "--to", "seurat", "-o", str(tmp_path / "x"))
    assert result.exit_code != 0
