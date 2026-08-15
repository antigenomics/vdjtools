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


def _table_header(result):
    """The TSV header line, ignoring any progress/info lines printed before it."""
    for line in result.output.splitlines():
        if line.startswith("cell_id\t"):
            return line.split("\t")
    raise AssertionError(f"no table in output:\n{result.output}")


def test_pgen_adds_the_three_probability_columns(airr_tsv):
    header = _table_header(_ok(_run("sc", "pgen", str(airr_tsv))))
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


# ------------------------------------------------------------------- format dispatch

def test_fmt_is_sniffed_from_the_header_not_the_filename(tmp_path, airr_tsv):
    """Filenames are not a contract; people rename exports."""
    from vdjtools.cli import _sniff_sc

    renamed = tmp_path / "some_contigs_thing.txt"
    renamed.write_bytes(airr_tsv.read_bytes())
    assert _sniff_sc(renamed) == "airr"
    _ok(_run("sc", "convert", str(renamed)))


def test_a_10x_csv_is_recognised_by_its_barcode_and_contig_id(tmp_path):
    from vdjtools.cli import _sniff_sc

    csv = tmp_path / "anything.csv"
    csv.write_text("barcode,contig_id,chain,cdr3,cdr3_nt,reads,umis,is_cell,high_confidence,"
                   "productive\nAAA-1,c1,TRB,CASSF,TGT,10,2,True,True,True\n")
    assert _sniff_sc(csv) == "10x"


def test_fmt_can_be_forced(airr_tsv):
    _ok(_run("sc", "convert", str(airr_tsv), "--fmt", "airr"))


def test_an_unknown_fmt_is_rejected(airr_tsv):
    assert _run("sc", "convert", str(airr_tsv), "--fmt", "bogus").exit_code != 0


def test_a_bulk_table_is_rejected_with_a_useful_message(tmp_path):
    bulk = tmp_path / "bulk.tsv"
    bulk.write_text("v_call\tjunction_aa\nTRBV9*01\tCASSF\n")
    result = _run("sc", "convert", str(bulk))
    assert result.exit_code != 0
    assert "not a single-cell table" in result.output


# --------------------------------------------------------------------------- options

def test_pair_no_resolve_keeps_every_chain_combination(airr_tsv):
    _ok(_run("sc", "pair", str(airr_tsv), "--no-resolve"))


def test_mispairing_options_require_flag_mispairing(airr_tsv):
    for opt in (["--drop-mispaired"], ["--max-slaves-per-master", "3"]):
        result = _run("sc", "pair", str(airr_tsv), *opt)
        assert result.exit_code != 0
        assert "need --flag-mispairing" in result.output


def test_pgen_reports_how_many_receptors_scored(airr_tsv):
    assert "scored 2/2 receptors" in _ok(_run("sc", "pgen", str(airr_tsv))).output


def test_pgen_no_resolve_genes_is_accepted(airr_tsv):
    _ok(_run("sc", "pgen", str(airr_tsv), "--no-resolve-genes"))


def test_pgen_loci_can_be_forced(airr_tsv):
    _ok(_run("sc", "pgen", str(airr_tsv), "--alpha-locus", "TRA", "--beta-locus", "TRB"))


def test_export_airr_cell_carries_the_repertoire_id(airr_tsv, tmp_path):
    out = tmp_path / "cells.yaml"
    _ok(_run("sc", "export", str(airr_tsv), "--to", "airr-cell", "-o", str(out),
             "--repertoire-id", "REP1"))
    assert "repertoire_id: REP1" in out.read_text()


def test_export_gex_only_applies_to_scirpy(airr_tsv, tmp_path):
    result = _run("sc", "export", str(airr_tsv), "--to", "airr", "-o", str(tmp_path / "a.tsv"),
                  "--gex", str(tmp_path / "nope.h5ad"))
    assert result.exit_code != 0 and "--gex applies to" in result.output


def test_export_with_gex_writes_a_mudata(airr_tsv, tmp_path):
    pytest.importorskip("scirpy")
    pytest.importorskip("mudata")
    import anndata as ad
    import numpy as np

    adata = sc.to_scirpy(sc.read_airr_cell(airr_tsv))
    gex = ad.AnnData(X=np.zeros((adata.n_obs, 3), dtype="float32"))
    gex.obs_names = adata.obs_names
    gex_path = tmp_path / "gex.h5ad"
    gex.write_h5ad(gex_path)

    out = tmp_path / "vdj.h5mu"
    _ok(_run("sc", "export", str(airr_tsv), "--to", "scirpy", "-o", str(out),
             "--gex", str(gex_path)))
    import mudata
    assert set(mudata.read_h5mu(out).mod) == {"gex", "airr"}
