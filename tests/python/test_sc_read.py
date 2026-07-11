"""10x / AIRR-Cell ingestion and AIRR Data File export."""
from __future__ import annotations

import hashlib

import polars as pl
import pytest

from vdjtools import sc
from vdjtools.sc.read import SC_COLUMNS

_ALL_CONTIG_HEADER = (
    "barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,"
    "c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,"
    "raw_consensus_id"
)


def _row(barcode, contig, chain, cdr3, *, is_cell="True", hc="True", prod="True",
         v="V", d="", j="J", reads=100, umis=10, clone="clone1", cons="cons1"):
    return (f"{barcode},{is_cell},{contig},{hc},300,{chain},{v},{d},{j},C,True,{prod},"
            f"{cdr3},TGTGCT,{reads},{umis},{clone},{cons}")


def _write(tmp_path, lines, name="all_contig.csv"):
    p = tmp_path / name
    p.write_text(_ALL_CONTIG_HEADER + "\n" + "\n".join(lines) + "\n")
    return p


def test_read_10x_basic_mapping(tmp_path):
    lines = [
        _row("bc1", "bc1_c1", "TRB", "CASSABC", v="TRBV20-1", d="TRBD1", j="TRBJ2-1",
             reads=100, umis=10, cons="cons_b"),
        _row("bc1", "bc1_c2", "TRA", "CAVLDS", v="TRAV1-2", j="TRAJ33",
             reads=80, umis=8, cons="cons_a"),
    ]
    df = sc.read_10x(_write(tmp_path, lines))
    assert df.columns == SC_COLUMNS
    assert df.height == 2
    assert set(df["locus"]) == {"TRA", "TRB"}
    assert df["cell_id"].to_list() == ["bc1", "bc1"]
    # 10x cdr3 IS the junction (anchors incl.) → mapped straight to cdr3_aa.
    assert set(df["cdr3_aa"]) == {"CASSABC", "CAVLDS"}
    assert df["duplicate_count"].to_list() == [100, 80]  # reads
    assert df["umi_count"].to_list() == [10, 8]          # umis
    assert df.filter(pl.col("locus") == "TRB")["v_call"][0] == "TRBV20-1"


def test_read_10x_filters(tmp_path):
    lines = [
        _row("keep", "k1", "TRB", "CASSKEEP"),
        _row("noncell", "n1", "TRB", "CASSNO", is_cell="False"),
        _row("lowconf", "l1", "TRB", "CASSLC", hc="False"),
        _row("nonprod", "p1", "TRB", "CASSNP", prod="False"),
        _row("multi", "m1", "Multi", "CASSMU"),
        _row("nocons", "x1", "TRB", "CASSNC", cons=""),   # null raw_consensus_id
    ]
    df = sc.read_10x(_write(tmp_path, lines))
    assert df["cell_id"].to_list() == ["keep"]


def test_read_10x_require_flags_off(tmp_path):
    lines = [
        _row("cellA", "a1", "TRB", "CASSA", is_cell="False", hc="False"),
    ]
    # With both flags off, is_cell / high_confidence are ignored (productive still filters).
    df = sc.read_10x(_write(tmp_path, lines), require_cell=False, require_high_conf=False)
    assert df["cell_id"].to_list() == ["cellA"]


def test_read_10x_consensus_join_overrides_calls(tmp_path):
    contig = _write(tmp_path, [
        _row("bc1", "bc1_c1", "TRB", "CASSABC", v="TRBV_CONTIG", j="TRBJ_CONTIG",
             clone="cloneA", cons="consA"),
    ])
    cons = tmp_path / "consensus.csv"
    cons.write_text(
        "clonotype_id,consensus_id,chain,cdr3,cdr3_nt,v_gene,d_gene,j_gene\n"
        "cloneA,consA,TRB,CASSABC,TGTGCT,TRBV_CONSENSUS,TRBD1,TRBJ_CONSENSUS\n"
    )
    df = sc.read_10x(contig, cons)
    assert df["v_call"][0] == "TRBV_CONSENSUS"   # consensus call overrides contig call
    assert df["j_call"][0] == "TRBJ_CONSENSUS"


def test_write_airr_cell_receptor_hash(tmp_path):
    df = pl.DataFrame(
        {
            "cell_id": ["bc1", "bc1"],
            "sequence_id": ["b1", "a1"],
            "locus": ["TRB", "TRA"],
            "v_call": ["TRBV1", "TRAV1"], "d_call": [None, None],
            "j_call": ["TRBJ1", "TRAJ1"], "c_call": [None, None],
            "cdr3_aa": ["CASSZZZ", "CAVAAA"], "cdr3_nt": ["TGT", "TGT"],
            "duplicate_count": [100, 80], "umi_count": [10, 8],
            "clone_id": ["c1", "c1"],
        }
    )
    out = tmp_path / "cells.yaml"
    sc.write_airr_cell(df, out)
    yaml = pytest.importorskip("yaml")
    doc = yaml.safe_load(out.read_text())
    assert len(doc["Cell"]) == 1
    assert doc["Cell"][0]["cell_id"] == "bc1"
    assert doc["Cell"][0]["virtual_pairing"] is False
    rec = doc["Receptor"][0]
    # domain_1 = heavy(beta) junction, domain_2 = light(alpha) junction.
    assert rec["receptor_variable_domain_1_aa"] == "CASSZZZ"
    assert rec["receptor_variable_domain_1_locus"] == "TRB"
    assert rec["receptor_variable_domain_2_aa"] == "CAVAAA"
    expect = hashlib.sha256(("CASSZZZ" + "CAVAAA").encode()).hexdigest()
    assert rec["receptor_hash"] == expect
    assert rec["receptor_type"] == "TCR"


def test_read_airr_cell_roundtrip(tmp_path):
    df = pl.DataFrame(
        {
            "cell_id": ["bc1"], "sequence_id": ["s1"], "locus": ["TRB"],
            "v_call": ["TRBV1"], "d_call": [None], "j_call": ["TRBJ1"], "c_call": [None],
            "cdr3_aa": ["CASS"], "cdr3_nt": ["TGT"],
            "duplicate_count": [5], "umi_count": [2], "clone_id": ["c1"],
        }
    )
    tsv = tmp_path / "cells.tsv"
    df.write_csv(tsv, separator="\t")
    back = sc.read_airr_cell(tsv)
    assert back.columns == SC_COLUMNS
    assert back["cell_id"].to_list() == ["bc1"]
    assert back["duplicate_count"].to_list() == [5]


def test_read_airr_cell_prefers_junction_over_imgt_cdr3(tmp_path):
    """When both junction_aa (anchors included) and IMGT cdr3_aa (excluded) are present,
    the reader takes the junction — the canonical cdr3_aa=junction convention."""
    tsv = tmp_path / "cells.tsv"
    pl.DataFrame({
        "cell_id": ["bc1"], "locus": ["TRB"], "v_call": ["TRBV1"], "j_call": ["TRBJ1"],
        "junction_aa": ["CASSLGQAYEQYF"], "cdr3_aa": ["ASSLGQAYEQY"],   # 13-mer vs 11-mer
        "junction": ["TGTGCC"], "cdr3": ["GCC"],
    }).write_csv(tsv, separator="\t")
    back = sc.read_airr_cell(tsv)
    assert back["cdr3_aa"].to_list() == ["CASSLGQAYEQYF"]   # junction, not the IMGT 11-mer
    assert back["cdr3_nt"].to_list() == ["TGTGCC"]
