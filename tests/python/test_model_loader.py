"""Phase 1a — OLGA model loader + parquet round-trip.

Proves our polars schema is a *lossless* representation of an OLGA model: we reconstruct
OLGA's own processed arrays from our long-format tables and assert equality, then round-trip
the model through parquet and assert the tables are unchanged.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from vdjtools.model import Event, EventKind, Model, from_olga
from vdjtools.model.events import validate_graph

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        str(Path(__file__).resolve().parent / "fixtures" / "olga" / "default_models"),
    )
)
olga = pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(
    not OLGA_MODELS.exists(), reason=f"OLGA bootstrap models not found at {OLGA_MODELS}"
)

TRB = OLGA_MODELS / "human_T_beta"
TRA = OLGA_MODELS / "human_T_alpha"


def _lookup(df: pl.DataFrame, keys: list[str]) -> dict:
    """{ (key cols...) : p } from a marginal table."""
    return {tuple(r[:-1]) if len(keys) > 1 else r[0]: r[-1] for r in df.select([*keys, "p"]).iter_rows()}


def test_event_graph_cycle_detected():
    good = {
        "a": Event("a", EventKind.GENE_CHOICE),
        "b": Event("b", EventKind.GENE_CHOICE, ("a",)),
    }
    validate_graph(good)  # no raise
    cyclic = {
        "a": Event("a", EventKind.GENE_CHOICE, ("b",)),
        "b": Event("b", EventKind.GENE_CHOICE, ("a",)),
    }
    with pytest.raises(ValueError, match="cycle"):
        validate_graph(cyclic)


def test_trb_vdj_shapes():
    m = from_olga(TRB, locus="TRB")  # validate() runs inside from_olga
    assert m.chain_type == "VDJ"
    counts = {name: df.height for name, df in m.tables.items()}
    assert counts == {
        "v_choice": 89,
        "j_choice": 15,
        "d_gene": 45,           # 15 J * 3 D
        "n_d": 1,
        "v_3_del": 89 * 21,
        "j_5_del": 15 * 23,
        "d_del": 3 * 21 * 21,
        "vd_ins": 31,
        "dj_ins": 31,
        "vd_dinucl": 16,
        "dj_dinucl": 16,
    }
    assert set(m.genomic) == {"genes_v", "genes_d", "genes_j"}


def test_trb_vdj_lossless_vs_olga():
    """Reconstruct OLGA's arrays from our tables and assert bit-faithfulness."""
    g = olga.GenomicDataVDJ()
    g.load_igor_genomic_data(
        str(TRB / "model_params.txt"),
        str(TRB / "V_gene_CDR3_anchors.csv"),
        str(TRB / "J_gene_CDR3_anchors.csv"),
    )
    om = olga.GenerativeModelVDJ()
    om.load_and_process_igor_model(str(TRB / "model_marginals.txt"))

    v_names = [x[0] for x in g.genV]
    d_names = [x[0] for x in g.genD]
    j_names = [x[0] for x in g.genJ]
    m = from_olga(TRB, locus="TRB")

    # PV
    pv = _lookup(m.tables["v_choice"], ["v_allele"])
    assert np.allclose([pv[n] for n in v_names], om.PV)

    # PDJ[d, j] = P(D|J) * P(J)
    pj = _lookup(m.tables["j_choice"], ["j_allele"])
    pdj_cond = _lookup(m.tables["d_gene"], ["j_allele", "d_allele"])
    pdj_rec = np.array([[pdj_cond[(j, d)] * pj[j] for j in j_names] for d in d_names])
    assert np.allclose(pdj_rec, om.PDJ)

    # PdelV_given_V[idx, v], idx = ndel + max_pal
    maxp = g.max_delV_palindrome
    dv = _lookup(m.tables["v_3_del"], ["v_allele", "ndel"])
    pdelv_rec = np.array([[dv[(v, idx - maxp)] for v in v_names] for idx in range(om.PdelV_given_V.shape[0])])
    assert np.allclose(pdelv_rec, om.PdelV_given_V)

    # PdelDldelDr_given_D[idx5, idx3, d]
    dd = _lookup(m.tables["d_del"], ["d_allele", "ndel5", "ndel3"])
    n5, n3, nd = om.PdelDldelDr_given_D.shape
    pdd_rec = np.array(
        [[[dd[(d_names[k], i5 - g.max_delDl_palindrome, i3 - g.max_delDr_palindrome)]
           for k in range(nd)] for i3 in range(n3)] for i5 in range(n5)]
    )
    assert np.allclose(pdd_rec, om.PdelDldelDr_given_D)

    # Rvd[next, prev] and PinsVD
    rvd = _lookup(m.tables["vd_dinucl"], ["from_nt", "to_nt"])
    rvd_rec = np.array([[rvd[(prev, nxt)] for prev in range(4)] for nxt in range(4)])
    assert np.allclose(rvd_rec, om.Rvd)
    ins = _lookup(m.tables["vd_ins"], ["length"])
    assert np.allclose([ins[k] for k in range(len(om.PinsVD))], om.PinsVD)

    # cut segments (palindrome-extended germline used by Pgen)
    cutv = dict(zip(m.genomic["genes_v"]["v_allele"], m.genomic["genes_v"]["cut_segment"]))
    assert [cutv[n] for n in v_names] == list(g.cutV_genomic_CDR3_segs)


def test_tra_vj_lossless_and_shapes():
    m = from_olga(TRA, locus="TRA")
    assert m.chain_type == "VJ"
    assert "d_gene" not in m.tables and "genes_d" not in m.genomic
    assert m.tables["v_choice"].height == 103
    assert m.tables["j_choice"].height == 103 * 68  # P(J|V) joint

    g = olga.GenomicDataVJ()
    g.load_igor_genomic_data(
        str(TRA / "model_params.txt"),
        str(TRA / "V_gene_CDR3_anchors.csv"),
        str(TRA / "J_gene_CDR3_anchors.csv"),
    )
    om = olga.GenerativeModelVJ()
    om.load_and_process_igor_model(str(TRA / "model_marginals.txt"))
    v_names = [x[0] for x in g.genV]
    j_names = [x[0] for x in g.genJ]

    pv = _lookup(m.tables["v_choice"], ["v_allele"])
    pjv = _lookup(m.tables["j_choice"], ["v_allele", "j_allele"])
    pvj_rec = np.array([[pjv[(v, j)] * pv[v] for j in j_names] for v in v_names])
    assert np.allclose(pvj_rec, om.PVJ)


@pytest.mark.parametrize("locus,sub", [("TRB", "human_T_beta"), ("TRA", "human_T_alpha")])
def test_parquet_roundtrip(locus, sub, tmp_path):
    m = from_olga(OLGA_MODELS / sub, locus=locus)
    m.save(tmp_path / "model")
    r = Model.load(tmp_path / "model")
    assert r.manifest.to_json() == m.manifest.to_json()
    assert set(r.tables) == set(m.tables)
    for name in m.tables:
        assert r.tables[name].equals(m.tables[name]), f"table {name} changed on round-trip"
    for name in m.genomic:
        assert r.genomic[name].equals(m.genomic[name])
    r.validate()
