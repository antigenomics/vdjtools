"""Tests for model.analyze — entropy, mutual information, and the Bayes-net DOT."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import pytest

from vdjtools.model import analyze, from_olga

OLGA = Path("/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models")
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA.exists(), reason="OLGA default models not available")


@pytest.fixture(scope="module")
def trb():
    return from_olga(OLGA / "human_T_beta", locus="TRB")


@pytest.fixture(scope="module")
def tra():
    return from_olga(OLGA / "human_T_alpha", locus="TRA")


def test_gene_marginals_normalized(trb):
    for seg in ("v", "j", "d"):
        assert abs(sum(analyze.gene_marginal(trb, seg).values()) - 1.0) < 1e-6


def test_entropy_nonnegative_and_conditioning_reduces(trb):
    ent = analyze.entropy_table(trb)
    assert ent.height == len(trb.manifest.events)
    assert (ent["H_bits"] >= -1e-9).all()
    # conditioning never increases entropy: H(X|parents) <= H(X)
    assert (ent["H_cond_bits"] <= ent["H_bits"] + 1e-9).all()
    # single-D OLGA model: n_d is degenerate → zero entropy
    assert ent.filter(pl.col("event") == "n_d")["H_bits"][0] == pytest.approx(0.0, abs=1e-9)


def test_mi_nonnegative_and_vj_independent_for_vdj(trb):
    mi = analyze.mutual_information(trb)
    assert (mi["mi_bits"] >= -1e-9).all()
    vj = mi.filter((pl.col("a") == "v_choice") & (pl.col("b") == "j_choice"))
    assert vj.height == 1 and vj["mi_bits"][0] == pytest.approx(0.0, abs=1e-9)


def test_mi_equals_entropy_drop(trb):
    """I(child;parent) reported by mutual_information == H(child) - H(child|parent)."""
    ent = {r["event"]: (r["H_bits"], r["H_cond_bits"]) for r in analyze.entropy_table(trb).to_dicts()}
    for r in analyze.mutual_information(trb).to_dicts():
        if r["a"] in ent and r["b"] in trb.manifest.events[r["a"]].given:
            h, hc = ent[r["a"]]
            assert r["mi_bits"] == pytest.approx(h - hc, abs=1e-3)


def test_vj_model_has_vj_coupling(tra):
    """A VJ locus (TRA) encodes P(J|V) → I(V;J) > 0 (unlike an independent-root VDJ)."""
    mi = analyze.mutual_information(tra)
    edge = mi.filter((pl.col("a") == "j_choice") & (pl.col("b") == "v_choice"))
    assert edge.height == 1 and edge["mi_bits"][0] > 0.0


def test_dot_contains_nodes_and_renders(trb, tmp_path):
    dot = analyze.bayes_net_dot(trb)
    for ev in trb.manifest.events:
        assert f'"{ev}"' in dot
    assert "I=" in dot and "H=" in dot
    import shutil
    if shutil.which("dot"):
        out = analyze.render_bayes_net(trb, tmp_path / "trb", fmt="png")
        assert out.exists() and out.stat().st_size > 0


def test_compare_entropy_wide(trb, tra):
    wide = analyze.compare_entropy({"TRB": trb, "TRA": tra})
    assert "TRB" in wide.columns and "TRA" in wide.columns
    # v_choice present in both; VDJ-only events (d_gene) null for the VJ model
    row = wide.filter(pl.col("event") == "d_gene")
    assert row["TRA"][0] is None and row["TRB"][0] is not None


def test_multiparent_event_unsupported(trb):
    """analyze handles single-parent factorizations only; a >=2-parent event must fail loudly."""
    from vdjtools.model import Event, EventKind, Manifest, Model
    ev = dict(trb.manifest.events)
    ev["d_gene"] = Event("d_gene", EventKind.GENE_CHOICE, ("j_choice", "v_choice"))
    man = Manifest(locus="TRB", organism="human", chain_type="VDJ", events=ev,
                   palindrome_max=trb.manifest.palindrome_max, source="x")
    m2 = Model(manifest=man, tables=trb.tables, genomic=trb.genomic)
    with pytest.raises(NotImplementedError, match="multi-parent"):
        analyze.mutual_information(m2)


def test_analyze_on_dd_model(trb):
    """Diagnostics on a tandem-D model: MI>=0, the d2_gene edge identity, the within-D2 coupling row."""
    from vdjtools.model.dd import to_dd
    m = to_dd(trb, p_nd2=0.05)
    mi = analyze.mutual_information(m)
    assert (mi["mi_bits"] >= -1e-9).all()
    ent = {r["event"]: (r["H_bits"], r["H_cond_bits"]) for r in analyze.entropy_table(m).to_dicts()}
    edge = mi.filter((pl.col("a") == "d2_gene") & (pl.col("b") == "d_gene"))
    h, hc = ent["d2_gene"]
    assert edge["mi_bits"][0] == pytest.approx(h - hc, abs=1e-3)
    assert mi.filter(pl.col("a") == "delD2_5").height == 1  # second-D within-coupling reported
    assert abs(sum(analyze.gene_marginal(m, "d2").values()) - 1.0) < 1e-6


def test_compare_entropy_dd_null(trb):
    """compare_entropy: D-D-only events are null for the single-D model, present for the D-D one."""
    from vdjtools.model.dd import to_dd
    wide = analyze.compare_entropy({"DD": to_dd(trb, p_nd2=0.05), "SD": trb})
    row = wide.filter(pl.col("event") == "d2_gene")
    assert row["SD"][0] is None and row["DD"][0] is not None
