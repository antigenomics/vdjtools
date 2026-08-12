"""Comparing two models: per-event divergence, usage, the comparison graph, and total entropy."""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model import load_bundled
from vdjtools.model.analyze import (
    compare_models,
    compare_net_dot,
    compare_usage,
    entropy_table,
    total_entropy,
)
from vdjtools.model.model import Model


def _rebuild(model, **tables):
    return Model(manifest=model.manifest, tables={**model.tables, **tables},
                 genomic=model.genomic)


def test_a_model_compared_with_itself_is_all_zero(toy_model):
    cmp = compare_models(toy_model, toy_model)
    assert set(cmp["status"]) == {"shared"}
    assert (cmp["tv"] == 0.0).all()
    assert (cmp["jsd_bits"] == 0.0).all()
    assert (cmp["support_only_a"] == 0).all()
    assert (cmp["support_only_b"] == 0).all()


def test_compare_columns(toy_model):
    cmp = compare_models(toy_model, toy_model)
    assert cmp.columns == ["event", "kind", "given", "status", "n_groups", "support_a",
                           "support_b", "support_shared", "support_only_a", "support_only_b",
                           "tv", "tv_max", "jsd_bits"]


def test_total_variation_is_exact_on_a_hand_made_difference(toy_model):
    """v_choice is uniform over 3 alleles; move it to (0.5, 0.25, 0.25) and check TV by hand."""
    moved = pl.DataFrame({"v_allele": ["TOYV1*01", "TOYV1*02", "TOYV2*01"],
                          "p": [0.5, 0.25, 0.25]})
    other = _rebuild(toy_model, v_choice=moved)
    row = compare_models(toy_model, other).filter(pl.col("event") == "v_choice").to_dicts()[0]
    third = 1.0 / 3.0
    expected = 0.5 * (abs(third - 0.5) + abs(third - 0.25) + abs(third - 0.25))
    assert row["tv"] == pytest.approx(expected, abs=1e-6)
    assert 0.0 < row["jsd_bits"] <= 1.0


def test_jsd_is_finite_on_disjoint_support(toy_model):
    """The reason JSD is the headline metric and KL is not reported at all."""
    a = _rebuild(toy_model, v_choice=pl.DataFrame(
        {"v_allele": ["TOYV1*01", "TOYV1*02", "TOYV2*01"], "p": [1.0, 0.0, 0.0]}))
    b = _rebuild(toy_model, v_choice=pl.DataFrame(
        {"v_allele": ["TOYV1*01", "TOYV1*02", "TOYV2*01"], "p": [0.0, 0.0, 1.0]}))
    row = compare_models(a, b).filter(pl.col("event") == "v_choice").to_dicts()[0]
    assert row["tv"] == pytest.approx(1.0)
    assert row["jsd_bits"] == pytest.approx(1.0)      # the maximum, and finite
    assert row["support_shared"] == 0
    assert row["support_only_a"] == 1 and row["support_only_b"] == 1


def test_tv_max_finds_the_one_broken_group(toy_model):
    """A single wrecked allele must be visible in tv_max even when the average stays small."""
    t = toy_model.tables["v_3_del"]
    spiked = t.with_columns(
        p=pl.when((pl.col("v_allele") == "TOYV2*01") & (pl.col("ndel") == 0))
        .then(1.0).otherwise(pl.when(pl.col("v_allele") == "TOYV2*01").then(0.0)
                             .otherwise(pl.col("p"))))
    row = compare_models(toy_model, _rebuild(toy_model, v_3_del=spiked)).filter(
        pl.col("event") == "v_3_del").to_dicts()[0]
    assert row["tv_max"] > row["tv"]


def test_events_missing_from_one_model_are_reported_not_dropped(toy_model, toy_model_vdj):
    cmp = compare_models(toy_model, toy_model_vdj)
    only_b = cmp.filter(pl.col("status") == "only_b")["event"].to_list()
    assert "d_gene" in only_b and "n_d" in only_b
    assert cmp.filter(pl.col("status") == "only_b")["tv"].null_count() == len(only_b)
    assert "vj_ins" in cmp.filter(pl.col("status") == "only_a")["event"].to_list()


def test_a_differently_factorized_event_is_flagged_not_joined(toy_model, toy_model_vdj):
    """j_choice is P(J|V) on a VJ locus and a root P(J) on a VDJ one — not comparable."""
    cmp = compare_models(toy_model, toy_model_vdj)
    row = cmp.filter(pl.col("event") == "j_choice").to_dicts()[0]
    assert row["status"] == "schema_differs"
    assert row["tv"] is None and row["jsd_bits"] is None
    assert row["support_a"] > 0 and row["support_b"] > 0


def test_gene_level_alignment_matches_differing_allele_namespaces(toy_model):
    """Two models whose alleles differ but whose genes agree only line up with by='gene'."""
    renamed_v = toy_model.genomic["genes_v"].with_columns(
        v_allele=pl.col("v_allele").str.replace(r"\*0\d$", "*09"))
    tables = {k: (v.with_columns(v_allele=pl.col("v_allele").str.replace(r"\*0\d$", "*09"))
                  if "v_allele" in v.columns else v)
              for k, v in toy_model.tables.items()}
    other = Model(manifest=toy_model.manifest, tables=tables,
                  genomic={**toy_model.genomic, "genes_v": renamed_v})
    by_allele = compare_models(toy_model, other, by="allele").filter(
        pl.col("event") == "v_choice").to_dicts()[0]
    by_gene = compare_models(toy_model, other, by="gene").filter(
        pl.col("event") == "v_choice").to_dicts()[0]
    assert by_allele["support_shared"] == 0            # nothing lines up allele-to-allele
    assert by_gene["tv"] == pytest.approx(0.0)          # but the genes are identical
    assert by_gene["support_shared"] > 0


def test_compare_usage(toy_model):
    other = _rebuild(toy_model, v_choice=pl.DataFrame(
        {"v_allele": ["TOYV1*01", "TOYV1*02", "TOYV2*01"], "p": [0.5, 0.25, 0.25]}))
    usage = compare_usage(toy_model, other, "v")
    assert usage.columns == ["name", "p_a", "p_b", "log2_ratio"]
    assert set(usage["name"]) == {"TOYV1", "TOYV2"}     # gene level by default
    assert usage["p_a"].sum() == pytest.approx(1.0)
    assert usage["p_b"].sum() == pytest.approx(1.0)
    allele_level = compare_usage(toy_model, other, "v", by="allele")
    assert set(allele_level["name"]) == set(toy_model.tables["v_choice"]["v_allele"])


def test_compare_net_dot_encodes_structure_and_divergence(toy_model, toy_model_vdj):
    dot = compare_net_dot(toy_model, toy_model_vdj, labels=("vj", "vdj"))
    assert dot.startswith("digraph compare {") and dot.rstrip().endswith("}")
    assert '"v_choice"' in dot and '"d_gene"' in dot
    assert "vj only" in dot and "vdj only" in dot     # both one-sided edge styles present
    assert "JSD=" in dot and "dH=" in dot


def test_total_entropy_sums_the_conditional_entropies(toy_model):
    te = total_entropy(toy_model)
    assert te.columns == ["event", "kind", "contribution_bits"]
    assert set(te["event"]) == set(toy_model.manifest.events)
    assert (te["contribution_bits"] >= 0).all()
    assert te["contribution_bits"].sum() > 0


def test_dinucleotide_contribution_scales_with_insertion_length(toy_model):
    """An N-region contributes per-step entropy times its expected length, not per-step alone."""
    base = total_entropy(toy_model)
    step = entropy_table(toy_model).filter(
        pl.col("event") == "vj_dinucl")["H_cond_bits"][0]
    ins = toy_model.tables["vj_ins"]
    mean_len = float((ins["length"].cast(pl.Float64) * ins["p"]).sum())
    got = base.filter(pl.col("event") == "vj_dinucl")["contribution_bits"][0]
    assert got == pytest.approx(step * mean_len, rel=1e-6)


def test_real_models_differ_most_in_v_usage():
    """OLGA is DNA-multiplex, `learned` is 5'RACE: V usage should be the biggest divergence."""
    cmp = compare_models(load_bundled("TRB", "olga"), load_bundled("TRB", "learned"), by="gene")
    shared = cmp.filter(pl.col("status") == "shared").drop_nulls("tv")
    assert shared.height > 5
    assert shared["tv"].max() == shared.filter(pl.col("event") == "v_choice")["tv"][0]
