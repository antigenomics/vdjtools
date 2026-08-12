"""Scoring sequences under a model: Pgen frames, free parameters, likelihood/BIC, diversity.

The free-parameter count is asserted against a number worked out by hand on the toy locus, because
that is the only way to catch the failure mode that matters: counting table *rows* instead of
occupied cells, which would inflate ``k`` several-fold and quietly wreck every BIC comparison.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest

from vdjtools.model import native
from vdjtools.model.generate import generate
from vdjtools.model.score import (
    compare_pgen,
    diversity,
    free_params,
    model_fit,
    pgen_frame,
    pgen_spectrum,
    pgen_summary,
)

# Worked out by hand from the toy VJ scaffold (see conftest), per normalization group:
#   v_choice   3 alleles, one group                     ->  3 - 1 =  2
#   j_choice   P(J|V): 3 V groups x 2 J alleles         -> 3 * (2 - 1) =  3
#   v_3_del    3 alleles x ndel -4..13 (18 bins)        -> 3 * (18 - 1) = 51
#   j_5_del    2 alleles x ndel -4..22 (27 bins)        -> 2 * (27 - 1) = 52
#   vj_ins     lengths 0..3, one group                  ->  4 - 1 =  3
#   vj_dinucl  column-stochastic 4x4, 4 groups          -> 4 * (4 - 1) = 12
TOY_FREE_PARAMS = 2 + 3 + 51 + 52 + 3 + 12


def test_free_params_matches_the_hand_count(toy_model):
    assert free_params(toy_model) == TOY_FREE_PARAMS


def test_free_params_by_event_sums_to_the_total(toy_model):
    per_event = free_params(toy_model, by_event=True)
    assert per_event["k"].sum() == TOY_FREE_PARAMS
    assert set(per_event["event"]) == set(toy_model.manifest.events)


def test_structural_zeros_are_not_parameters(toy_model):
    """Zeroing cells must lower k: an empty cell parameterizes nothing."""
    before = free_params(toy_model)
    tables = dict(toy_model.tables)
    # Knock out three deletion bins for one allele and renormalize that allele's group.
    t = tables["v_3_del"].with_columns(
        p=pl.when((pl.col("v_allele") == "TOYV1*01") & (pl.col("ndel") > 10))
        .then(0.0).otherwise(pl.col("p")))
    t = t.with_columns(p=pl.col("p") / pl.col("p").sum().over("v_allele"))
    tables["v_3_del"] = t
    zeroed = type(toy_model)(manifest=toy_model.manifest, tables=tables, genomic=toy_model.genomic)
    assert free_params(zeroed) == before - 3


def test_unreachable_conditionals_are_excluded(toy_model):
    """A group whose parent has zero marginal is not estimable and must not be counted."""
    tables = dict(toy_model.tables)
    tables["v_choice"] = tables["v_choice"].with_columns(
        p=pl.when(pl.col("v_allele") == "TOYV2*01").then(0.0).otherwise(pl.col("p")))
    tables["v_choice"] = tables["v_choice"].with_columns(p=pl.col("p") / pl.col("p").sum())
    dead = type(toy_model)(manifest=toy_model.manifest, tables=tables, genomic=toy_model.genomic)
    # TOYV2*01 loses its own v_choice cell (-1), its 18-bin deletion group (-17), and the
    # P(J|V=TOYV2*01) group (-1).
    assert free_params(dead, reachable_only=True) == TOY_FREE_PARAMS - 1 - 17 - 1
    assert free_params(dead, reachable_only=False) == TOY_FREE_PARAMS - 1


def test_pgen_frame_matches_per_sequence_pgen(toy_model):
    seqs = generate(toy_model, 20, seed=3)["junction_nt"].to_list()
    frame = pgen_frame(toy_model, seqs, use_calls=False)
    expected = [native.pgen_nt(toy_model, s) for s in seqs]
    assert frame["pgen"].to_list() == expected
    assert frame["kind"].to_list() == ["nt"] * len(seqs)


def test_threaded_and_serial_nt_pgen_agree(small_model):
    """The thread pool must be an optimization only — bitwise-identical, thread-count-invariant."""
    seqs = generate(small_model, 80, seed=5)["junction_nt"].to_list()
    serial = pgen_frame(small_model, seqs, use_calls=False, threads=1)["pgen"].to_numpy()
    threaded = pgen_frame(small_model, seqs, use_calls=False, threads=4)["pgen"].to_numpy()
    assert np.array_equal(serial, threaded)


def test_unscoreable_sequence_is_null_not_negative_infinity(toy_model):
    """A sequence the model cannot generate must be counted, never turned into -inf."""
    frame = pgen_frame(toy_model, ["TGTGCCAGCAACTATGGCTATACCTTT", "TTTTTTTTTTTTTTT"],
                       use_calls=False)
    assert frame["scoreable"].to_list() == [True, False]
    assert frame["pgen"][1] == 0.0
    assert frame["log_pgen"][1] is None


def test_model_fit_reports_scoreable_share_and_finite_loglik(toy_model):
    seqs = generate(toy_model, 40, seed=7)["junction_nt"].to_list() + ["TTTTTTTTTTTTTTT"]
    fit = model_fit(toy_model, seqs, use_calls=False).to_dicts()[0]
    assert fit["n"] == 41
    assert fit["n_scoreable"] == 40
    assert fit["frac_scoreable"] == pytest.approx(40 / 41)
    assert math.isfinite(fit["loglik_sum"])
    assert fit["kind"] == "nt"


def test_bic_identity(toy_model):
    seqs = generate(toy_model, 30, seed=8)["junction_nt"].to_list()
    fit = model_fit(toy_model, seqs, use_calls=False).to_dicts()[0]
    assert fit["k"] == TOY_FREE_PARAMS
    assert fit["bic"] == pytest.approx(fit["k"] * math.log(fit["n_scoreable"]) - 2 * fit["loglik_sum"])
    assert fit["aic"] == pytest.approx(2 * fit["k"] - 2 * fit["loglik_sum"])


def test_weights_equal_repeating_the_sequences(toy_model):
    seqs = generate(toy_model, 12, seed=9)["junction_nt"].to_list()
    weights = [3] * len(seqs)
    weighted = model_fit(toy_model, seqs, weights=weights, use_calls=False).to_dicts()[0]
    repeated = model_fit(toy_model, seqs * 3, use_calls=False).to_dicts()[0]
    assert weighted["loglik_sum"] == pytest.approx(repeated["loglik_sum"])
    assert weighted["n"] == repeated["n"]


def test_the_generating_model_explains_its_own_sequences_best(toy_model):
    """A model must fit its own draws better than one whose insertion length is wrong."""
    seqs = generate(toy_model, 300, seed=11)["junction_nt"].to_list()
    tables = dict(toy_model.tables)
    tables["vj_ins"] = tables["vj_ins"].with_columns(p=pl.Series([0.01, 0.02, 0.07, 0.90]))
    wrong = type(toy_model)(manifest=toy_model.manifest, tables=tables, genomic=toy_model.genomic)
    truth = model_fit(toy_model, seqs, use_calls=False)["loglik_sum"][0]
    assert truth > model_fit(wrong, seqs, use_calls=False)["loglik_sum"][0]


def test_gene_level_call_resolves_when_unambiguous(small_model):
    """Bundled models are gene-collapsed, so a bare gene name must resolve rather than raise."""
    gen = generate(small_model, 5, seed=13)
    genes = [c.split("*")[0] for c in gen["v_call"].to_list()]
    frame = pgen_frame(small_model, gen["junction_nt"].to_list(), v=genes,
                       j=gen["j_call"].to_list())
    assert frame["v_call"].to_list() == gen["v_call"].to_list()


def test_unknown_call_raises_unless_marginalizing(small_model):
    seqs = generate(small_model, 3, seed=15)["junction_nt"].to_list()
    with pytest.raises(KeyError, match="not in this model"):
        pgen_frame(small_model, seqs, v=["NOSUCHV*01"] * 3)
    frame = pgen_frame(small_model, seqs, v=["NOSUCHV*01"] * 3, on_unknown="marginalize")
    assert frame["v_call"].to_list() == [None] * 3


def test_mixed_nt_and_aa_input_is_rejected(toy_model):
    with pytest.raises(ValueError, match="mixed"):
        pgen_frame(toy_model, ["TGTGCCAGC", "CASSLF"], use_calls=False)


def test_compare_pgen_of_a_model_with_itself_is_exact(toy_model):
    seqs = generate(toy_model, 60, seed=17)["junction_nt"].to_list()
    cmp = compare_pgen(toy_model, toy_model, seqs, use_calls=False)
    assert (cmp["delta_log10"].drop_nulls().abs() < 1e-12).all()
    summary = pgen_summary(cmp).to_dicts()[0]
    assert summary["only_a_scoreable"] == 0 and summary["only_b_scoreable"] == 0
    assert summary["spearman_log10"] == pytest.approx(1.0)
    assert summary["ks_stat"] == pytest.approx(0.0)


def test_pgen_summary_counts_one_sided_coverage(toy_model):
    """The headline number: sequences one model can score and the other cannot."""
    seqs = generate(toy_model, 30, seed=19)["junction_nt"].to_list()
    tables = dict(toy_model.tables)
    tables["v_choice"] = tables["v_choice"].with_columns(
        p=pl.when(pl.col("v_allele") == "TOYV1*01").then(1.0).otherwise(0.0))
    narrow = type(toy_model)(manifest=toy_model.manifest, tables=tables, genomic=toy_model.genomic)
    summary = pgen_summary(compare_pgen(toy_model, narrow, seqs, use_calls=False)).to_dicts()[0]
    assert summary["only_a_scoreable"] > 0
    assert summary["only_b_scoreable"] == 0


def test_diversity_reports_both_hill_numbers(toy_model):
    d = diversity(toy_model, n=400, seed=21).to_dicts()[0]
    assert d["scenario_entropy_bits"] > 0
    assert d["sequence_entropy_bits"] > 0
    # Different scenarios can yield the same junction, so the sequence distribution can never
    # carry more information than the scenario it came from.
    assert d["sequence_entropy_bits"] <= d["scenario_entropy_bits"] + 1e-9
    assert d["scenario_diversity"] == pytest.approx(2.0 ** d["scenario_entropy_bits"])
    assert d["diversity_shannon"] == pytest.approx(2.0 ** d["sequence_entropy_bits"])
    assert d["diversity_simpson"] == pytest.approx(1.0 / d["pgen_mean"])
    # Hill numbers decrease in q: the coincidence-based diversity weights common sequences more.
    assert d["diversity_simpson"] <= d["diversity_shannon"]


def test_diversity_is_reproducible(toy_model):
    a = diversity(toy_model, n=200, seed=23)
    b = diversity(toy_model, n=200, seed=23)
    assert a.equals(b)


def test_narrower_model_is_less_diverse(toy_model):
    tables = dict(toy_model.tables)
    tables["v_choice"] = tables["v_choice"].with_columns(
        p=pl.when(pl.col("v_allele") == "TOYV1*01").then(1.0).otherwise(0.0))
    narrow = type(toy_model)(manifest=toy_model.manifest, tables=tables, genomic=toy_model.genomic)
    wide = diversity(toy_model, n=400, seed=25)["sequence_entropy_bits"][0]
    assert diversity(narrow, n=400, seed=25)["sequence_entropy_bits"][0] < wide


def test_pgen_spectrum_is_a_normalized_histogram(toy_model):
    spec = pgen_spectrum(toy_model, n=300, seed=27, bins=8)
    assert spec.height == 8
    assert spec["count"].sum() == 300
    assert spec["frac"].sum() == pytest.approx(1.0)
    assert (spec["bin_left"] < spec["bin_right"]).all()
