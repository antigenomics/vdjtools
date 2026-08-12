"""Table export/import, model-directory formats, the training log, and allele-library extension."""
from __future__ import annotations

import json

import polars as pl
import pytest

from vdjtools.model import load_bundled, load_germline
from vdjtools.model.analyze import gene_marginal
from vdjtools.model.check import check_model
from vdjtools.model.generate import generate
from vdjtools.model.infer import extend_alleles, infer_frame, infer_native, training_frame
from vdjtools.model.io import load_model, marginals_frame, save_model, set_marginals


def _gene_usage(model, seg):
    out = {}
    for allele, p in gene_marginal(model, seg).items():
        out[allele.split("*")[0]] = out.get(allele.split("*")[0], 0.0) + p
    return out


# --- flat marginal export / import -------------------------------------------------------------

def test_marginals_frame_covers_every_event(toy_model):
    flat = marginals_frame(toy_model)
    assert set(flat["event"]) == set(toy_model.manifest.events)
    assert flat.height == sum(t.height for t in toy_model.tables.values())
    assert flat.columns[:3] == ["event", "kind", "given"]
    assert flat.columns[-1] == "p"


def test_marginals_round_trip_exactly(toy_model):
    rebuilt = set_marginals(toy_model, marginals_frame(toy_model))
    for name, table in toy_model.tables.items():
        assert rebuilt.tables[name].equals(table), name


def test_marginals_round_trip_through_tsv(tmp_path, toy_model):
    """The practical claim: edit the probabilities in any tool, read them back as a model."""
    path = tmp_path / "marginals.tsv"
    marginals_frame(toy_model).write_csv(path, separator="\t")
    reread = pl.read_csv(path, separator="\t")
    rebuilt = set_marginals(toy_model, reread)
    for name, table in toy_model.tables.items():
        assert rebuilt.tables[name].schema == table.schema, name
        assert rebuilt.tables[name].equals(table), name


def test_hand_edited_marginals_take_effect(toy_model):
    flat = marginals_frame(toy_model).with_columns(
        p=pl.when((pl.col("event") == "v_choice") & (pl.col("v_allele") == "TOYV1*01"))
        .then(0.5).when(pl.col("event") == "v_choice").then(0.25).otherwise(pl.col("p")))
    edited = set_marginals(toy_model, flat)
    assert dict(zip(edited.tables["v_choice"]["v_allele"],
                    edited.tables["v_choice"]["p"]))["TOYV1*01"] == pytest.approx(0.5)


def test_set_marginals_rejects_a_broken_frame(toy_model):
    flat = marginals_frame(toy_model)
    with pytest.raises(ValueError, match="no rows for event"):
        set_marginals(toy_model, flat.filter(pl.col("event") != "vj_ins"))
    with pytest.raises(ValueError, match="not in the manifest"):
        set_marginals(toy_model, flat.with_columns(
            event=pl.when(pl.col("event") == "vj_ins").then(pl.lit("nonsense"))
            .otherwise(pl.col("event"))))
    with pytest.raises(ValueError, match="at least an 'event'"):
        set_marginals(toy_model, flat.drop("event"))


def test_unnormalized_edit_is_caught(toy_model):
    flat = marginals_frame(toy_model).with_columns(
        p=pl.when(pl.col("event") == "v_choice").then(pl.col("p") * 2).otherwise(pl.col("p")))
    with pytest.raises(ValueError, match="sum to neither"):
        set_marginals(toy_model, flat)


# --- model directory formats -------------------------------------------------------------------

@pytest.mark.parametrize("fmt", ["parquet", "tsv", "csv"])
def test_model_directory_round_trip(tmp_path, toy_model, fmt):
    path = tmp_path / fmt
    save_model(toy_model, path, fmt=fmt)
    reloaded = load_model(path, validate=True)
    for name, table in toy_model.tables.items():
        assert reloaded.tables[name].equals(table), name
    assert reloaded.manifest.locus == toy_model.locus
    assert reloaded.chain_type == toy_model.chain_type


def test_unknown_format_is_rejected(tmp_path, toy_model):
    with pytest.raises(ValueError, match="fmt must be"):
        save_model(toy_model, tmp_path / "x", fmt="xlsx")


def test_load_model_does_not_validate_by_default(tmp_path, toy_model):
    """A broken model must still load, or check_model could never diagnose it."""
    broken = type(toy_model)(
        manifest=toy_model.manifest,
        tables={**toy_model.tables,
                "vj_ins": toy_model.tables["vj_ins"].with_columns(p=pl.col("p") * 3)},
        genomic=toy_model.genomic)
    save_model(broken, tmp_path / "broken")
    load_model(tmp_path / "broken")                      # no raise
    with pytest.raises(ValueError):
        load_model(tmp_path / "broken", validate=True)


# --- training log ------------------------------------------------------------------------------

def test_bundled_models_have_no_training_log_and_still_load():
    m = load_bundled("TRG", "olga")
    assert m.training is None
    assert training_frame(m).is_empty()


def test_em_populates_and_persists_the_training_log(tmp_path, toy_model):
    seqs = generate(toy_model, 120, seed=1)["junction_nt"].to_list()
    fitted, report = infer_native(toy_model, seqs, max_iter=2)
    assert fitted.training is not None
    runs = fitted.training["runs"]
    assert len(runs) == 1
    assert runs[0]["n_sequences"] == 120 and runs[0]["native"] is True
    assert runs[0]["max_iter"] == 2 and runs[0]["finished_at"]

    save_model(fitted, tmp_path / "fitted")
    assert json.loads((tmp_path / "fitted" / "training.json").read_text())["runs"]
    reloaded = load_model(tmp_path / "fitted")
    assert reloaded.training == fitted.training
    log = training_frame(reloaded)
    assert log.columns == ["run", "iter", "loglik", "n_scoreable", "rel_change"]
    assert log.height == len(report.loglik)


def test_a_warm_start_appends_a_second_run(toy_model):
    seqs = generate(toy_model, 80, seed=2)["junction_nt"].to_list()
    once, _ = infer_native(toy_model, seqs, max_iter=2)
    twice, _ = infer_native(once, seqs, max_iter=2, init="template")
    assert len(twice.training["runs"]) == 2
    assert set(training_frame(twice)["run"]) == {0, 1}


def test_report_to_frame_and_dict_round_trip(toy_model):
    seqs = generate(toy_model, 60, seed=3)["junction_nt"].to_list()
    _, report = infer_native(toy_model, seqs, max_iter=2)
    assert report.to_frame().height == len(report.loglik)
    from vdjtools.model.infer import InferenceReport

    assert InferenceReport.from_dict(report.to_dict()).loglik == report.loglik
    # An older log missing newer fields must still load.
    assert InferenceReport.from_dict({"loglik": [1.0], "unknown_future_field": 7}).loglik == [1.0]


def test_infer_frame_uses_the_calls(toy_model):
    gen = generate(toy_model, 150, seed=4)
    clones = gen.rename({"junction_nt": "junction"}).select(["junction", "v_call", "j_call"])
    fitted, report = infer_frame(toy_model, clones, max_iter=2)
    assert report.n_iter == 2
    assert fitted.training["runs"][0]["n_sequences"] == clones.height
    fitted.validate()


def test_infer_frame_needs_a_junction_column(toy_model):
    with pytest.raises(ValueError, match="no nucleotide junction column"):
        infer_frame(toy_model, pl.DataFrame({"nope": ["A"]}))


def test_em_log_likelihood_is_monotone(toy_model):
    """EM must never go backwards on its own objective — the property `rel_change` assumes.

    The convergence test is a *relative log-likelihood improvement*, which is only a sound stopping
    rule if the sequence is monotone; a dip would let it stop early on a spurious small change.
    Checked step by step, not just first-vs-last.
    """
    seqs = generate(toy_model, 400, seed=5)["junction_nt"].to_list()
    _, report = infer_native(toy_model, seqs, max_iter=6, tol=0.0)
    assert len(report.loglik) == 6
    for prev, cur in zip(report.loglik, report.loglik[1:]):
        assert cur >= prev - 1e-12, f"log-likelihood decreased: {report.loglik}"


def test_converged_flag_and_iteration_cap(toy_model):
    """`converged` must mean "the tolerance was met", not "the cap was hit"."""
    seqs = generate(toy_model, 400, seed=5)["junction_nt"].to_list()
    _, loose = infer_native(toy_model, seqs, max_iter=20, tol=1e-2)
    assert loose.converged and loose.n_iter < 20
    _, capped = infer_native(toy_model, seqs, max_iter=2, tol=0.0)
    assert not capped.converged and capped.n_iter == 2


# --- extending the allele library ---------------------------------------------------------------

def test_extend_adds_alleles_and_a_new_gene(toy_model, toy_germline_extended):
    extended = extend_alleles(toy_model, toy_germline_extended)
    alleles = set(extended.genomic["genes_v"]["v_allele"])
    assert {"TOYV1*03", "TOYV3*01"} <= alleles
    assert alleles >= set(toy_model.genomic["genes_v"]["v_allele"])
    extended.validate()


def test_new_alleles_of_known_genes_leave_gene_usage_untouched(toy_model, toy_germline_more_alleles):
    """The core invariant: alleles of a gene are alternatives, so a gene's total must not grow.

    Without this, extending human TRB from one to ~3 alleles per gene moved gene-level V usage by
    up to 6 percentage points and silently reweighted every Pgen through those genes.
    """
    extended = extend_alleles(toy_model, toy_germline_more_alleles)
    assert "TOYV1*03" in set(extended.genomic["genes_v"]["v_allele"])
    before, after = _gene_usage(toy_model, "v"), _gene_usage(extended, "v")
    assert set(after) == set(before)
    for gene, p in before.items():
        assert after[gene] == pytest.approx(p, abs=1e-12)


def test_a_brand_new_gene_gets_a_floor_and_rescales_the_rest_proportionally(
        toy_model, toy_germline_extended):
    """A gene with no evidence takes a floor share; everything else keeps its RELATIVE usage."""
    extended = extend_alleles(toy_model, toy_germline_extended)
    before, after = _gene_usage(toy_model, "v"), _gene_usage(extended, "v")
    assert after["TOYV3"] > 0
    assert after["TOYV3"] < min(before.values())        # a floor, not an estimate
    assert sum(after.values()) == pytest.approx(1.0)
    # Pre-existing genes are only diluted by the floor, never reordered or reweighted amongst
    # themselves: every pairwise ratio survives exactly.
    ratio = {g: after[g] / before[g] for g in before}
    assert max(ratio.values()) == pytest.approx(min(ratio.values()))


def test_extend_is_idempotent(toy_model, toy_germline_extended):
    once = extend_alleles(toy_model, toy_germline_extended)
    twice = extend_alleles(once, toy_germline_extended)
    assert twice.genomic["genes_v"].height == once.genomic["genes_v"].height
    for name, table in once.tables.items():
        assert twice.tables[name].height == table.height, name


def test_extend_never_changes_an_existing_allele_germline(toy_model, toy_germline_extended):
    extended = extend_alleles(toy_model, toy_germline_extended)
    old = dict(zip(toy_model.genomic["genes_v"]["v_allele"],
                   toy_model.genomic["genes_v"]["cut_segment"]))
    new = dict(zip(extended.genomic["genes_v"]["v_allele"],
                   extended.genomic["genes_v"]["cut_segment"]))
    for allele, seq in old.items():
        assert new[allele] == seq


def test_extended_model_stays_clean_and_generates(toy_model, toy_germline_extended):
    extended = extend_alleles(toy_model, toy_germline_extended)
    assert check_model(extended, germline="none").filter(
        pl.col("severity") == "error").is_empty()
    forced = type(extended)(
        manifest=extended.manifest,
        tables={**extended.tables, "v_choice": extended.tables["v_choice"].with_columns(
            p=pl.when(pl.col("v_allele") == "TOYV3*01").then(1.0).otherwise(0.0))},
        genomic=extended.genomic)
    gen = generate(forced, 20, seed=6)
    assert set(gen["v_call"]) == {"TOYV3*01"}          # the new gene is really usable


# --- the shipped real-read examples -------------------------------------------------------------

@pytest.mark.parametrize("chain,min_rows", [("TRB", 50_000), ("TRA", 20_000)])
def test_prepared_examples_load(chain, min_rows):
    """The arda-mapped TRA/TRB examples that ship in the tree — no network, no arda, no mmseqs2."""
    from vdjtools.model.data import load_prepared

    clones = load_prepared("human", chain)
    assert clones.height >= min_rows
    assert clones.columns == ["junction", "v_call", "j_call", "d_call", "d2_call", "count"]
    assert clones["junction"].null_count() == 0
    assert clones["v_call"].str.starts_with(f"{chain[:2]}").all()
    assert clones["j_call"].str.contains(f"{chain}J").all()
    assert clones["count"].min() >= 1


def test_prepared_example_trains_a_model():
    """The end the examples exist for: real junctions, straight into EM."""
    from vdjtools.model.data import load_prepared

    clones = load_prepared("human", "TRA").head(2000)
    template = load_bundled("TRA", "arda", collapse=False)
    fitted, report = infer_frame(template, clones, max_iter=2)
    assert report.n_iter == 2
    assert report.loglik[-1] >= report.loglik[0]
    fitted.validate()


def _with_ambiguous(toy_model, n=60):
    gen = generate(toy_model, n, seed=8).rename({"junction_nt": "junction"})
    return pl.concat([
        gen.select(["junction", "v_call", "j_call"]),
        pl.DataFrame({"junction": ["TGTNCCAGC"], "v_call": ["TOYV1*01"], "j_call": ["TOYJ1*01"]}),
    ])


def test_ambiguous_bases_are_substituted_by_default(toy_model):
    """Real reads carry 'N' and the encoder knows only ACGT — substitute, keep the clonotype."""
    with_n = _with_ambiguous(toy_model)
    with pytest.warns(UserWarning, match="substituting"):
        fitted, report = infer_frame(toy_model, with_n, max_iter=1)
    assert report.n_scoreable[0] == with_n.height          # nothing was thrown away
    fitted.validate()


def test_ambiguous_bases_can_be_dropped_instead(toy_model):
    with_n = _with_ambiguous(toy_model)
    with pytest.warns(UserWarning, match="dropped"):
        fitted, report = infer_frame(toy_model, with_n, max_iter=1, ambiguous=None)
    assert report.n_scoreable[0] == with_n.height - 1
    fitted.validate()


def test_sanitize_junctions_substitutes_every_iupac_code():
    from vdjtools.model.infer import sanitize_junctions

    df = pl.DataFrame({"junction": ["acgTN", "ACGRY", "ACGT"]})
    with pytest.warns(UserWarning, match="substituting"):
        out = sanitize_junctions(df, "junction")
    assert out["junction"].to_list() == ["ACGTA", "ACGAA", "ACGT"]
    with pytest.raises(ValueError, match="ambiguous must be"):
        sanitize_junctions(df, "junction", ambiguous="N")


def test_prepared_fasta_round_trip(tmp_path):
    """write_prepared -> load_prepared is exact for the fields a FASTA header carries."""
    from vdjtools.model.data import load_prepared, write_prepared

    clones = pl.DataFrame({
        "junction": ["TGTGCCAGC", "TGTGCTTCC"],
        "v_call": ["TRBV1*01", "TRBV2*01,TRBV3*01"],
        "j_call": ["TRBJ1*01", "TRBJ2*01"],
        "d_call": ["TRBD1*01", None],
        "d2_call": [None, None],
        "count": [7, 1],
    })
    path = write_prepared(clones, tmp_path / "x.fa.gz")
    back = load_prepared(path=path)
    assert back.to_dicts() == clones.to_dicts()


def test_extend_a_real_model_with_the_full_arda_library():
    """The realistic case: a gene-collapsed model meeting the whole IMGT allele set."""
    m = load_bundled("TRG", "learned")
    extended = extend_alleles(m, load_germline("TRG", "human"))
    assert extended.genomic["genes_v"].height > m.genomic["genes_v"].height
    extended.validate()
    before, after = _gene_usage(m, "v"), _gene_usage(extended, "v")
    for gene, p in before.items():
        assert after[gene] == pytest.approx(p, abs=5e-3)
