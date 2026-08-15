"""Paired-chain generation probability (vdjtools.sc.paired_pgen)."""
import polars as pl

from vdjtools import sc
from vdjtools.model import load_bundled
from vdjtools.model.generate import generate


def _paired_frame(n=8, seed=0):
    """Build a paired α(TRA)/β(TRB) frame from generated (allele-level) clonotypes."""
    a = generate(load_bundled("TRA", "olga"), n, seed=seed, productive_only=True)
    b = generate(load_bundled("TRB", "olga"), n, seed=seed + 1, productive_only=True)
    return pl.DataFrame({
        "cell_id": [f"c{i}" for i in range(n)],
        "alpha_v_call": a["v_call"], "alpha_j_call": a["j_call"], "alpha_junction_aa": a["junction_aa"],
        "beta_v_call": b["v_call"], "beta_j_call": b["j_call"], "beta_junction_aa": b["junction_aa"],
    })


def test_paired_pgen_product_and_columns():
    df = sc.paired_pgen(_paired_frame())
    for c in ("pgen_alpha", "pgen_beta", "pgen_paired"):
        assert c in df.columns
    for r in df.iter_rows(named=True):
        assert r["pgen_alpha"] > 0 and r["pgen_beta"] > 0
        assert abs(r["pgen_paired"] - r["pgen_alpha"] * r["pgen_beta"]) < 1e-30


def test_paired_pgen_conditioning_reduces_pgen():
    """Conditioning on the specific V/J allele gives a smaller Pgen than marginalising."""
    df = _paired_frame()
    cond = sc.paired_pgen(df, condition_vj=True)
    marg = sc.paired_pgen(df, condition_vj=False)
    # marginal (over all V/J) is >= the allele-conditioned value for every cell.
    assert (marg["pgen_beta"] >= cond["pgen_beta"] - 1e-30).all()
    assert (cond["pgen_beta"] < marg["pgen_beta"]).any()


def test_paired_pgen_missing_chain_is_null():
    df = _paired_frame(4)
    df = df.with_columns(
        pl.when(pl.col("cell_id") == "c0").then(None).otherwise(pl.col("alpha_junction_aa"))
        .alias("alpha_junction_aa")
    )
    out = sc.paired_pgen(df)
    r0 = out.filter(pl.col("cell_id") == "c0").row(0, named=True)
    assert r0["pgen_alpha"] is None and r0["pgen_paired"] is None
    assert r0["pgen_beta"] is not None  # β chain still present


def test_paired_pgen_null_v_call_yields_null_chain():
    """A chain whose V-call column is entirely null infers no locus → that chain's Pgen is null."""
    df = _paired_frame(4).with_columns(pl.lit(None, dtype=pl.Utf8).alias("alpha_v_call"))
    out = sc.paired_pgen(df)
    assert out["pgen_alpha"].to_list() == [None] * 4
    assert out["pgen_paired"].to_list() == [None] * 4
    assert out["pgen_beta"].null_count() == 0            # β still scored


def test_paired_pgen_unscoreable_junction_is_null():
    """A junction the native model cannot score yields null, not a crash."""
    df = _paired_frame(4).with_columns(pl.lit(1, dtype=pl.Int64).alias("beta_junction_aa"))
    out = sc.paired_pgen(df)
    assert out["pgen_beta"].to_list() == [None] * 4
    assert out["pgen_paired"].to_list() == [None] * 4


def _strip_alleles(df):
    """Gene-level V/J calls, as CellRanger reports them (TRBV10-3, not TRBV10-3*01)."""
    return df.with_columns(
        pl.col(c).str.split("*").list.first()
        for c in ("alpha_v_call", "alpha_j_call", "beta_v_call", "beta_j_call")
    )


def test_gene_level_calls_are_scored_not_silently_nulled():
    """Regression: CellRanger emits GENE-level calls, and every row used to score null.

    `native.pgen_aa` deliberately raises on a gene name (marginalising silently once gave a
    Pgen 2.38x too high), and paired_pgen swallowed that -- so on the single most common
    real input, all three columns came back null with no error. The gene is now resolved to
    a representative allele before scoring.
    """
    out = sc.paired_pgen(_strip_alleles(_paired_frame()))
    assert out["pgen_paired"].null_count() == 0
    assert all(p > 0 for p in out["pgen_paired"])


def test_resolve_genes_false_refuses_gene_level_calls():
    """Opting out must give nulls, never a silently marginalised (larger) Pgen."""
    out = sc.paired_pgen(_strip_alleles(_paired_frame()), resolve_genes=False)
    assert out["pgen_paired"].null_count() == out.height


def test_an_all_null_locus_warns_rather_than_shipping_a_silent_column():
    import pytest

    bogus = _paired_frame(4).with_columns(pl.lit("TRBV9000").alias("beta_v_call"))
    with pytest.warns(UserWarning, match="every TRB chain scored null"):
        sc.paired_pgen(bogus, resolve_genes=False)
