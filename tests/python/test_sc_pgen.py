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
