"""Collapse a V(D)J model to one allele per gene — the default working resolution.

IMGT alleles of one gene (``TRBV12-4*01``, ``TRBV12-4*02``) differ by a handful of SNPs, usually
in framework rather than the CDR3 region, and short-read aligners cannot reliably tell them apart
(arda splits usage across ``*01``/``*03`` as mapping noise). So the default working model sums each
gene's sub-variants into a single ``*01`` representative, and Pgen collapses a clonotype's allele to
that representative too — the allele suffix stops carrying (spurious) information. Pass
``collapse=False`` (or keep the un-collapsed model) when you genuinely want allele resolution, e.g.
the exact-OLGA-Pgen fidelity check.

The collapse is a proper marginalisation, not a truncation:

* a **choice** probability (``P(V=a)``) becomes the plain sum over the gene's alleles — exact;
* a **conditional** (``P(delV | V=a)``, ``P(D | J=ja)``) becomes the usage-weighted average over
  the gene's alleles, ``P(x | gene) = Σ_a [P(a)/P(gene)] · P(x | a)`` — the correct allele-marginal;
* the **germline** of the highest-usage allele is kept as the gene's representative and relabelled
  ``gene*01``. This is the one lossy step: where alleles differ *inside* the CDR3 region the mixture
  cannot be reproduced by a single germline, so collapsed Pgen is approximate there (exact wherever
  the CDR3-region germline is allele-invariant, which is the common case).
"""
from __future__ import annotations

import polars as pl

from .model import Model


def _gene(col: str) -> pl.Expr:
    return pl.col(col).str.split("*").list.first()


def _rep(gene: str) -> str:
    return f"{gene}*01"


def _allele_weights(choice: pl.DataFrame, allele_col: str) -> pl.DataFrame:
    """``(allele, gene, w)`` where ``w = P(allele) / P(gene)`` — the within-gene mixing weight.

    A gene whose alleles are all zero (never used) gets uniform weights so its collapsed
    conditional is defined (harmless: the collapsed choice mass is zero, so it is never sampled).
    """
    c = choice.select(allele_col, "p").with_columns(_gene(allele_col).alias("_g"))
    tot = pl.col("p").sum().over("_g")
    n = pl.len().over("_g")
    return c.with_columns(w=pl.when(tot > 0).then(pl.col("p") / tot).otherwise(1.0 / n)).select(
        pl.col(allele_col).alias("_a"), pl.col("_g"), "w")


def _collapse_choice(choice: pl.DataFrame, allele_col: str) -> pl.DataFrame:
    """Marginal choice at gene resolution, relabelled ``gene*01``."""
    return (choice.with_columns(_gene(allele_col).alias("_g"))
            .group_by("_g").agg(pl.col("p").sum())
            .with_columns(pl.col("_g").map_elements(_rep, return_dtype=pl.Utf8).alias(allele_col))
            .drop("_g").select(allele_col, "p"))


def _collapse_conditional(tbl: pl.DataFrame, allele_col: str, weights: pl.DataFrame,
                          bin_cols: list[str]) -> pl.DataFrame:
    """Usage-weighted average of a per-allele conditional, keyed ``gene*01`` + its bins."""
    w = weights.rename({"_a": allele_col})
    out = (tbl.join(w, on=allele_col, how="inner")
           .with_columns(pw=pl.col("p") * pl.col("w"))
           .group_by(["_g", *bin_cols]).agg(pl.col("pw").sum().alias("p"))
           .with_columns(pl.col("_g").map_elements(_rep, return_dtype=pl.Utf8).alias(allele_col))
           .drop("_g"))
    return out.select(allele_col, *bin_cols, "p")


def _marginal(choice: pl.DataFrame, allele_col: str) -> pl.DataFrame:
    """P over the raw alleles of a root choice — just the choice table, renamed for joining."""
    return choice.select(pl.col(allele_col).alias("_a"), pl.col("p"))


def collapse_alleles(model: Model) -> Model:
    """Return a copy of ``model`` with every gene reduced to a single ``*01`` allele.

    Args:
        model: A V(D)J :class:`Model` (OLGA-derived or EM-learned).

    Returns:
        A new :class:`Model` at gene resolution: one allele (``gene*01``) per gene in every table
        and in the germline. Marginal usage is preserved exactly; conditionals are the correct
        usage-weighted allele averages; germline is the top-usage allele's, relabelled ``*01``.

    Example:
        >>> m = collapse_alleles(load_bundled("TRB", "learned"))
        >>> m.tables["v_choice"].height        # one row per V gene, not per allele
    """
    vdj = model.chain_type == "VDJ"
    t = model.tables
    wv = _allele_weights(t["v_choice"], "v_allele")
    new: dict[str, pl.DataFrame] = {}

    # ---- V ----
    new["v_choice"] = _collapse_choice(t["v_choice"], "v_allele")
    new["v_3_del"] = _collapse_conditional(t["v_3_del"], "v_allele", wv, ["ndel"])

    # ---- J ----
    if vdj:
        # j_choice is the marginal P(J); weights for collapsing J downstream come from it.
        new["j_choice"] = _collapse_choice(t["j_choice"], "j_allele")
        wj = _allele_weights(t["j_choice"], "j_allele")
        new["j_5_del"] = _collapse_conditional(t["j_5_del"], "j_allele", wj, ["ndel"])
    else:
        # VJ: j_choice is P(J|V) keyed (v_allele, j_allele). Collapse J within each V allele
        # (sum over J sub-alleles), then collapse V (weighted average over V sub-alleles).
        jc = (t["j_choice"].with_columns(_gene("v_allele").alias("_gv"), _gene("j_allele").alias("_gj"))
              .group_by(["v_allele", "_gv", "_gj"]).agg(pl.col("p").sum()))   # sum J alleles
        jc = (jc.join(wv.rename({"_a": "v_allele", "_g": "_gv"}), on=["v_allele", "_gv"])
              .with_columns(pw=pl.col("p") * pl.col("w"))
              .group_by(["_gv", "_gj"]).agg(pl.col("pw").sum().alias("p"))    # weighted avg over V alleles
              .with_columns(pl.col("_gv").map_elements(_rep, return_dtype=pl.Utf8).alias("v_allele"),
                            pl.col("_gj").map_elements(_rep, return_dtype=pl.Utf8).alias("j_allele")))
        new["j_choice"] = jc.select("v_allele", "j_allele", "p")
        # j_5_del keyed by j_allele: weight by the J marginal P(J) = Σ_V P(J|V)P(V).
        pv = _marginal(t["v_choice"], "v_allele").rename({"_a": "v_allele"})
        jmarg = (t["j_choice"].join(pv, on="v_allele").with_columns(pj=pl.col("p") * pl.col("p_right"))
                 .group_by("j_allele").agg(pl.col("pj").sum().alias("p")))
        wjm = _allele_weights(jmarg, "j_allele")
        new["j_5_del"] = _collapse_conditional(t["j_5_del"], "j_allele", wjm, ["ndel"])

    # ---- D (VDJ only) ----
    if vdj:
        # d_gene is P(D | J), keyed (j_allele, d_allele). Collapse D (sum sub-alleles) then J
        # (weighted average over J sub-alleles, by j_choice usage).
        wj_full = _allele_weights(t["j_choice"], "j_allele")
        dg = (t["d_gene"].with_columns(_gene("j_allele").alias("_gj"), _gene("d_allele").alias("_gd"))
              .group_by(["j_allele", "_gj", "_gd"]).agg(pl.col("p").sum()))     # sum D alleles
        dg = (dg.join(wj_full.rename({"_a": "j_allele", "_g": "_gj"}), on=["j_allele", "_gj"])
              .with_columns(pw=pl.col("p") * pl.col("w"))
              .group_by(["_gj", "_gd"]).agg(pl.col("pw").sum().alias("p"))       # weighted avg over J alleles
              .with_columns(pl.col("_gj").map_elements(_rep, return_dtype=pl.Utf8).alias("j_allele"),
                            pl.col("_gd").map_elements(_rep, return_dtype=pl.Utf8).alias("d_allele")))
        new["d_gene"] = dg.select("j_allele", "d_allele", "p")
        # d_del keyed by d_allele: weight by the D marginal P(D) = Σ_J P(D|J)P(J).
        pj = _marginal(t["j_choice"], "j_allele").rename({"_a": "j_allele"})
        dmarg = (t["d_gene"].join(pj, on="j_allele").with_columns(pd=pl.col("p") * pl.col("p_right"))
                 .group_by("d_allele").agg(pl.col("pd").sum().alias("p")))
        wdm = _allele_weights(dmarg, "d_allele")
        bins = [c for c in t["d_del"].columns if c not in ("d_allele", "p")]
        new["d_del"] = _collapse_conditional(t["d_del"], "d_allele", wdm, bins)

        # ---- tandem D (d2_gene P(D2|D1), d2_del P(delD2|D2)) — same D allele set ----
        if "d2_gene" in t:
            wd1 = _allele_weights(dmarg, "d_allele")                     # weights over the parent D1
            d2 = (t["d2_gene"].with_columns(_gene("d_allele").alias("_gd1"), _gene("d2_allele").alias("_gd2"))
                  .group_by(["d_allele", "_gd1", "_gd2"]).agg(pl.col("p").sum()))          # sum D2 sub-alleles
            d2 = (d2.join(wd1.rename({"_a": "d_allele", "_g": "_gd1"}), on=["d_allele", "_gd1"])
                  .with_columns(pw=pl.col("p") * pl.col("w"))
                  .group_by(["_gd1", "_gd2"]).agg(pl.col("pw").sum().alias("p"))            # weighted avg over D1
                  .with_columns(pl.col("_gd1").map_elements(_rep, return_dtype=pl.Utf8).alias("d_allele"),
                                pl.col("_gd2").map_elements(_rep, return_dtype=pl.Utf8).alias("d2_allele")))
            new["d2_gene"] = d2.select("d_allele", "d2_allele", "p")
        if "d2_del" in t:
            # weight by the D2 marginal P(D2) = Σ_D1 P(D2|D1)P(D1)
            pd1 = dmarg.select(pl.col("d_allele").alias("_a"), "p")
            d2marg = (t["d2_gene"].join(pd1.rename({"_a": "d_allele"}), on="d_allele")
                      .with_columns(pd2=pl.col("p") * pl.col("p_right"))
                      .group_by("d2_allele").agg(pl.col("pd2").sum().alias("p")))
            wd2 = _allele_weights(d2marg, "d2_allele")
            b2 = [c for c in t["d2_del"].columns if c not in ("d2_allele", "p")]
            new["d2_del"] = _collapse_conditional(t["d2_del"], "d2_allele", wd2, b2)

    # ---- remaining tables must be allele-independent; guard rather than silently pass through ----
    for name in t:
        if name in new:
            continue
        allele_cols = [c for c in t[name].columns if c.endswith("allele")]
        if allele_cols:
            raise NotImplementedError(
                f"collapse_alleles: table {name!r} is allele-keyed ({allele_cols}) but has no "
                f"collapse rule — passing it through would leave the model inconsistent"
            )
        new[name] = t[name]

    # ---- germline: keep the top-usage allele per gene, relabel *01 ----
    pv_all = dict(zip(t["v_choice"]["v_allele"].to_list(), t["v_choice"]["p"].to_list()))
    genomic = {}
    for gname, g in model.genomic.items():
        seg = gname.split("_")[1][0]                    # genes_v -> 'v'
        acol = f"{seg}_allele"
        usage = pv_all if seg == "v" else None
        rows = []
        for gene, sub in g.group_by(_gene(acol).alias("_g"), maintain_order=True):
            key = gene[0] if isinstance(gene, tuple) else gene
            # representative = highest-usage allele (V), else the first (stable) — usually *01
            if usage is not None:
                best = max(sub[acol].to_list(), key=lambda a: usage.get(a, 0.0))
            else:
                best = sorted(sub[acol].to_list())[0]
            rep = sub.filter(pl.col(acol) == best).with_columns(pl.lit(_rep(key)).alias(acol),
                                                                pl.lit(key).alias("gene"))
            rows.append(rep)
        genomic[gname] = pl.concat(rows)

    return Model(manifest=model.manifest, tables=new, genomic=genomic)
