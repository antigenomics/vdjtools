"""Seed list: clonotypes in >=N COVID donors and ZERO controls -> re-test.

The full 1.9M-candidate scan spends its whole FDR budget on clonotypes that were never
plausible. Seeding on a hard incidence prior -- present in >=5 COVID donors, absent from every
control -- cuts the tested set by orders of magnitude, so BH has budget left for real hits.

Scopes: exact | aa (CDR3aa+-1mm) | aa_v (V + CDR3aa+-1mm).

THREE units, because depth is a genuine confound and the first two normalize it differently:

  donor  : Fisher on presence/absence. Ignores depth entirely -- and depth is NOT balanced
           here (controls carry ~40% (TRA) / ~52% (TRB) more rearrangements per donor), so a
           control is more likely to carry ANY clonotype. This biases the donor test AGAINST
           COVID enrichment.

  rows   : conditional binomial with p0 = N_pos/(N_pos+N_neg) -- the ROW ratio (0.5315 TRA).
           Normalizes total sequencing effort. BUT: for a clonotype whose row count does not
           scale with depth (public, carried ~once per donor), the correct null is the DONOR
           ratio 761/1240 = 0.6137. Those differ by 1.40x, so every such clonotype gets a
           spurious 1.4x "enrichment" -- which is why this test enriched VDJdb CMV as hard as
           SARS-CoV-2. Reported here with BOTH nulls so the artifact is visible, not argued.

  freq   : per-donor frequency f_i = n_rearr_i / N_i, compared across arms (Mann-Whitney U on
           donors, zeros included). This is the synthesis: the REARRANGEMENT is the unit (so
           depth is normalized, per donor rather than per arm) while the DONOR stays the
           independent replicate (so no pseudoreplication and no row/donor-ratio mismatch).

VDJdb SARS-CoV-2 vs CMV decides which unit is real. 2026-07-16.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import binom, fisher_exact, mannwhitneyu

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19")
RES = Path("/projects/biomarkers/results")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")
MIN_READS = 10_000


def oracle(chain, pattern):
    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    return set(db.filter((pl.col("gene") == chain) & (pl.col("species") == "HomoSapiens")
                         & pl.col("antigen.species").str.contains(pattern))["cdr3"]
               .unique().to_list())


def enrich(sig, rest, truth):
    a = len(sig & truth); b = len(sig) - a
    c = len(rest & truth); d = len(rest) - c
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, orr, p


def load(meta, chain):
    frames = []
    for r in meta.filter(pl.col("locus") == chain).iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{chain}.*"))
        if hits:
            frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                          .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    return pl.concat(frames, how="vertical_relaxed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRA")
    ap.add_argument("--scope", default="aa_v", choices=["exact", "aa", "aa_v"])
    ap.add_argument("--min-covid", type=int, default=5)
    ap.add_argument("--max-control", type=int, default=0)
    args = ap.parse_args()
    print(f"=== SEED {args.chain} scope={args.scope}: >={args.min_covid} COVID donors, "
          f"<={args.max_control} control donors ===", flush=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    meta = meta.with_columns(pl.col("reads").cast(pl.Int64, strict=False)).filter(
        pl.col("reads") >= MIN_READS)
    design = (meta.filter(pl.col("locus") == args.chain)
              .select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                      pl.when(pl.col("COVID_status") == "COVID").then(True)
                        .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                        .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))

    t0 = time.perf_counter()
    cohort = load(meta, args.chain)
    vgene = pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene")
    rows = (cohort.lazy()
            .filter(pl.col(S.JUNCTION_AA).is_not_null()
                    & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .with_columns(vgene)
            .group_by([S.JUNCTION_AA, "v_gene", "sample_id"])
            .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
            .join(design.lazy(), on="sample_id", how="inner")
            .collect(engine="streaming"))
    print(f"[ingest {time.perf_counter()-t0:.0f}s]", flush=True)

    depth = (rows.group_by(["sample_id", "_pos"]).agg(pl.col("n_rearr").sum().alias("N_i")))
    n_pos = int(depth["_pos"].sum()); n_neg = depth.height - n_pos
    N_pos = int(depth.filter(pl.col("_pos"))["N_i"].sum())
    N_neg = int(depth.filter(~pl.col("_pos"))["N_i"].sum())
    p0_rows = N_pos / (N_pos + N_neg)
    p0_donor = n_pos / (n_pos + n_neg)
    print(f"  donors {n_pos}+/{n_neg}-;  rows {N_pos:,}+/{N_neg:,}-", flush=True)
    print(f"  depth/donor: COVID+ {N_pos/n_pos:,.0f}  control {N_neg/n_neg:,.0f}  "
          f"(control/{'COVID':5s} = {(N_neg/n_neg)/(N_pos/n_pos):.2f}x)", flush=True)
    print(f"  NULLS: p0_rows={p0_rows:.4f}  p0_donor={p0_donor:.4f}  "
          f"ratio={p0_donor/p0_rows:.3f}  <- the row test's bias for non-depth-scaling clonotypes",
          flush=True)

    # ---- per-candidate donor sets & rows, under the chosen scope -------------------------
    if args.scope == "exact":
        j = rows.rename({S.JUNCTION_AA: "cand"})
        gcols = ["cand"]
    else:
        import vdjmatch.cluster as vc
        universe = rows[S.JUNCTION_AA].unique().sort().to_list()
        t0 = time.perf_counter()
        pairs = vc.overlap(universe, universe, scope="1,0,0,1", threads=0)
        print(f"[1mm search {time.perf_counter()-t0:.0f}s] pairs {pairs.height:,}", flush=True)
        umap = pl.DataFrame({"idx": np.arange(len(universe), dtype=np.int64),
                             "aa": universe})
        nb = (pairs.join(umap.rename({"idx": "a_idx", "aa": "cand"}), on="a_idx")
              .join(umap.rename({"idx": "b_idx", "aa": S.JUNCTION_AA}), on="b_idx")
              .select("cand", S.JUNCTION_AA))
        if args.scope == "aa":
            j = nb.join(rows, on=S.JUNCTION_AA, how="inner")
            gcols = ["cand"]
        else:
            # candidate = (aa, V); a donor matches only with the SAME V gene
            cv = rows.select(pl.col(S.JUNCTION_AA).alias("cand"), "v_gene").unique()
            j = (nb.join(cv, on="cand", how="inner")
                 .join(rows, on=[S.JUNCTION_AA, "v_gene"], how="inner"))
            gcols = ["cand", "v_gene"]

    agg = (j.group_by(gcols).agg(
        pl.col("sample_id").filter(pl.col("_pos")).n_unique().alias("a"),
        pl.col("sample_id").filter(~pl.col("_pos")).n_unique().alias("b"),
        pl.col("n_rearr").filter(pl.col("_pos")).sum().alias("k_pos"),
        pl.col("n_rearr").filter(~pl.col("_pos")).sum().alias("k_neg")))

    # ---- THE SEED FILTER ------------------------------------------------------------------
    seed = agg.filter((pl.col("a") >= args.min_covid) & (pl.col("b") <= args.max_control))
    print(f"\n  SEED LIST: {seed.height:,} clonotypes "
          f"(from {agg.height:,} tested = {seed.height/max(agg.height,1):.4%})", flush=True)
    if not seed.height:
        return

    a = seed["a"].to_numpy().astype(np.int64)
    b = seed["b"].to_numpy().astype(np.int64)
    k_pos = seed["k_pos"].fill_null(0).to_numpy().astype(np.int64)
    k_neg = seed["k_neg"].fill_null(0).to_numpy().astype(np.int64)
    c, d = n_pos - a, n_neg - b

    p_donor = stats.fisher_p(a, b, c, d, alternative="greater")
    k = k_pos + k_neg
    p_rows = binom.sf(k_pos - 1, k, p0_rows)          # row-ratio null (what I used before)
    p_rows_dn = binom.sf(k_pos - 1, k, p0_donor)      # donor-ratio null (conservative)

    out = seed.with_columns(
        pl.Series("p_donor", p_donor), pl.Series("q_donor", stats.fdr_bh(p_donor)),
        pl.Series("p_rows", p_rows), pl.Series("q_rows", stats.fdr_bh(p_rows)),
        pl.Series("p_rows_donornull", p_rows_dn),
        pl.Series("q_rows_donornull", stats.fdr_bh(p_rows_dn)))

    print(f"\n  --- FDR on the seed list ({seed.height:,} tests, not {agg.height:,}) ---",
          flush=True)
    for lab, col in (("donor  (Fisher)", "q_donor"), ("rows   (p0=row ratio)", "q_rows"),
                     ("rows   (p0=donor ratio)", "q_rows_donornull")):
        print(f"    {lab:26s} q<0.01 {out.filter(pl.col(col) < 0.01).height:6,}   "
              f"q<0.05 {out.filter(pl.col(col) < 0.05).height:6,}", flush=True)

    # ---- VDJdb: which unit is REAL -------------------------------------------------------
    cov = oracle(args.chain, "SARS-CoV-2")
    cmv = oracle(args.chain, "CMV")
    allc = set(agg["cand"].to_list())
    seedset = set(out["cand"].to_list())
    print(f"\n  --- VDJdb (SARS-CoV-2 {len(cov):,} | CMV {len(cmv):,} control) ---", flush=True)
    aa_, orr, pv = enrich(seedset, allc - seedset, cov)
    a2, or2, p2 = enrich(seedset, allc - seedset, cmv)
    print(f"    SEED LIST itself       : CoV2 {aa_:4d} OR={orr:7.2f} p={pv:9.2e}  |  "
          f"CMV {a2:4d} OR={or2:7.2f} p={p2:9.2e}", flush=True)
    for col, lab in (("q_donor", "donor"), ("q_rows", "rows(rowNull)"),
                     ("q_rows_donornull", "rows(donorNull)")):
        sig = set(out.filter(pl.col(col) < 0.01)["cand"].to_list())
        if not sig:
            print(f"    {lab:22s}: 0 significant", flush=True)
            continue
        aa_, orr, pv = enrich(sig, allc - sig, cov)
        a2, or2, p2 = enrich(sig, allc - sig, cmv)
        spec = (aa_ / max(len(sig & set(allc)), 1)) / max(a2 / max(len(sig), 1), 1e-9)
        print(f"    {lab:22s}: n={len(sig):5,}  CoV2 {aa_:4d} OR={orr:7.2f}  |  "
              f"CMV {a2:4d} OR={or2:7.2f}   CoV2/CMV OR ratio={orr/max(or2,1e-9):5.2f}",
              flush=True)

    out.sort("p_donor").write_parquet(RES / f"seed_{args.chain}_{args.scope}.parquet")
    print(f"\n  top 15 of the seed list by donor p:", flush=True)
    print(out.sort("p_donor").head(15).select(gcols + ["a", "b", "k_pos", "k_neg", "p_donor",
                                                       "q_donor", "q_rows"]), flush=True)
    print(f"\n-> {RES}/seed_{args.chain}_{args.scope}.parquet", flush=True)


if __name__ == "__main__":
    main()
