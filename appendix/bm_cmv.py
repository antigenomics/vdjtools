"""A2/H — CMV association on airr_hip (Emerson 2017), plain and HLA-restricted.

CMV is the positive control the whole framework should pass: unlike convalescent COVID, a CMV+
donor carries a large, persistent, HLA-restricted memory compartment at steady state. If the
pipeline cannot find CMV-associated TCRb here it cannot find anything.

Arms:
  status      : CMV+ vs CMV-  (Emerson's 340/421)
  hla         : same contrast restricted to carriers of each top-N HLA allele. hip types HLA at
                1-field / 2-digit only (HLA-A*02, not A*02:01), comma-separated in one column.

Emerson's own threshold was NOMINAL P<1e-4 with a permutation-estimated FDR of 0.14 -- not BH
q<0.05. Both are reported: BH is the stricter modern default, P<1e-4 is what the 164-TCRb
published figure was produced with, so it is the comparable number.

Units: donor (Fisher on incidence) and rows (conditional binomial on rearrangements) under both
nulls, same as the COVID arms. Key = V + CDR3aa±1mm. 2026-07-16.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import binom, fisher_exact

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/hip")
RES = Path("/projects/biomarkers/results")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")


def oracle(pattern):
    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    return set(db.filter((pl.col("gene") == "TRB") & (pl.col("species") == "HomoSapiens")
                         & pl.col("antigen.species").str.contains(pattern))["cdr3"]
               .unique().to_list())


def enrich(sig, rest, truth):
    a = len(sig & truth); b = len(sig) - a
    c = len(rest & truth); d = len(rest) - c
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, orr, p


def run_contrast(rows, nb, label, cov, cmv_ep=None):
    dep = rows.group_by(["sample_id", "_pos"]).agg(pl.col("n_rearr").sum().alias("N_i"))
    n_pos = int(dep["_pos"].sum()); n_neg = dep.height - n_pos
    if not (n_pos and n_neg):
        print(f"  {label}: degenerate arms, skip", flush=True)
        return None
    N_pos = int(dep.filter(pl.col("_pos"))["N_i"].sum())
    N_neg = int(dep.filter(~pl.col("_pos"))["N_i"].sum())
    p0_rows = N_pos / (N_pos + N_neg)
    p0_donor = n_pos / (n_pos + n_neg)

    cv = rows.select(pl.col(S.JUNCTION_AA).alias("cand"), "v_gene").unique()
    j = (nb.join(cv, on="cand", how="inner")
         .join(rows, on=[S.JUNCTION_AA, "v_gene"], how="inner"))
    agg = (j.group_by(["cand", "v_gene"]).agg(
        pl.col("sample_id").filter(pl.col("_pos")).n_unique().alias("a"),
        pl.col("sample_id").filter(~pl.col("_pos")).n_unique().alias("b"),
        pl.col("n_rearr").filter(pl.col("_pos")).sum().alias("k_pos"),
        pl.col("n_rearr").filter(~pl.col("_pos")).sum().alias("k_neg"))
        .filter(pl.col("a") + pl.col("b") >= 2))
    a = agg["a"].to_numpy().astype(np.int64)
    b = agg["b"].to_numpy().astype(np.int64)
    kp = agg["k_pos"].fill_null(0).to_numpy().astype(np.int64)
    kn = agg["k_neg"].fill_null(0).to_numpy().astype(np.int64)
    p_d = stats.fisher_p(a, b, n_pos - a, n_neg - b, alternative="greater")
    p_r = binom.sf(kp - 1, kp + kn, p0_rows)
    p_rd = binom.sf(kp - 1, kp + kn, p0_donor)
    res = agg.with_columns(
        pl.lit(label).alias("stratum"),
        pl.Series("odds_ratio", stats.odds_ratio(a, b, n_pos - a, n_neg - b)),
        pl.Series("p_donor", p_d), pl.Series("q_donor", stats.fdr_bh(p_d)),
        pl.Series("p_rows", p_r), pl.Series("q_rows", stats.fdr_bh(p_r)),
        pl.Series("p_rows_dn", p_rd), pl.Series("q_rows_dn", stats.fdr_bh(p_rd)))

    n_p4 = res.filter(pl.col("p_donor") < 1e-4).height
    print(f"\n  --- {label}: {n_pos}+/{n_neg}-;  {res.height:,} tested; "
          f"depth {N_pos/n_pos:,.0f} vs {N_neg/n_neg:,.0f}", flush=True)
    print(f"      donor: P<1e-4 {n_p4:6,} (Emerson's threshold; his TRb figure = 164)  |  "
          f"q<0.01 {res.filter(pl.col('q_donor') < 0.01).height:6,}  "
          f"q<0.05 {res.filter(pl.col('q_donor') < 0.05).height:6,}", flush=True)
    print(f"      rows : q<0.01 rowNull {res.filter(pl.col('q_rows') < 0.01).height:6,}  |  "
          f"donorNull {res.filter(pl.col('q_rows_dn') < 0.01).height:6,}", flush=True)
    allc = set(res["cand"].to_list())
    for col, thr, nm in (("p_donor", 1e-4, "donor P<1e-4"), ("q_donor", 0.01, "donor q<0.01"),
                         ("q_rows_dn", 0.01, "rows(dn) q<0.01")):
        sig = set(res.filter(pl.col(col) < thr)["cand"].to_list())
        if sig:
            h, orr, pv = enrich(sig, allc - sig, cov)
            print(f"      VDJdb CMV among {nm:16s}: {h:4d}/{len(sig):<6,} OR={orr:7.2f} "
                  f"p={pv:9.2e}", flush=True)
    if n_p4:
        print(res.filter(pl.col("p_donor") < 1e-4).sort("p_donor").head(6)
              .select("cand", "v_gene", "a", "b", "odds_ratio", "p_donor"), flush=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="status", choices=["status", "hla"])
    ap.add_argument("--top", type=int, default=5)
    args = ap.parse_args()
    print(f"=== hip CMV {args.arm}  key=V+CDR3aa±1mm ===", flush=True)

    meta = pl.read_csv(ROOT / "metadata.txt", separator="\t", infer_schema_length=0)
    design = (meta.select(pl.col("sample_id").cast(pl.String),
                          pl.col("file_name"), pl.col("hla"),
                          pl.when(pl.col("cmv") == "+").then(True)
                            .when(pl.col("cmv") == "-").then(False)
                            .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))
    print(f"  CMV-typed donors: {design.height} "
          f"({int(design['_pos'].sum())}+ / {design.height-int(design['_pos'].sum())}-)",
          flush=True)

    t0 = time.perf_counter()
    frames = []
    for r in design.iter_rows(named=True):
        f = ROOT / r["file_name"]
        if f.exists():
            frames.append(vio.read(str(f), fmt="vdjtools")
                          .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    cohort = pl.concat(frames, how="vertical_relaxed")
    rows = (cohort.lazy()
            .filter(pl.col(S.JUNCTION_AA).is_not_null()
                    & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
            .group_by([S.JUNCTION_AA, "v_gene", "sample_id"])
            .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
            .join(design.select("sample_id", "_pos").lazy(), on="sample_id", how="inner")
            .collect(engine="streaming"))
    print(f"[ingest {time.perf_counter()-t0:.0f}s] rows {rows.height:,}", flush=True)

    # hip is immunoSEQ-deep (~190k rows/donor vs covid's ~14k): an all-vs-all 1mm search over
    # every CDR3aa OOMs at 370G. Restrict to CDR3aa seen in >=2 donors -- the paper's own
    # candidate rule. APPROXIMATION, stated not hidden: a private neighbour (1 donor) can no
    # longer contribute its donor to a candidate's fuzzy incidence, so incidence is a slight
    # undercount. It cannot manufacture an association: the loss applies to both arms.
    pub = (rows.group_by(S.JUNCTION_AA).agg(pl.col("sample_id").n_unique().alias("d"))
           .filter(pl.col("d") >= 2))
    rows = rows.join(pub.select(S.JUNCTION_AA), on=S.JUNCTION_AA, how="semi")
    print(f"  public prefilter (>=2 donors): {pub.height:,} CDR3aa; rows -> {rows.height:,}",
          flush=True)

    import vdjmatch.cluster as vc
    universe = rows[S.JUNCTION_AA].unique().sort().to_list()
    t0 = time.perf_counter()
    pairs = vc.overlap(universe, universe, scope="1,0,0,1", threads=0)
    umap = pl.DataFrame({"idx": np.arange(len(universe), dtype=np.int64), "aa": universe})
    nb = (pairs.join(umap.rename({"idx": "a_idx", "aa": "cand"}), on="a_idx")
          .join(umap.rename({"idx": "b_idx", "aa": S.JUNCTION_AA}), on="b_idx")
          .select("cand", S.JUNCTION_AA))
    print(f"[1mm search {time.perf_counter()-t0:.0f}s] pairs {nb.height:,}", flush=True)

    cov = oracle("CMV")
    print(f"  VDJdb CMV TRB: {len(cov):,}", flush=True)
    out = []
    if args.arm == "status":
        r = run_contrast(rows, nb, "CMV+ vs CMV- (all donors)", cov)
        if r is not None:
            out.append(r.filter((pl.col("p_donor") < 1e-4) | (pl.col("q_rows_dn") < 0.05)))
    else:
        # hip HLA: 2-digit, comma-separated in one column
        h = (design.filter(pl.col("hla").is_not_null() & (pl.col("hla") != "")
                           & (pl.col("hla") != "NA"))
             .with_columns(pl.col("hla").str.split(",").alias("_al"))
             .explode("_al").with_columns(pl.col("_al").str.strip_chars().alias("allele"))
             .filter(pl.col("allele") != "").unique(["sample_id", "allele"]))
        top = (h.group_by("allele").agg(pl.len().alias("carriers"),
                                        pl.col("_pos").sum().alias("cmv_pos"))
               .with_columns((pl.col("carriers") - pl.col("cmv_pos")).alias("cmv_neg"))
               .filter((pl.col("cmv_pos") >= 20) & (pl.col("cmv_neg") >= 20))
               .sort("carriers", descending=True).head(args.top))
        print("  HLA strata (hip is 1-field):", flush=True)
        print(top, flush=True)
        for row in top.iter_rows(named=True):
            who = set(h.filter(pl.col("allele") == row["allele"])["sample_id"].to_list())
            r = run_contrast(rows.filter(pl.col("sample_id").is_in(who)), nb,
                             f"CMV | {row['allele']}", cov)
            if r is not None:
                out.append(r.filter((pl.col("p_donor") < 1e-4) | (pl.col("q_rows_dn") < 0.05)))
    if out:
        allres = pl.concat(out, how="vertical_relaxed")
        f = RES / f"cmv_hip_{args.arm}.parquet"
        allres.write_parquet(f)
        print(f"\n-> {f} ({allres.height:,} rows)", flush=True)


if __name__ == "__main__":
    main()
