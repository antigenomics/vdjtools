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
import os
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


def _keys_of(path: str):
    """Pass 1 worker: unique (CDR3aa, V) in ONE repertoire. Dedup is per-donor -- a key counts
    once per donor, which is exactly what an incidence pool needs."""
    df = (vio.read(path, fmt="vdjtools")
          .filter(pl.col(S.JUNCTION_AA).is_not_null()
                  & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
          .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
          .select(S.JUNCTION_AA, "v_gene").unique())
    return list(zip(df[S.JUNCTION_AA].to_list(), df["v_gene"].to_list()))


def _rows_of(args):
    """Pass 2 worker: ONE repertoire, candidate keys only."""
    path, sid, keep = args
    df = (vio.read(path, fmt="vdjtools")
          .filter(pl.col(S.JUNCTION_AA).is_not_null()
                  & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
          .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
          .group_by([S.JUNCTION_AA, "v_gene"])
          .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
          .join(keep, on=[S.JUNCTION_AA, "v_gene"], how="semi"))
    return df.with_columns(pl.lit(sid).alias("sample_id")) if df.height else None


def pool_incidence(paths, threads):
    """Pass 1 -- POOLING: key -> #donors carrying it. Streamed per file, parallel.

    Memory is O(unique keys), NOT O(total rows). hip is 144.8M (CDR3aa,V,donor) rows over 761
    immunoSEQ repertoires (~190k clonotypes/donor). Concatenating them and THEN calling
    group_by blows 370G and stalls before vdjmatch is ever reached (measured: jobs 1417286 OOM,
    1417291 stalled 41min). Only the pooled counter stays resident: each file is read, deduped
    within the donor, folded in, and dropped.
    """
    from collections import Counter
    from concurrent.futures import ProcessPoolExecutor
    total = Counter()
    with ProcessPoolExecutor(max_workers=threads) as ex:
        for i, keys in enumerate(ex.map(_keys_of, paths, chunksize=4), 1):
            total.update(keys)
            if i % 100 == 0:
                print(f"    pooled {i}/{len(paths)} donors; {len(total):,} unique keys",
                      flush=True)
    return total


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
    ap.add_argument("--min-donors", type=int, default=2,
                    help="pool cutoff: a key must be seen in >=N donors")
    ap.add_argument("--scope", default="1mm", choices=["exact", "1mm"],
                    help="1mm is the default. `exact` additionally reproduces Emerson faithfully "
                         "-- his 164 TCRb came from exact V+CDR3aa incidence; 1mm is Vlasova's "
                         "addition -- so it is worth running as its own arm, not as a fallback.")
    args = ap.parse_args()
    key = "V+CDR3aa" if args.scope == "exact" else "V+CDR3aa±1mm"
    print(f"=== hip CMV {args.arm}  key={key} ===", flush=True)

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

    files = [(str(ROOT / r["file_name"]), str(r["sample_id"]))
             for r in design.iter_rows(named=True) if (ROOT / r["file_name"]).exists()]
    nthreads = os.cpu_count() or 8

    # PASS 1 -- pool incidence per key, streamed. Never materialize the cohort.
    t0 = time.perf_counter()
    total = pool_incidence([p for p, _ in files], nthreads)
    keep = pl.DataFrame(
        {S.JUNCTION_AA: [k[0] for k, n in total.items() if n >= args.min_donors],
         "v_gene": [k[1] for k, n in total.items() if n >= args.min_donors]})
    print(f"[pass1 pool {time.perf_counter()-t0:.0f}s] {len(total):,} unique keys -> "
          f"{keep.height:,} public (>={args.min_donors} donors)", flush=True)
    del total

    # PASS 2 -- per-donor rows, candidates only (a small fraction of the cohort).
    t0 = time.perf_counter()
    from concurrent.futures import ProcessPoolExecutor
    parts = []
    with ProcessPoolExecutor(max_workers=nthreads) as ex:
        for d in ex.map(_rows_of, [(p, s, keep) for p, s in files], chunksize=4):
            if d is not None:
                parts.append(d)
    rows = (pl.concat(parts, how="vertical_relaxed")
            .join(design.select("sample_id", "_pos"), on="sample_id", how="inner"))
    del parts
    print(f"[pass2 rows {time.perf_counter()-t0:.0f}s] rows {rows.height:,}", flush=True)
    # NOTE the >=min_donors pool is the paper's own candidate rule, but it is also an
    # APPROXIMATION of fuzzy incidence, stated not hidden: a PRIVATE neighbour (1 donor) can no
    # longer contribute its donor to a candidate's ball. It cannot manufacture an association --
    # the loss applies to both arms -- but it does slightly under-count incidence.

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
