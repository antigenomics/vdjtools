"""A3 — vaccination timepoint (airr_covid19_vacc): before vs after.

Owner's caveat, which sets the prior: the nominal timepoint is d15-20, but donors can be 1-2
MONTHS out -- i.e. past the effector peak, like the convalescent COVID cohort. A weak effect is
expected. This arm is worth running anyway because it is PAIRED where the COVID arm is not: the
same donor before and after, so donor-level confounders (HLA, depth, batch, prior exposure)
cancel within the pair.

Two contrasts:
  timepoint : before vs after vaccination (all donors)
  vaccine   : GamCOVIDVac vs CoviVac, restricted to the after-vaccination samples

Paired structure is exploited where the metadata supports it: a clonotype gained by a donor
between timepoints is far stronger evidence than a clonotype merely more common in the "after"
arm, because the same repertoire is its own control. Reported as `n_gained` (absent before,
present after) vs `n_lost`, tested with McNemar -- the paired analogue of Fisher, and the right
test when each subject contributes both arms.

Key = V + CDR3aa±1mm; units = donor and rows (donor-ratio null), as everywhere else. 2026-07-17.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import binom, binomtest

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19_vacc")
RES = Path("/projects/biomarkers/results")
MIN_READS = 10_000


def _keys_of(path):
    df = (vio.read(path, fmt="vdjtools")
          .filter(pl.col(S.JUNCTION_AA).is_not_null()
                  & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
          .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
          .select(S.JUNCTION_AA, "v_gene").unique())
    return list(zip(df[S.JUNCTION_AA].to_list(), df["v_gene"].to_list()))


def _rows_of(a):
    path, sid, keep = a
    df = (vio.read(path, fmt="vdjtools")
          .filter(pl.col(S.JUNCTION_AA).is_not_null()
                  & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
          .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
          .group_by([S.JUNCTION_AA, "v_gene"])
          .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
          .join(keep, on=[S.JUNCTION_AA, "v_gene"], how="semi"))
    return df.with_columns(pl.lit(sid).alias("sample_id")) if df.height else None


def pool_incidence(paths, threads):
    """Streamed per-file pooling; THREADS not processes (polars drops the GIL on read, and
    forked workers each spin their own thread pool -> ~1600 threads and a stall)."""
    from collections import Counter
    from concurrent.futures import ThreadPoolExecutor
    total = Counter()
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for i, k in enumerate(ex.map(_keys_of, paths), 1):
            total.update(k)
            if i % 200 == 0:
                print(f"    pooled {i}/{len(paths)}; {len(total):,} keys "
                      f"[{time.perf_counter()-t0:.0f}s]", flush=True)
    return total


def main():
    import os
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    ap.add_argument("--arm", default="timepoint", choices=["timepoint", "vaccine"])
    ap.add_argument("--min-donors", type=int, default=2)
    args = ap.parse_args()
    print(f"=== vacc {args.chain} {args.arm}  key=V+CDR3aa±1mm ===", flush=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    cols = meta.columns
    print(f"  metadata: {meta.height} rows; cols={cols[:12]}", flush=True)
    m = meta.filter(pl.col("file_name").str.contains(args.chain)) if "file_name" in cols else meta
    print(f"  {args.chain} rows: {m.height}", flush=True)
    print("  timepoint:", m.group_by("timepoint").len().sort("len", descending=True).to_dicts(),
          flush=True)
    print("  vaccine:", m.group_by("vaccine").len().sort("len", descending=True).to_dicts(),
          flush=True)

    if args.arm == "timepoint":
        design = (m.select(pl.col("full_id").cast(pl.String).alias("sample_id"),
                           pl.col("name").cast(pl.String).alias("donor"),
                           pl.when(pl.col("timepoint").str.contains("(?i)after")).then(True)
                             .when(pl.col("timepoint").str.contains("(?i)before")).then(False)
                             .otherwise(None).alias("_pos"))
                  .drop_nulls("_pos").unique(subset="sample_id"))
    else:
        design = (m.filter(pl.col("timepoint").str.contains("(?i)after"))
                  .select(pl.col("full_id").cast(pl.String).alias("sample_id"),
                          pl.col("name").cast(pl.String).alias("donor"),
                          pl.when(pl.col("vaccine") == "GamCOVIDVac").then(True)
                            .when(pl.col("vaccine") == "CoviVac").then(False)
                            .otherwise(None).alias("_pos"))
                  .drop_nulls("_pos").unique(subset="sample_id"))
    n_pos = int(design["_pos"].sum())
    print(f"  design: {n_pos}+ / {design.height-n_pos}-", flush=True)
    if not (n_pos and design.height - n_pos):
        print("  degenerate arms -- stop", flush=True)
        return

    # map sample_id -> file
    fmap = {}
    for r in m.iter_rows(named=True):
        for cand in (r.get("full_id"), r.get("name")):
            if cand:
                hits = sorted(ROOT.glob(f"*{cand}*.{args.chain}.*"))
                if hits:
                    fmap[str(r["full_id"])] = str(hits[0])
                    break
    files = [(fmap[s], s) for s in design["sample_id"].to_list() if s in fmap]
    print(f"  files resolved: {len(files)}/{design.height}", flush=True)
    if len(files) < 20:
        print("  too few files resolved -- check the file_name/full_id map", flush=True)
        return

    nthreads = os.cpu_count() or 8
    total = pool_incidence([p for p, _ in files], nthreads)
    keep = pl.DataFrame(
        {S.JUNCTION_AA: [k[0] for k, n in total.items() if n >= args.min_donors],
         "v_gene": [k[1] for k, n in total.items() if n >= args.min_donors]})
    print(f"  pooled {len(total):,} keys -> {keep.height:,} public", flush=True)
    del total

    from concurrent.futures import ThreadPoolExecutor
    parts = []
    with ThreadPoolExecutor(max_workers=nthreads) as ex:
        for d in ex.map(_rows_of, [(p, s, keep) for p, s in files]):
            if d is not None:
                parts.append(d)
    rows = pl.concat(parts, how="vertical_relaxed").join(design, on="sample_id", how="inner")
    del parts
    print(f"  rows {rows.height:,}", flush=True)

    import vdjmatch.cluster as vc
    universe = rows[S.JUNCTION_AA].unique().sort().to_list()
    pairs = vc.overlap(universe, universe, scope="1,0,0,1", threads=0)
    um = pl.DataFrame({"idx": np.arange(len(universe), dtype=np.int64), "aa": universe})
    nb = (pairs.join(um.rename({"idx": "a_idx", "aa": "cand"}), on="a_idx")
          .join(um.rename({"idx": "b_idx", "aa": S.JUNCTION_AA}), on="b_idx")
          .select("cand", S.JUNCTION_AA))
    print(f"  1mm pairs {nb.height:,}", flush=True)

    dep = rows.group_by(["sample_id", "_pos"]).agg(pl.col("n_rearr").sum().alias("N_i"))
    np_, nn = int(dep["_pos"].sum()), dep.height - int(dep["_pos"].sum())
    N_pos = int(dep.filter(pl.col("_pos"))["N_i"].sum())
    N_neg = int(dep.filter(~pl.col("_pos"))["N_i"].sum())
    p0_donor = np_ / (np_ + nn)
    print(f"  analysed {np_}+/{nn}-;  depth {N_pos/np_:,.0f} vs {N_neg/nn:,.0f}", flush=True)

    cv = rows.select(pl.col(S.JUNCTION_AA).alias("cand"), "v_gene").unique()
    j = nb.join(cv, on="cand", how="inner").join(rows, on=[S.JUNCTION_AA, "v_gene"], how="inner")
    agg = (j.group_by(["cand", "v_gene"]).agg(
        pl.col("sample_id").filter(pl.col("_pos")).n_unique().alias("a"),
        pl.col("sample_id").filter(~pl.col("_pos")).n_unique().alias("b"),
        pl.col("n_rearr").filter(pl.col("_pos")).sum().alias("k_pos"),
        pl.col("n_rearr").filter(~pl.col("_pos")).sum().alias("k_neg"),
        pl.col("donor").filter(pl.col("_pos")).unique().alias("d_pos"),
        pl.col("donor").filter(~pl.col("_pos")).unique().alias("d_neg"))
        .filter(pl.col("a") + pl.col("b") >= 2))
    a = agg["a"].to_numpy().astype(np.int64)
    b = agg["b"].to_numpy().astype(np.int64)
    kp = agg["k_pos"].fill_null(0).to_numpy().astype(np.int64)
    kn = agg["k_neg"].fill_null(0).to_numpy().astype(np.int64)
    p_d = stats.fisher_p(a, b, np_ - a, nn - b, alternative="greater")
    p_rd = binom.sf(kp - 1, kp + kn, p0_donor)
    res = agg.with_columns(
        pl.Series("odds_ratio", stats.odds_ratio(a, b, np_ - a, nn - b)),
        pl.Series("p_donor", p_d), pl.Series("q_donor", stats.fdr_bh(p_d)),
        pl.Series("p_rows_dn", p_rd), pl.Series("q_rows_dn", stats.fdr_bh(p_rd)))

    if args.arm == "timepoint":
        # PAIRED: gained = donor has it after but not before. McNemar on discordant pairs only,
        # which is the point -- the donor is its own control, so HLA/depth/prior exposure cancel.
        gained = res.with_columns(
            pl.col("d_pos").list.set_difference(pl.col("d_neg")).list.len().alias("n_gained"),
            pl.col("d_neg").list.set_difference(pl.col("d_pos")).list.len().alias("n_lost"))
        g = gained["n_gained"].to_numpy()
        l = gained["n_lost"].to_numpy()
        mc = np.array([binomtest(int(gi), int(gi + li), 0.5, alternative="greater").pvalue
                       if (gi + li) >= 5 else 1.0 for gi, li in zip(g, l)])
        res = gained.with_columns(pl.Series("p_mcnemar", mc),
                                  pl.Series("q_mcnemar", stats.fdr_bh(mc)))
        print(f"\n  PAIRED (McNemar, gained vs lost within donor):", flush=True)
        for lab, c in (("q<0.01", 0.01), ("q<0.05", 0.05)):
            print(f"    {lab}: {res.filter(pl.col('q_mcnemar') < c).height:,}", flush=True)

    print(f"\n  features tested {res.height:,}", flush=True)
    for lab, c in (("donor q<0.01", "q_donor"), ("rows(dn) q<0.01", "q_rows_dn")):
        print(f"    {lab}: {res.filter(pl.col(c) < 0.01).height:,}", flush=True)
    show = [c for c in ("cand", "v_gene", "a", "b", "n_gained", "n_lost", "odds_ratio",
                        "p_donor", "q_donor", "p_mcnemar", "q_mcnemar") if c in res.columns]
    print("\n  top 12 by donor p:", flush=True)
    print(res.sort("p_donor").head(12).select(show), flush=True)
    res.drop([c for c in ("d_pos", "d_neg") if c in res.columns]).write_parquet(
        RES / f"vacc_{args.chain}_{args.arm}.parquet")
    print(f"\n-> {RES}/vacc_{args.chain}_{args.arm}.parquet", flush=True)


if __name__ == "__main__":
    main()
