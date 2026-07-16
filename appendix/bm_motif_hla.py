"""Which HLA actually restricts our top hits? Ask the motif, not the stratum.

The stratified scan says the COVID alpha signal lives in A*02 carriers. VDJdb says our top
pair (CAMRDW.LTGGGNKLTF / TRAV14DV4 + CASSL..PQETQYF / TRBV27) is a SARS-CoV-2 receptor against
TTDPSFLGRY restricted by HLA-A*01:01. Both cannot be right as stated, and a stratified hit
count is a weak instrument for this: A*02 is carried by ~half the cohort, so ANY signal leaks
into every other stratum through co-carriage, and BH boundaries move whole related families at
once.

So test each motif directly: among COVID+ donors ONLY, is carriage of the motif associated with
carriage of each HLA allele? Restricting to COVID+ removes the infection-status confound -- we
are asking "given you were infected, does presenting allele X predict this clonotype", which is
the actual definition of HLA restriction.

Also reports what VDJdb says each significant hit is restricted BY, so the stratum result and
the database annotation are compared on the same clonotypes rather than by eyeball. 2026-07-16.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19")
RES = Path("/projects/biomarkers/results")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")
MIN_READS = 10_000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRA")
    ap.add_argument("--top", type=int, default=40)
    args = ap.parse_args()
    print(f"=== HLA restriction of our top {args.chain} hits (within COVID+ donors only) ===",
          flush=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    meta = meta.with_columns(pl.col("reads").cast(pl.Int64, strict=False)).filter(
        pl.col("reads") >= MIN_READS)
    m = meta.filter(pl.col("locus") == args.chain)
    covid = (m.filter(pl.col("COVID_status") == "COVID")
             .select(pl.col("sample_id").cast(pl.String)).unique())

    hcols = [c for c in m.columns if c.startswith("HLA-")]
    hla = (m.select(pl.col("sample_id").cast(pl.String), *[pl.col(c) for c in hcols])
           .unpivot(index="sample_id", on=hcols, variable_name="slot", value_name="allele")
           .filter(pl.col("allele").is_not_null() & (pl.col("allele") != "")
                   & (pl.col("allele") != "nan"))
           .join(covid, on="sample_id", how="inner")
           .unique(["sample_id", "allele"]))
    typed = sorted(set(hla["sample_id"].to_list()))
    print(f"  HLA-typed COVID+ donors: {len(typed):,}", flush=True)

    bm = pl.read_parquet(RES / f"bm_{args.chain}_aa_v.parquet").sort("p_donor").head(args.top)
    cands = bm.select("cand", "cand_v")
    print(f"  testing top {cands.height} hits", flush=True)

    frames = []
    for r in m.iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{args.chain}.*"))
        if hits and str(r["sample_id"]) in set(typed):
            frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                          .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    cohort = pl.concat(frames, how="vertical_relaxed")
    rows = (cohort.lazy()
            .filter(pl.col(S.JUNCTION_AA).is_not_null()
                    & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
            .select(S.JUNCTION_AA, "v_gene", "sample_id").unique()
            .collect(engine="streaming"))

    import vdjmatch.cluster as vc
    universe = rows[S.JUNCTION_AA].unique().sort().to_list()
    cq = cands["cand"].to_list()
    pairs = vc.overlap(cq, universe, scope="1,0,0,1", threads=0)
    qmap = pl.DataFrame({"a_idx": np.arange(len(cq), dtype=np.int64), "cand": cq,
                         "cand_v": cands["cand_v"].to_list()})
    umap = pl.DataFrame({"b_idx": np.arange(len(universe), dtype=np.int64),
                         S.JUNCTION_AA: universe})
    carry = (pairs.join(qmap, on="a_idx").join(umap, on="b_idx")
             .join(rows, left_on=[S.JUNCTION_AA, "cand_v"],
                   right_on=[S.JUNCTION_AA, "v_gene"], how="inner")
             .select("cand", "cand_v", "sample_id").unique())

    # alleles with enough carriers among COVID+ typed donors
    ac = hla.group_by("allele").agg(pl.len().alias("n")).filter(pl.col("n") >= 30)
    alleles = ac["allele"].to_list()
    print(f"  alleles with >=30 COVID+ carriers: {len(alleles)}", flush=True)
    who = {al: set(hla.filter(pl.col("allele") == al)["sample_id"].to_list()) for al in alleles}
    ntot = len(typed)

    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    vd = (db.filter((db["gene"] == args.chain) & (db["species"] == "HomoSapiens")
                    & db["antigen.species"].str.contains("SARS-CoV-2"))
          .select("cdr3", "antigen.epitope", "mhc.a").unique())
    vmap = {r["cdr3"]: (r["antigen.epitope"], r["mhc.a"]) for r in vd.iter_rows(named=True)}

    print(f"\n  {'clonotype':22s} {'V':10s} {'nCOVID':>6s}  best HLA (Fisher within COVID+)"
          f"          VDJdb says", flush=True)
    hits = []
    for row in cands.iter_rows(named=True):
        c, v = row["cand"], row["cand_v"]
        car = set(carry.filter((pl.col("cand") == c) & (pl.col("cand_v") == v))["sample_id"]
                  .to_list()) & set(typed)
        if len(car) < 5:
            continue
        best = None
        for al in alleles:
            w = who[al]
            a = len(car & w); b = len(car) - a
            cc = len(w) - a; d = ntot - a - b - cc
            orr, p = fisher_exact([[a, b], [cc, d]], alternative="greater")
            if best is None or p < best[2]:
                best = (al, orr, p, a, len(car))
        vinfo = vmap.get(c, ("", ""))
        hits.append({"cand": c, "v": v, "n": len(car), "allele": best[0], "or": best[1],
                     "p": best[2], "n_carriers_with_allele": best[3],
                     "vdjdb_epitope": vinfo[0], "vdjdb_mhc": vinfo[1]})
        print(f"  {c[:22]:22s} {v[:10]:10s} {len(car):6d}  {best[0]:12s} OR={best[1]:6.2f} "
              f"p={best[2]:8.1e}   {vinfo[1]:12s} {vinfo[0]}", flush=True)

    if hits:
        h = pl.DataFrame(hits)
        h = h.with_columns(pl.Series("q", stats.fdr_bh(h["p"].to_numpy())))
        sig = h.filter(pl.col("q") < 0.05)
        print(f"\n  hits with a SIGNIFICANT HLA association (q<0.05): {sig.height}/{h.height}",
              flush=True)
        if sig.height:
            print(sig.group_by("allele").agg(pl.len().alias("n_clonotypes"))
                  .sort("n_clonotypes", descending=True), flush=True)
        known = h.filter(pl.col("vdjdb_mhc") != "")
        if known.height:
            print(f"\n  where VDJdb knows the restriction, do we agree?", flush=True)
            print(known.select("cand", "allele", "p", "vdjdb_mhc", "vdjdb_epitope"), flush=True)
        h.write_parquet(RES / f"motif_hla_{args.chain}.parquet")


if __name__ == "__main__":
    main()
