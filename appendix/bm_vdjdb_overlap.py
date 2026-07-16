"""V1 — VDJdb as the oracle for our COVID association ranking.

The published biomarker list is one lab's output on one cohort. VDJdb is independent: if our
ranking carries real SARS-CoV-2 signal, VDJdb's SARS-CoV-2 CDR3s should enrich at the top of
it -- and VDJdb's CMV CDR3s should NOT (the built-in negative control that makes the test
mean something; a ranking that enriches both is enriching "public/annotated", not "COVID").

Reported per (chain x match) for our exact and fuzzy runs:
  * enrichment of VDJdb SARS-CoV-2 vs CMV in top-K of our ranking (Fisher, one-tailed)
  * the same for THEIR published biomarker list, so the two oracles are on one scale

Matching is exact CDR3aa (VDJdb <-> cohort). Species human. 2026-07-16.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import fisher_exact

RES = Path("/projects/biomarkers/results")
ROOT = Path("/projects/biomarkers/raw/covid19")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")


def vdjdb_cdr3(chain: str, pattern: str) -> set[str]:
    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    gene = "TRA" if chain == "TRA" else "TRB"
    f = ((pl.col("gene") == gene)
         & (pl.col("species") == "HomoSapiens")
         & pl.col("antigen.species").str.contains(pattern))
    return set(db.filter(f)["cdr3"].unique().to_list())


def enrich(rank_col: pl.DataFrame, key: str, truth: set[str], label: str, total: int):
    r = rank_col.sort("p_value").with_row_index("rank")
    inset = r.filter(pl.col(key).is_in(list(truth)))
    n_truth = inset.height
    if not n_truth:
        print(f"   {label:16s} 0 of the oracle's CDR3s are in our tested set", flush=True)
        return
    rk = inset["rank"].to_numpy()
    out = []
    for K in (100, 1000, 10000, 100000):
        a = int((rk < K).sum())
        b = n_truth - a
        c = K - a
        d = (total - K) - b
        orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        out.append(f"top-{K:<6} {a:5d}/{n_truth:<6} OR={orr:6.2f} p={p:9.2e}")
    print(f"   {label:16s} (n={n_truth})", flush=True)
    for o in out:
        print(f"       {o}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    args = ap.parse_args()
    chain_lab = "beta" if args.chain == "TRB" else "alpha"

    cov = vdjdb_cdr3(args.chain, "SARS-CoV-2")
    cmv = vdjdb_cdr3(args.chain, "CMV")
    print(f"=== VDJdb oracle {args.chain}: SARS-CoV-2 {len(cov)} | CMV {len(cmv)} (control) ===",
          flush=True)

    pub = pl.read_csv(ROOT / "covid_associated_clonotypes.csv", infer_schema_length=0)
    theirs = set(pub.filter((pl.col("chain") == chain_lab)
                            & (pl.col("has_covid_association") == "True"))["cdr3"].to_list())
    print(f"    published biomarkers {args.chain}: {len(theirs)}", flush=True)
    print(f"    VDJdb SARS-CoV-2 ∩ published: {len(cov & theirs)}", flush=True)

    for tag, key in ((f"a1_covid_{args.chain}_exact", "junction_aa"),
                     (f"a1_fuzzy_{args.chain}", "cand")):
        p = RES / f"{tag}.parquet"
        if not p.exists():
            print(f"\n### {tag}: MISSING", flush=True)
            continue
        res = pl.read_parquet(p)
        k = key if key in res.columns else ("cand" if "cand" in res.columns else "junction_aa")
        print(f"\n### {tag}  ({res.height} tested)", flush=True)
        enrich(res, k, cov, "VDJdb SARS-CoV-2", res.height)
        enrich(res, k, cmv, "VDJdb CMV (ctrl)", res.height)
        enrich(res, k, theirs, "published list", res.height)


if __name__ == "__main__":
    main()
