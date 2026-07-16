"""Are OUR biomarkers real? VDJdb SARS-CoV-2 vs CMV, on our own list.

The row-based test buys a large power boost. A boost is only good news if it lands on the
right clonotypes -- an inflated variance would boost noise just as happily. So this asks an
oracle that had no part in the discovery:

  * VDJdb SARS-CoV-2 CDR3s should enrich among our significant biomarkers.
  * VDJdb CMV CDR3s should NOT. This is the control that makes the test mean anything: a
    ranking that enriches both is enriching "public and well-annotated", not "COVID".

Reported for each unit (donor / rows) at its own FDR cut, plus a rank sweep so the two units
are on one scale. The published list is included last as a courtesy check, NOT as the target
-- our list is the deliverable. 2026-07-16.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl
from scipy.stats import fisher_exact

RES = Path("/projects/biomarkers/results")
ROOT = Path("/projects/biomarkers/raw/covid19")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")


def oracle(chain: str, pattern: str) -> set[str]:
    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    return set(db.filter((pl.col("gene") == chain)
                         & (pl.col("species") == "HomoSapiens")
                         & pl.col("antigen.species").str.contains(pattern))["cdr3"]
               .unique().to_list())


def enrich_set(sig: set[str], rest: set[str], truth: set[str]):
    a = len(sig & truth)
    b = len(sig) - a
    c = len(rest & truth)
    d = len(rest) - c
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, len(sig), orr, p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB")
    ap.add_argument("--scope", default="aa")
    args = ap.parse_args()
    tag = f"{args.chain}_{args.scope}"
    p = RES / f"bm_{tag}.parquet"
    if not p.exists():
        print(f"{p} MISSING")
        return
    res = pl.read_parquet(p)
    cov = oracle(args.chain, "SARS-CoV-2")
    cmv = oracle(args.chain, "CMV")
    print(f"=== validate bm_{tag}  ({res.height:,} tested) ===", flush=True)
    print(f"    VDJdb {args.chain}: SARS-CoV-2 {len(cov):,} | CMV {len(cmv):,} (control)",
          flush=True)

    allc = set(res["cand"].to_list())
    print(f"    our tested CDR3aa n VDJdb: SARS-CoV-2 {len(allc & cov):,} | "
          f"CMV {len(allc & cmv):,}", flush=True)

    for qcol, unit in (("q_donor", "DONOR (Fisher, incidence)"),
                       ("q_rows", "ROWS  (binomial, rearrangements)")):
        for thr in (0.01, 0.05):
            sig = set(res.filter(pl.col(qcol) < thr)["cand"].to_list())
            if not sig:
                print(f"\n  {unit} q<{thr}: 0 significant", flush=True)
                continue
            rest = allc - sig
            print(f"\n  {unit} q<{thr}: {len(sig):,} significant", flush=True)
            for name, truth in (("SARS-CoV-2", cov), ("CMV (ctrl)", cmv)):
                a, n, orr, pv = enrich_set(sig, rest, truth)
                print(f"      {name:12s} {a:4d}/{n:<7,} hits  OR={orr:6.2f}  p={pv:9.2e}",
                      flush=True)

    # rank sweep -- puts the two units on one scale regardless of their FDR behaviour
    for pcol, unit in (("p_donor", "donor"), ("p_rows", "rows")):
        r = res.sort(pcol).with_row_index("rank")
        print(f"\n  --- rank sweep by {unit} p ---", flush=True)
        for K in (100, 1000, 10000):
            top = set(r.filter(pl.col("rank") < K)["cand"].to_list())
            rest = allc - top
            line = []
            for name, truth in (("CoV2", cov), ("CMV", cmv)):
                a, n, orr, pv = enrich_set(top, rest, truth)
                line.append(f"{name} {a:3d} OR={orr:6.2f} p={pv:8.1e}")
            print(f"      top-{K:<6} {'  |  '.join(line)}", flush=True)

    # courtesy: overlap with the published list (not the target)
    pub = pl.read_csv(ROOT / "covid_associated_clonotypes.csv", infer_schema_length=0)
    lab = "beta" if args.chain == "TRB" else "alpha"
    theirs = set(pub.filter(pl.col("chain") == lab)["cdr3"].to_list())
    for qcol, unit in (("q_donor", "donor"), ("q_rows", "rows")):
        sig = set(res.filter(pl.col(qcol) < 0.01)["cand"].to_list())
        if sig:
            a, n, orr, pv = enrich_set(sig, allc - sig, theirs)
            print(f"\n  [courtesy] published {lab} list vs our {unit} q<0.01: "
                  f"{a}/{n:,} OR={orr:.2f} p={pv:.2e}  (their list n={len(theirs):,})",
                  flush=True)


if __name__ == "__main__":
    main()
