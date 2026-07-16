"""Diagnostic: what does OUR pipeline say about the PAPER'S published biomarkers?

We get 7 alpha / 0 beta at q=0.01; the paper reports 4,393 / 567. Guessing at the gap is
cheap and wrong. This asks the data directly, using their own published list
(covid_associated_clonotypes.csv: 4,393 alpha / 567 beta, has_covid_association flag):

  1. Are their biomarkers even present in our cohort?  (if not -> data/QC mismatch)
  2. What incidence do they have in our arms?           (if low -> our cohort differs)
  3. What p/q do WE assign them?                        (if small p but big q -> multiple
                                                          testing / candidate-set size)
  4. Are they enriched among OUR low-p features?        (rank concordance -> same signal,
                                                          different threshold)

That distinguishes "we can't see the signal" from "we see it but call it non-significant",
which are entirely different bugs. 2026-07-16.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path("/projects/biomarkers/raw/covid19")
RES = Path("/projects/biomarkers/results")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    args = ap.parse_args()
    chain_lab = "beta" if args.chain == "TRB" else "alpha"

    pub = pl.read_csv(ROOT / "covid_associated_clonotypes.csv", infer_schema_length=0)
    pub_c = pub.filter((pl.col("chain") == chain_lab))
    pub_true = pub_c.filter(pl.col("has_covid_association") == "True")
    print(f"=== published list: {pub.height} rows; {chain_lab}: {pub_c.height} "
          f"({pub_true.height} has_covid_association=True) ===", flush=True)

    res = pl.read_parquet(RES / f"a1_fuzzy_{args.chain}.parquet")
    print(f"our fuzzy A1 {args.chain}: {res.height} tested", flush=True)

    ours = set(res["cand"].to_list())
    theirs = set(pub_c["cdr3"].to_list())
    theirs_t = set(pub_true["cdr3"].to_list())
    print(f"\n1. PRESENCE: {len(theirs & ours)} of {len(theirs)} published {chain_lab} CDR3s are "
          f"in our tested set ({len(theirs & ours)/max(len(theirs),1):.1%})", flush=True)
    print(f"   (has_covid_association=True: {len(theirs_t & ours)} of {len(theirs_t)})", flush=True)

    hit = res.filter(pl.col("cand").is_in(list(theirs_t)))
    if not hit.height:
        print("   none of their biomarkers are in our tested set — data mismatch, stop here.")
        return
    print(f"\n2. INCIDENCE of their biomarkers in OUR arms (n={hit.height}):", flush=True)
    for c in ("incidence", "n_pos_present", "n_neg_present", "odds_ratio"):
        v = hit[c].to_numpy()
        print(f"   {c:15s} median {np.median(v):8.2f}  p90 {np.percentile(v,90):8.2f} "
              f" max {v.max():8.2f}", flush=True)

    print("\n3. OUR p/q for THEIR biomarkers:", flush=True)
    for thr, lab in ((1e-4, "p<1e-4"), (0.01, "q<0.01"), (0.05, "q<0.05")):
        col = "p_value" if lab.startswith("p") else "q_value"
        n = hit.filter(pl.col(col) < thr).height
        print(f"   {lab:8s}: {n:5d} of {hit.height} ({n/hit.height:.1%})", flush=True)
    print(f"   median p {hit['p_value'].median():.3g}   min p {hit['p_value'].min():.3g}",
          flush=True)

    # 4. rank concordance: are their biomarkers enriched at the top of OUR ranking?
    r = res.sort("p_value").with_row_index("rank")
    rk = r.filter(pl.col("cand").is_in(list(theirs_t)))["rank"].to_numpy()
    from scipy.stats import fisher_exact
    for top_n in (1000, 10000, 100000):
        a = int((rk < top_n).sum())
        b = len(rk) - a
        c = top_n - a
        d = (res.height - top_n) - b
        orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        print(f"\n4. top-{top_n:>6} of our ranking holds {a} of their {len(rk)} biomarkers "
              f"-> OR={orr:.1f} p={p:.2e}", flush=True)


if __name__ == "__main__":
    main()
