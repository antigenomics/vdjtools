"""A1 — COVID association with the paper's ACTUAL 1-mismatch semantics.

Vlasova 2026 Methods p.5-6: *"we modified our method to allow a single mismatch in the CDR3
amino acid sequence WHEN SEARCHING FOR A CLONOTYPE OF INTEREST IN THE REPERTOIRE."*

That is a fuzzy **SEARCH**, not clustering:

    incidence(c) = # donors whose repertoire contains ANY CDR3aa within Hamming-1 of c

Each candidate keeps its identity and GAINS incidence, which is what buys the power that takes
TCRb from "no significant hits" (their exact result, and ours: 1 hit at q<0.05) to 567 at q=0.01.

vdjtools' `match="1mm"` does the OTHER thing — `metaclonotypes()` union-find CLUSTERING, which
MERGES candidates into groups and tests the group. On this cohort that collapsed 1.92M candidates
into 179k and diluted the signal: the strongest exact hit (CASSSAHTYTEAFF, 41 COVID / 1 control,
OR=18.4, p=2.3e-8) vanished into a 35-member metaclonotype at OR=2.9. Grouping is a legitimate
operation (it is what the paper does LATER, to build metaclonotypes from the biomarker list) but
it is not the biomarker search.

Unit is still the DONOR (Fisher on presence/absence); candidates are still the paper's
">=2 unique nt variants OR >=2 samples" rearrangement filter. 2026-07-16.
"""
from __future__ import annotations

import argparse
import resource
import time
from pathlib import Path

import numpy as np
import polars as pl

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19")


def _rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def t(t0, m):
    print(f"  [{m}] {time.perf_counter()-t0:7.1f}s  peakRSS {_rss():.1f} GB", flush=True)


def load(meta, chain):
    frames, miss = [], 0
    for r in meta.filter(pl.col("locus") == chain).iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{chain}.*"))
        if not hits:
            miss += 1
            continue
        frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                      .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    if miss:
        print(f"  {miss} rows without a file", flush=True)
    return pl.concat(frames, how="vertical_relaxed").lazy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    ap.add_argument("--scope", default="1,0,0,1", help="1 substitution, no indels")
    ap.add_argument("--min-samples", type=int, default=2)
    ap.add_argument("--max-candidates", type=int, default=0, help="0 = all")
    ap.add_argument("--out", default="/projects/biomarkers/results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    design = (meta.filter(pl.col("locus") == args.chain)
              .select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                      pl.when(pl.col("COVID_status") == "COVID").then(True)
                        .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                        .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))
    print(f"=== A1-fuzzy covid19 {args.chain} (paper's 1mm SEARCH, scope={args.scope}) ===",
          flush=True)

    t0 = time.perf_counter()
    cohort = load(meta, args.chain).collect().lazy()
    t(t0, "ingest")

    # (cdr3aa, donor) presence, restricted to labelled+observed donors — the association universe.
    t0 = time.perf_counter()
    inc = (cohort.filter(pl.col(S.JUNCTION_AA).is_not_null()
                         & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
           .select(S.JUNCTION_AA, "sample_id").unique()
           .join(design.lazy(), on="sample_id", how="inner")
           .collect(engine="streaming"))
    n_pos = int(design.join(inc.select("sample_id").unique(), on="sample_id",
                            how="semi")["_pos"].sum())
    seen = design.join(inc.select("sample_id").unique(), on="sample_id", how="semi")
    n_pos, n_neg = int(seen["_pos"].sum()), seen.height - int(seen["_pos"].sum())
    t(t0, "incidence table")
    print(f"  ANALYSED arms: {n_pos} COVID / {n_neg} control; "
          f"(cdr3aa,donor) rows: {inc.height}", flush=True)

    # Paper's candidate filter: >=2 unique nt variants OR >=2 samples. The REARRANGEMENT
    # (unique nt row) is the unit -- not reads.
    t0 = time.perf_counter()
    cand = (cohort.filter(pl.col(S.JUNCTION_AA).is_not_null()
                          & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .group_by(S.JUNCTION_AA)
            .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_nt_variants"),
                 pl.col("sample_id").n_unique().alias("exact_incidence"))
            .filter((pl.col("n_nt_variants") >= 2)
                    | (pl.col("exact_incidence") >= args.min_samples))
            .collect(engine="streaming")
            .sort([S.JUNCTION_AA]))                       # deterministic
    if args.max_candidates:
        cand = cand.head(args.max_candidates)
    t(t0, "candidates")
    print(f"  candidates: {cand.height}", flush=True)

    # THE FUZZY SEARCH: candidate -> every cohort CDR3aa within Hamming-1 (incl. itself).
    import vdjmatch.cluster as vc
    universe = inc[S.JUNCTION_AA].unique().sort().to_list()
    cq = cand[S.JUNCTION_AA].to_list()
    t0 = time.perf_counter()
    pairs = vc.overlap(cq, universe, scope=args.scope, threads=0)
    t(t0, "1mm search")
    print(f"  neighbour pairs: {pairs.height}  ({pairs.height/max(len(cq),1):.1f}/candidate)",
          flush=True)

    # Fuzzy incidence: donors carrying ANY neighbour of the candidate.
    t0 = time.perf_counter()
    qmap = pl.DataFrame({"a_idx": np.arange(len(cq), dtype=np.int64), "cand": cq})
    umap = pl.DataFrame({"b_idx": np.arange(len(universe), dtype=np.int64),
                         S.JUNCTION_AA: universe})
    hit = (pairs.join(qmap, on="a_idx", how="inner").join(umap, on="b_idx", how="inner")
           .join(inc, on=S.JUNCTION_AA, how="inner")
           .group_by("cand")
           .agg(pl.col("sample_id").n_unique().alias("present"),
                pl.col("sample_id").filter(pl.col("_pos")).n_unique().alias("a")))
    t(t0, "fuzzy incidence")

    a = hit["a"].to_numpy().astype(np.int64)
    present = hit["present"].to_numpy().astype(np.int64)
    b = present - a
    c, d = n_pos - a, n_neg - b
    p = stats.fisher_p(a, b, c, d, alternative="greater")
    res = hit.with_columns(
        pl.Series("incidence", present), pl.Series("n_pos_present", a),
        pl.Series("n_neg_present", b),
        pl.Series("odds_ratio", stats.odds_ratio(a, b, c, d)),
        pl.Series("log2_or", np.log2(stats.odds_ratio(a, b, c, d))),
        pl.Series("p_value", p), pl.Series("q_value", stats.fdr_bh(p))).sort("p_value")

    print(f"\nfeatures tested : {res.height}", flush=True)
    for lab, n in (("q<0.01 (paper's final)", res.filter(pl.col("q_value") < 0.01).height),
                   ("q<0.05 (paper's explor.)", res.filter(pl.col("q_value") < 0.05).height),
                   ("P<1e-4 (Emerson nominal)", res.filter(pl.col("p_value") < 1e-4).height)):
        print(f"  {lab:26s}: {n}", flush=True)
    print("\npaper (TCRb): 567 at q=0.01;  (TCRa): 4,393", flush=True)
    print("\ntop 12:", flush=True)
    print(res.head(12).select("cand", "incidence", "n_pos_present", "n_neg_present",
                              "odds_ratio", "p_value", "q_value"), flush=True)

    res.write_parquet(out / f"a1_fuzzy_{args.chain}.parquet")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = res["log2_or"].to_numpy()
    y = np.clip(-np.log10(np.clip(res["p_value"].to_numpy(), 1e-320, 1)), 0, 320)
    sig = res["q_value"].to_numpy() < 0.01
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    ax.scatter(x[~sig], y[~sig], s=5, c="#c8ccd4", alpha=.4, lw=0, label="ns")
    ax.scatter(x[sig], y[sig], s=13, c="#D55E00", alpha=.8, lw=0,
               label=f"q<0.01 (n={int(sig.sum())})")
    ax.axhline(4, ls=":", c="#0072B2", lw=.8, label="P<1e-4 (Emerson)")
    ax.axvline(0, ls="-", c="k", lw=.6, alpha=.3)
    ax.set(xlabel="log2 odds ratio  (enriched in COVID →)", ylabel="−log10 p",
           title=f"COVID-associated {args.chain} — CDR3aa±1mm search ({n_pos}+/{n_neg}−)")
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"a1_fuzzy_volcano_{args.chain}.png", dpi=140)
    print(f"\n-> {out}", flush=True)


if __name__ == "__main__":
    main()
