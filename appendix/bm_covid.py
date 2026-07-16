"""A1 — COVID-19 association on the FULL FMBA cohort, the paper's method.

Vlasova 2026 (Genome Medicine 18:20), Methods pp.5-6, reproduced faithfully:

  * unit      : the DONOR. Fisher exact on presence/absence. Never weight by reads —
                this is gDNA multiplex with no UMIs, so `count` is reads and cells are
                unrecoverable; weighting the table by reads is pseudoreplication
                (Hurlbert 1984), and Emerson 2017 tested template-weighted abundance
                head-to-head and it lost.
  * candidates: "rearranged at least twice — i.e. supported by at least two unique
                nucleotide variants OR found in at least two samples". The REARRANGEMENT
                (a unique nt row) is the countable unit: one CDR3aa reached by 3 distinct
                nt variants is 3 independent recombination events.
  * key       : CDR3aa ALONE, Hamming <= 1. No V, no J. The paper is explicit that EXACT
                matching yields ZERO significant TCRb hits — so exact is run here as a
                negative control, not as the primary.
  * threshold : BH q=0.01 (paper's final) AND nominal P<1e-4 (Emerson's field-standard
                line, whose FDR is estimated by permuting labels rather than assumed).

Data: /projects/biomarkers/raw/covid19 (the FULL 1,258-donor HF repo — NOT
/projects/fmba_covid/data, which is a 573-donor subset with 36 controls). Sample ids come
from the sheet's own file_name->sample_id map, never from the filename.  2026-07-16.
"""
from __future__ import annotations

import argparse
import resource
import time
from pathlib import Path

import numpy as np
import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S
from vdjtools.biomarker import association

ROOT = Path("/projects/biomarkers/raw/covid19")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")


def _rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def t(t0, msg):
    print(f"  [{msg}] {time.perf_counter()-t0:7.1f}s  peakRSS {_rss():.1f} GB", flush=True)


def load(meta: pl.DataFrame, chain: str) -> pl.LazyFrame:
    """Cohort for one chain, sample ids driven by the metadata's file_name -> sample_id map."""
    m = meta.filter(pl.col("locus") == chain)
    frames, miss = [], 0
    for r in m.iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{chain}.*"))
        if not hits:
            miss += 1
            continue
        frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                      .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    if miss:
        print(f"  {miss} {chain} metadata rows without a file", flush=True)
    return pl.concat(frames, how="vertical_relaxed").lazy()


def rearrangement_candidates(cohort: pl.LazyFrame, min_rearr=2, min_samples=2) -> pl.DataFrame:
    """The paper's filter: >=2 unique nt variants OR present in >=2 samples.

    `n_rearrangements` counts DISTINCT (junction_nt, sample_id) — a unique nt row in a
    repertoire is one recombination event; the same nt in two donors is two events
    (convergent recombination). This is the unit Pgen is defined on, and it is NOT reads.
    """
    return (cohort.filter(pl.col(S.JUNCTION_AA).is_not_null()
                          & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .group_by(S.JUNCTION_AA)
            .agg(pl.struct(S.JUNCTION_NT, "sample_id").n_unique().alias("n_rearrangements"),
                 pl.col(S.JUNCTION_NT).n_unique().alias("n_nt_variants"),
                 pl.col("sample_id").n_unique().alias("incidence"))
            .filter((pl.col("n_nt_variants") >= min_rearr) | (pl.col("incidence") >= min_samples))
            .collect(engine="streaming"))


def volcano(res, out: Path, title: str, q=0.01):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = res["log2_or"].to_numpy()
    y = np.clip(-np.log10(np.clip(res["p_value"].to_numpy(), 1e-320, 1)), 0, 320)
    sig = res["q_value"].to_numpy() < q
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    ax.scatter(x[~sig], y[~sig], s=5, c="#c8ccd4", alpha=.4, lw=0, label="ns")
    ax.scatter(x[sig], y[sig], s=13, c="#D55E00", alpha=.8, lw=0, label=f"q<{q} (n={int(sig.sum())})")
    ax.axhline(4, ls=":", c="#0072B2", lw=.8, label="P<1e-4 (Emerson)")
    ax.axvline(0, ls="-", c="k", lw=.6, alpha=.3)
    ax.set(xlabel="log2 odds ratio  (enriched in COVID →)", ylabel="−log10 p", title=title)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"  volcano -> {out}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    ap.add_argument("--match", default="1mm", choices=["exact", "1mm"])
    ap.add_argument("--with-v", action="store_true", help="key = CDR3aa + V (paper's secondary)")
    ap.add_argument("--min-samples", type=int, default=2)
    ap.add_argument("--out", default="/projects/biomarkers/results")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    key = (S.JUNCTION_AA, S.V_CALL) if args.with_v else (S.JUNCTION_AA,)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    design = (meta.filter(pl.col("locus") == args.chain)
              .select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                      pl.when(pl.col("COVID_status") == "COVID").then(True)
                        .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                        .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))
    npos = int(design["_pos"].sum())
    print(f"=== A1 covid19 {args.chain}  key={key}  match={args.match} ===", flush=True)
    print(f"design: {npos} COVID / {design.height-npos} control (healthy+precovid)", flush=True)

    t0 = time.perf_counter()
    cohort = load(meta, args.chain).collect().lazy()
    t(t0, "ingest")

    t0 = time.perf_counter()
    cand = rearrangement_candidates(cohort, min_samples=args.min_samples)
    t(t0, "candidates")
    print(f"  candidate CDR3aa (>=2 nt variants OR >={args.min_samples} samples): {cand.height}",
          flush=True)
    print(f"  rearrangements: total {int(cand['n_rearrangements'].sum())}, "
          f"median/clonotype {cand['n_rearrangements'].median()}", flush=True)

    t0 = time.perf_counter()
    res = association(cohort, design, test="fisher", key=key, match=args.match,
                      min_incidence=args.min_samples,
                      candidates=cand.select(list(key)) if set(key) <= set(cand.columns)
                      else cand.select(S.JUNCTION_AA),
                      alternative="greater")
    t(t0, f"association ({args.match})")

    n_pos, n_neg = int(res["n_pos"][0]), int(res["n_neg"][0])
    q01 = res.filter(pl.col("q_value") < 0.01).height
    q05 = res.filter(pl.col("q_value") < 0.05).height
    p4 = res.filter(pl.col("p_value") < 1e-4).height
    print(f"\nANALYSED arms: {n_pos} COVID / {n_neg} control", flush=True)
    print(f"features tested : {res.height}", flush=True)
    print(f"  q<0.01 (paper's final)   : {q01}", flush=True)
    print(f"  q<0.05 (paper's explor.) : {q05}", flush=True)
    print(f"  P<1e-4 (Emerson nominal) : {p4}", flush=True)
    print("\npaper: 4,393 alpha / 567 beta at q=0.01 (CDR3aa+-1mm)", flush=True)
    print("\ntop 12:", flush=True)
    cols = [c for c in ("junction_aa", "v_call", "n_members", "incidence", "n_pos_present",
                        "n_neg_present", "odds_ratio", "p_value", "q_value") if c in res.columns]
    print(res.sort("p_value").head(12).select(cols), flush=True)

    tag = f"{args.chain}_{args.match}{'_V' if args.with_v else ''}"
    res.write_parquet(outdir / f"a1_covid_{tag}.parquet")
    volcano(res, outdir / f"a1_covid_volcano_{tag}.png",
            f"COVID-associated {args.chain} ({args.match}, CDR3aa{'+V' if args.with_v else ''}; "
            f"{n_pos}+/{n_neg}−)")
    print(f"\n-> {outdir}", flush=True)


if __name__ == "__main__":
    main()
