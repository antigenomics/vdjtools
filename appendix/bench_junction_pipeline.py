#!/usr/bin/env python
"""Benchmark the junction pipeline against real nucleotide truth -- 2026-09-30.

`vdjtools.model.annotate_junctions` answers four questions off an amino-acid junction alone. Three
of them have an external answer in `isalgo/airr_control`, whose rows are REAL nucleotide
rearrangements with V/D/J called from the sequence and `VEnd`/`DStart`/`DEnd`/`JStart` beside them:

* the inferred nucleotide junction -- against the observed one, exactly;
* the D gene by nucleotide alignment (stage 3) -- against the nucleotide caller's;
* the D gene by the length-and-prior posterior (`posterior_d_batch`) -- against the same truth,
  which is what says whether the probabilistic route is still needed now the nucleotide one exists.

The input to the pipeline is only `(cdr3aa, v gene, j gene)`: the nucleotides are hidden from it and
used only for scoring.

    python appendix/bench_junction_pipeline.py --locus trb --n 5000

Nothing is cached: the sample is redrawn from the source table on every run.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import polars as pl

SRC = Path.home() / "hf" / "airr_control"


def sample(locus: str, n: int, min_donors: int, seed: int) -> pl.DataFrame:
    """`n` rearrangements with a called D, seen in >= `min_donors` donors, one row per junction."""
    src = SRC / f"human.{locus}.ntvj.vdjtools.tsv.gz"
    lf = pl.scan_csv(src, separator="\t", has_header=True,
                     schema_overrides={"v": pl.String, "j": pl.String, "d": pl.String,
                                       "cdr3nt": pl.String, "cdr3aa": pl.String,
                                       "count": pl.Int64, "incidence": pl.Int32,
                                       "VEnd": pl.Int32, "DStart": pl.Int32, "DEnd": pl.Int32,
                                       "JStart": pl.Int32})
    df = (lf.select("cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart",
                    "count", "incidence")
          .filter(pl.col("incidence") >= min_donors)
          .filter(pl.col("cdr3aa").str.starts_with("C") & pl.col("d").is_not_null()
                  & (pl.col("DStart") >= 0) & (pl.col("DEnd") >= 0))
          # One row per junction, the most abundant reading of it. Sorted on the junction too, so
          # `unique(keep="first")` is deterministic and a recorded number is reproducible -- a
          # streaming collect does not promise input order.
          .sort(["count", "cdr3aa"], descending=[True, False])
          .unique(subset=["cdr3aa"], keep="first")
          .collect(engine="streaming"))
    # Sort AFTER the collect: a streaming engine does not promise to hand rows back in plan order, so
    # sorting inside the plan still let `sample` draw a different 4,000 on every run -- which showed
    # up as +-2 pp of "difference" between two builds of the same code.
    return df.sort("cdr3aa").sample(n=min(n, df.height), seed=seed)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--locus", default="trb")
    ap.add_argument("--n", type=int, default=5000)
    ap.add_argument("--min-donors", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--n-best", type=int, default=8)
    args = ap.parse_args()

    from vdjtools.model import annotate_junctions

    t = sample(args.locus, args.n, args.min_donors, args.seed)
    print(f"truth: {t.height:,} human {args.locus.upper()} rearrangements with a nucleotide D call")

    cdr3, v, j = t["cdr3aa"].to_list(), t["v"].to_list(), t["j"].to_list()
    annotate_junctions(cdr3[:10], v[:10], j[:10], species="human")     # warm the caches
    t0 = time.perf_counter()
    got = annotate_junctions(cdr3, v, j, species="human", n_best=args.n_best)
    dt = time.perf_counter() - t0
    print(f"pipeline: {dt:.2f}s = {dt / t.height * 1e6:.0f} us/junction "
          f"-> {dt / t.height * 190_000:.0f}s for a 190k-key corpus")

    df = pl.concat([t, got.select(pl.exclude("cdr3_aa"))], how="horizontal")
    gene = lambda c: pl.col(c).str.split("*").list.first()        # noqa: E731
    d_true = gene("d")

    # ---- the nucleotide guess
    have_nt = df.filter(pl.col("cdr3_nt").is_not_null())
    exact = have_nt.filter(pl.col("cdr3_nt") == pl.col("cdr3nt"))
    print(f"\nnucleotide junction: inferred for {have_nt.height:,}/{df.height:,}, "
          f"EXACT on {exact.height:,} = {exact.height / max(have_nt.height, 1):.2%}")
    if have_nt.height:
        ham = (have_nt.with_columns(
            pl.struct("cdr3_nt", "cdr3nt").map_elements(
                lambda s: sum(a != b for a, b in zip(s["cdr3_nt"], s["cdr3nt"])),
                return_dtype=pl.Int64).alias("mismatches")))
        print(f"  mismatching nucleotides per junction: median "
              f"{ham['mismatches'].median():.0f}, mean {ham['mismatches'].mean():.2f} "
              f"of {have_nt['cdr3nt'].str.len_chars().mean():.1f} nt")

    # ---- the two D answers, against the same truth
    for label, col in (("alignment on the inferred nt (stage 3)", "d_call"),
                       ("length-and-prior posterior", "d_posterior_call")):
        sub = df.filter(pl.col(col).is_not_null())
        ok = sub.filter(gene(col).str.contains(d_true, literal=True)
                        if False else gene(col) == d_true)
        print(f"\nD gene, {label}: called {sub.height:,}/{df.height:,} "
              f"({sub.height / df.height:.1%}), correct {ok.height:,} = "
              f"{ok.height / max(sub.height, 1):.2%} of called, "
              f"{ok.height / df.height:.2%} of all")

    both = df.filter(pl.col("d_call").is_not_null() & pl.col("d_posterior_call").is_not_null())
    if both.height:
        a = both.filter(gene("d_call") == d_true).height
        p = both.filter(gene("d_posterior_call") == d_true).height
        agree = both.filter(gene("d_call") == gene("d_posterior_call")).height
        print(f"\nhead to head on the {both.height:,} where both answer: alignment {a} correct, "
              f"posterior {p} correct, they agree on {agree} ({agree / both.height:.2%})")
        # what does the posterior add where the alignment declines?
        only_p = df.filter(pl.col("d_call").is_null() & pl.col("d_posterior_call").is_not_null())
        if only_p.height:
            print(f"posterior-only rows (alignment declined): {only_p.height:,}, correct "
                  f"{only_p.filter(gene('d_posterior_call') == d_true).height:,} = "
                  f"{only_p.filter(gene('d_posterior_call') == d_true).height / only_p.height:.2%}")

    # ---- the combination the pipeline actually reports
    sub = df.filter(pl.col("d_best").is_not_null())
    ok = sub.filter(pl.col("d_best") == d_true)
    print(f"\nD gene, `d_best` (alignment, posterior as fallback): called {sub.height:,}/"
          f"{df.height:,}, correct {ok.height:,} = {ok.height / df.height:.2%} of all")
    for src in ("alignment", "posterior"):
        part = sub.filter(pl.col("d_best_source") == src)
        if part.height:
            print(f"  from the {src:<9}: {part.height:,} rows, correct "
                  f"{part.filter(pl.col('d_best') == d_true).height / part.height:.2%}")

    # ---- where the D sits, in nucleotides (truth coordinates are 0-based half-open)
    pos = df.filter(pl.col("d_start_nt").is_not_null() & (gene("d_call") == d_true))
    if pos.height:
        pos = pos.with_columns((pl.col("d_start_nt") - 1 - pl.col("DStart")).alias("ds"),
                               (pl.col("d_end_nt") - 1 - pl.col("DEnd")).alias("de"))
        print(f"\nD position on the {pos.height:,} correctly called: d_start exact "
              f"{pos.filter(pl.col('ds') == 0).height / pos.height:.2%}, within 1 nt "
              f"{pos.filter(pl.col('ds').abs() <= 1).height / pos.height:.2%}; d_end exact "
              f"{pos.filter(pl.col('de') == 0).height / pos.height:.2%}, within 1 nt "
              f"{pos.filter(pl.col('de').abs() <= 1).height / pos.height:.2%}")


if __name__ == "__main__":
    main()
