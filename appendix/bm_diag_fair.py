"""Is our 1mm incidence FAIR, and does it track Pgen the way the paper claims?

Three questions, no hand-waving:

Q1 FAIRNESS. `incidence = n_unique(donor)` over the exploded (candidate -> neighbour ->
   donor) join is a UNION, not a sum -- 24 neighbours x 40 donors is NOT 960 donors, it may
   be 40. Verify the union independently: recompute a few candidates' donor sets in plain
   Python from the raw incidence table and assert equality with the polars aggregate. Also
   report how many DISTINCT neighbours actually carry each donor -- if the median is 1, the
   ball is doing nothing but exact matching; if it is large, one hub sequence may dominate.

Q2 EXACT vs FUZZY on THEIR biomarkers. We already know their 256 beta biomarkers get
   OR~0.95 under our fuzzy search. Look the SAME CDR3s up in the exact-match run:
     - exact OR high / p low  -> the fuzzy step destroys the contrast; the ball is the bug.
     - exact OR ~1 too        -> they are not associated in our data at all; the difference
                                 is upstream (batch correction / resampling / QC), not search.
   This is the decisive split and it costs one parquet read.

Q3 PGEN. The paper justifies 1mm by claiming *"the estimate of clonotype population
   frequency computed this way is in perfect agreement with theoretical baseline V(D)J
   rearrangement probability"*. We have EXACT native Hamming-1 ball Pgen, so this is
   directly testable. Detection model (Pogorelyy 2018): a donor with N_i rearrangements
   carries the ball w.p. 1-(1-P_ball)^N_i ~ 1-exp(-N_i*P_ball). So
        E[incidence] = sum_i (1 - exp(-N_i * P_ball))
   Compare observed fuzzy incidence to that. Agreement => our ball is fair and the paper's
   claim holds. Systematic excess => we over-match.

2026-07-16.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19")
RES = Path("/projects/biomarkers/results")


def load(meta, chain):
    frames = []
    for r in meta.filter(pl.col("locus") == chain).iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{chain}.*"))
        if hits:
            frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                          .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    return pl.concat(frames, how="vertical_relaxed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    ap.add_argument("--n-pgen", type=int, default=300, help="biomarkers to Pgen-check")
    args = ap.parse_args()
    chain_lab = "beta" if args.chain == "TRB" else "alpha"

    pub = pl.read_csv(ROOT / "covid_associated_clonotypes.csv", infer_schema_length=0)
    theirs = (pub.filter((pl.col("chain") == chain_lab)
                         & (pl.col("has_covid_association") == "True"))["cdr3"]
              .unique().to_list())
    print(f"=== their {chain_lab} biomarkers (has_covid_association=True): {len(theirs)} ===",
          flush=True)

    # ---------------------------------------------------------------- Q2 (cheap, decisive)
    print("\n########## Q2: THEIR biomarkers under OUR EXACT match ##########", flush=True)
    ex = pl.read_parquet(RES / f"a1_covid_{args.chain}_exact.parquet")
    kcol = "junction_aa" if "junction_aa" in ex.columns else "cand"
    hit_ex = ex.filter(pl.col(kcol).is_in(theirs))
    print(f"exact run: {ex.height} tested; {hit_ex.height} of {len(theirs)} of their "
          f"biomarkers present", flush=True)
    if hit_ex.height:
        for c in ("incidence", "n_pos_present", "n_neg_present", "odds_ratio", "p_value"):
            if c in hit_ex.columns:
                v = hit_ex[c].to_numpy().astype(float)
                print(f"   EXACT {c:15s} median {np.median(v):10.4g}  p90 {np.percentile(v,90):10.4g}"
                      f"  max {v.max():10.4g}", flush=True)
        print(f"   EXACT p<1e-4: {hit_ex.filter(pl.col('p_value') < 1e-4).height} "
              f"of {hit_ex.height}", flush=True)
        # rate in each arm -- the thing that actually decides the OR
        npos, nneg = int(ex["n_pos"][0]), int(ex["n_neg"][0])
        r_pos = hit_ex["n_pos_present"].to_numpy() / npos
        r_neg = hit_ex["n_neg_present"].to_numpy() / nneg
        print(f"   EXACT carrier rate: COVID median {np.median(r_pos):.4f} | "
              f"control median {np.median(r_neg):.4f}   (arms {npos}+/{nneg}-)", flush=True)

    fz = pl.read_parquet(RES / f"a1_fuzzy_{args.chain}.parquet")
    hit_fz = fz.filter(pl.col("cand").is_in(theirs))
    r_pos = hit_fz["n_pos_present"].to_numpy() / 761
    r_neg = hit_fz["n_neg_present"].to_numpy() / 479
    print(f"\n   FUZZY carrier rate: COVID median {np.median(r_pos):.4f} | "
          f"control median {np.median(r_neg):.4f}", flush=True)
    print(f"   FUZZY incidence median {hit_fz['incidence'].median():.1f}; "
          f"OR median {hit_fz['odds_ratio'].median():.3f}", flush=True)

    # ---------------------------------------------------------------- rebuild incidence
    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    design = (meta.filter(pl.col("locus") == args.chain)
              .select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                      pl.when(pl.col("COVID_status") == "COVID").then(True)
                        .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                        .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))
    t0 = time.perf_counter()
    cohort = load(meta, args.chain)
    print(f"\n[ingest {time.perf_counter()-t0:.0f}s]", flush=True)

    inc = (cohort.lazy()
           .filter(pl.col(S.JUNCTION_AA).is_not_null()
                   & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
           .select(S.JUNCTION_AA, "sample_id").unique()
           .join(design.lazy(), on="sample_id", how="inner")
           .collect(engine="streaming"))
    # per-donor rearrangement count N_i = unique nt rows in that repertoire (the Pgen unit)
    depth = (cohort.lazy()
             .filter(pl.col(S.JUNCTION_AA).is_not_null()
                     & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
             .group_by("sample_id").agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
             .collect(engine="streaming")
             .join(design, on="sample_id", how="semi"))
    N = depth["n_rearr"].to_numpy().astype(float)
    print(f"donors {depth.height}; rearrangements/donor: median {np.median(N):,.0f} "
          f"mean {N.mean():,.0f} total {N.sum():,.0f}", flush=True)

    # ---------------------------------------------------------------- Q1 fairness
    print("\n########## Q1: is the fuzzy incidence a FAIR union? ##########", flush=True)
    import vdjmatch.cluster as vc
    universe = inc[S.JUNCTION_AA].unique().sort().to_list()
    probe = [c for c in theirs if c in set(universe)][:40]
    pairs = vc.overlap(probe, universe, scope="1,0,0,1", threads=0)
    umap = pl.DataFrame({"b_idx": np.arange(len(universe), dtype=np.int64),
                         S.JUNCTION_AA: universe})
    qmap = pl.DataFrame({"a_idx": np.arange(len(probe), dtype=np.int64), "cand": probe})
    j = (pairs.join(qmap, on="a_idx").join(umap, on="b_idx").join(inc, on=S.JUNCTION_AA))

    agg = j.group_by("cand").agg(
        pl.col("sample_id").n_unique().alias("union_donors"),
        pl.col(S.JUNCTION_AA).n_unique().alias("n_neighbours"),
        pl.len().alias("n_rows"))
    # independent recomputation in plain python -- must match exactly
    inc_map: dict[str, set] = {}
    for aa, sid in zip(inc[S.JUNCTION_AA].to_list(), inc["sample_id"].to_list()):
        inc_map.setdefault(aa, set()).add(sid)
    nb_map: dict[str, list] = {}
    for cand, aa in zip(j["cand"].to_list(), j[S.JUNCTION_AA].to_list()):
        nb_map.setdefault(cand, []).append(aa)
    bad = 0
    for cand, nbs in nb_map.items():
        ref = set()
        for aa in set(nbs):
            ref |= inc_map.get(aa, set())
        got = int(agg.filter(pl.col("cand") == cand)["union_donors"][0])
        if len(ref) != got:
            bad += 1
            print(f"   MISMATCH {cand}: python union {len(ref)} vs polars {got}", flush=True)
    print(f"   union cross-check on {len(nb_map)} candidates: "
          f"{'ALL MATCH — the aggregate is a fair union' if not bad else f'{bad} MISMATCH'}",
          flush=True)
    print(f"   neighbours/candidate : median {agg['n_neighbours'].median():.0f} "
          f"max {agg['n_neighbours'].max()}", flush=True)
    print(f"   union donors         : median {agg['union_donors'].median():.0f} "
          f"max {agg['union_donors'].max()}", flush=True)
    # how many distinct neighbours does a matched donor actually carry?
    per = (j.group_by(["cand", "sample_id"]).agg(pl.col(S.JUNCTION_AA).n_unique().alias("k")))
    print(f"   distinct neighbours carried per (candidate,donor): median "
          f"{per['k'].median():.0f}  p90 {per['k'].quantile(0.9):.0f}  max {per['k'].max()}",
          flush=True)
    # is one hub sequence carrying the union? report the most public neighbour per candidate
    hub = (j.group_by(["cand", S.JUNCTION_AA]).agg(pl.col("sample_id").n_unique().alias("d"))
           .group_by("cand").agg(pl.col("d").max().alias("top_nb_donors"))
           .join(agg.select("cand", "union_donors"), on="cand")
           .with_columns((pl.col("top_nb_donors") / pl.col("union_donors")).alias("hub_frac")))
    print(f"   single most-public neighbour explains what fraction of the union: median "
          f"{hub['hub_frac'].median():.3f}  p90 {hub['hub_frac'].quantile(0.9):.3f}", flush=True)

    # ---------------------------------------------------------------- Q3 Pgen agreement
    print("\n########## Q3: does fuzzy incidence track the 1mm-ball Pgen? ##########",
          flush=True)
    from vdjtools.model import bundled, native
    m = bundled.load_bundled(args.chain, "olga")
    pm = native.pack(m)
    probe2 = [c for c in theirs if c in set(fz["cand"].to_list())][:args.n_pgen]
    t0 = time.perf_counter()
    ball = native.pgen_aa_batch(pm, probe2, mismatches=1, threads=0)
    pt = native.pgen_aa_batch(pm, probe2, mismatches=0, threads=0)
    print(f"[pgen {len(probe2)} seqs, {time.perf_counter()-t0:.1f}s]", flush=True)
    ball = np.asarray(ball, dtype=float)
    pt = np.asarray(pt, dtype=float)
    # E[incidence] = sum_i 1 - exp(-N_i * P_ball)
    exp_inc = np.array([np.sum(1.0 - np.exp(-N * p)) for p in ball])
    obs = (pl.DataFrame({"cand": probe2})
           .join(fz.select("cand", "incidence"), on="cand", how="left")["incidence"]
           .to_numpy().astype(float))
    ok = np.isfinite(exp_inc) & np.isfinite(obs) & (ball > 0)
    print(f"   point Pgen : median {np.median(pt[ok]):.3g}", flush=True)
    print(f"   ball  Pgen : median {np.median(ball[ok]):.3g}  "
          f"(ball/point median {np.median(ball[ok]/np.maximum(pt[ok],1e-300)):.0f}x)", flush=True)
    print(f"   OBSERVED fuzzy incidence : median {np.median(obs[ok]):8.1f}", flush=True)
    print(f"   EXPECTED from ball Pgen  : median {np.median(exp_inc[ok]):8.1f}", flush=True)
    ratio = obs[ok] / np.maximum(exp_inc[ok], 1e-9)
    print(f"   obs/exp ratio: median {np.median(ratio):.2f}  p10 {np.percentile(ratio,10):.2f}"
          f"  p90 {np.percentile(ratio,90):.2f}", flush=True)
    lo = np.log10(np.maximum(obs[ok], .5))
    le = np.log10(np.maximum(exp_inc[ok], .5))
    if lo.std() > 0 and le.std() > 0:
        print(f"   r(log10 obs, log10 exp) = {np.corrcoef(lo, le)[0,1]:.3f}", flush=True)
    print("\n   VERDICT: obs/exp ~1 => the ball is fair and the paper's Pgen claim holds;\n"
          "            obs >> exp  => we over-match.", flush=True)


if __name__ == "__main__":
    main()
