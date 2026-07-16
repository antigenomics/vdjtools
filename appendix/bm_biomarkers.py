"""COVID biomarker discovery on the full FMBA cohort — OUR list, two units, two scopes.

WHAT THIS IS. A list of individual biomarker CLONOTYPES and how many pass FDR. The 1-mismatch
is only a trick to estimate each candidate's incidence more precisely -- the candidate keeps its
identity throughout. Nothing is merged into anything. (Hamming graphs / metaclonotypes /
classifiers built ON this list are a separate downstream step and are out of scope here.)

SCOPES
  aa    : candidate = CDR3aa;      a donor matches if it carries any CDR3aa within Hamming-1.
  aa_v  : candidate = (CDR3aa, V); a donor matches if it carries any CDR3aa within Hamming-1
          AND the same V gene. More precise -- V pins the germline half of the contact surface,
          so the 1mm ball stops absorbing unrelated rearrangements that merely look similar.

UNITS -- the same candidate scored two ways:

  donor    : Fisher exact on presence/absence across donors (Emerson 2017 / the paper).
             Counts are tiny (<=1,240) so exactness is cheap and correct here.

  rows     : the REARRANGEMENT unit. k_pos = total matching rows summed over COVID+ donors,
             out of N_pos = all rows in the COVID+ part of the dataset; same for COVID-.
             A row = a unique nt rearrangement in a repertoire = one recombination event, so
             one CDR3aa reached by 3 nt variants is 3 events (convergent recombination). This
             is the unit Pgen is defined on -- and it is NOT reads.

             N is ~1.7e7, so an exact hypergeometric (factorials) is the wrong tool. Two smooth
             tests instead:
               * conditional binomial (the standard two-sample Poisson test, Przyborowski &
                 Wilenski 1940): given k = k_pos + k_neg, k_pos ~ Binom(k, N_pos/(N_pos+N_neg)).
                 Smooth (regularized incomplete beta), valid at small AND large k.
               * G-test (likelihood ratio, 2*sum O*ln(O/E) ~ chi2_1) as a cross-check.

             CAVEAT, reported not hidden: rows within a donor are NOT independent (clonal
             expansion + donor effects), so the row-based p is anticonservative relative to the
             donor test. `overdispersion` below quantifies it per candidate (what fraction of
             k_pos comes from its single biggest donor). Both units are reported side by side.

2026-07-16.
"""
from __future__ import annotations

import argparse
import resource
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import binom, chi2

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19")
MIN_READS = 10_000


def _rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2


def t(t0, m):
    print(f"  [{m}] {time.perf_counter()-t0:7.1f}s  peakRSS {_rss():.1f} GB", flush=True)


def binom_p(k_pos, k_neg, n_pos_rows, n_neg_rows):
    """Conditional binomial (two-sample Poisson) test, one-tailed for enrichment in +.

    Given the total k = k_pos + k_neg matching rows, under H0 each falls in the + arm with
    probability p0 = N_pos/(N_pos+N_neg). Smooth: binom.sf is a regularized incomplete beta.
    """
    k = k_pos + k_neg
    p0 = n_pos_rows / (n_pos_rows + n_neg_rows)
    return binom.sf(k_pos - 1, k, p0)


def gtest_p(a, b, c, d):
    """G-test (likelihood-ratio) on the 2x2 -- smooth, no factorials. Two-sided chi2_1."""
    a, b, c, d = (np.asarray(x, float) for x in (a, b, c, d))
    n = a + b + c + d
    obs = np.stack([a, b, c, d])
    exp = np.stack([(a + b) * (a + c), (a + b) * (b + d), (c + d) * (a + c), (c + d) * (b + d)]) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(obs > 0, obs * np.log(np.where(exp > 0, obs / exp, 1.0)), 0.0)
    return chi2.sf(2.0 * term.sum(axis=0), 1)


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
        print(f"  {miss} metadata rows without a file", flush=True)
    return pl.concat(frames, how="vertical_relaxed")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRB", choices=["TRA", "TRB"])
    ap.add_argument("--scope", default="aa", choices=["aa", "aa_v"])
    ap.add_argument("--min-samples", type=int, default=2)
    ap.add_argument("--out", default="/projects/biomarkers/results")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    tag = f"{args.chain}_{args.scope}"
    print(f"=== BIOMARKERS covid19 {args.chain}  scope={args.scope}  (QC reads>={MIN_READS:,}) ===",
          flush=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    meta = meta.with_columns(pl.col("reads").cast(pl.Int64, strict=False))
    n_all = meta.filter(pl.col("locus") == args.chain).height
    meta = meta.filter(pl.col("reads") >= MIN_READS)
    print(f"  QC: {args.chain} {n_all} -> {meta.filter(pl.col('locus')==args.chain).height} "
          f"samples pass reads>={MIN_READS:,}", flush=True)

    # healthy == precovid: both are controls.
    design = (meta.filter(pl.col("locus") == args.chain)
              .select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                      pl.when(pl.col("COVID_status") == "COVID").then(True)
                        .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                        .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))

    t0 = time.perf_counter()
    cohort = load(meta, args.chain)
    t(t0, "ingest")

    # ---- rows: one row per (CDR3aa, V, donor) carrying its rearrangement count -------------
    t0 = time.perf_counter()
    vgene = pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene")
    rows = (cohort.lazy()
            .filter(pl.col(S.JUNCTION_AA).is_not_null()
                    & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .with_columns(vgene)
            .group_by([S.JUNCTION_AA, "v_gene", "sample_id"])
            .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
            .join(design.lazy(), on="sample_id", how="inner")
            .collect(engine="streaming"))
    t(t0, "rows table")

    seen = design.join(rows.select("sample_id").unique(), on="sample_id", how="semi")
    n_pos, n_neg = int(seen["_pos"].sum()), seen.height - int(seen["_pos"].sum())
    # total rows (rearrangements) in each arm -- the row-test denominators
    arm = rows.group_by("_pos").agg(pl.col("n_rearr").sum().alias("N"))
    N_pos = int(arm.filter(pl.col("_pos"))["N"][0])
    N_neg = int(arm.filter(~pl.col("_pos"))["N"][0])
    print(f"  ANALYSED donors: {n_pos} COVID / {n_neg} control (healthy+precovid)", flush=True)
    print(f"  TOTAL ROWS (rearrangements): COVID+ {N_pos:,} / COVID- {N_neg:,} "
          f"= {N_pos+N_neg:,}", flush=True)

    # ---- candidates: the paper's filter, >=2 nt variants OR >=2 samples -------------------
    t0 = time.perf_counter()
    keycols = [S.JUNCTION_AA] + (["v_gene"] if args.scope == "aa_v" else [])
    cand = (rows.lazy().group_by(keycols)
            .agg(pl.col("n_rearr").sum().alias("n_nt_variants"),
                 pl.col("sample_id").n_unique().alias("exact_incidence"))
            .filter((pl.col("n_nt_variants") >= 2)
                    | (pl.col("exact_incidence") >= args.min_samples))
            .collect(engine="streaming")
            .sort(keycols))
    t(t0, "candidates")
    print(f"  candidates ({'+'.join(keycols)}): {cand.height:,}", flush=True)

    # ---- the 1mm search: candidate CDR3aa -> every cohort CDR3aa within Hamming-1 ---------
    import vdjmatch.cluster as vc
    universe = rows[S.JUNCTION_AA].unique().sort().to_list()
    cq = cand[S.JUNCTION_AA].unique().sort().to_list()
    t0 = time.perf_counter()
    pairs = vc.overlap(cq, universe, scope="1,0,0,1", threads=0)
    t(t0, "1mm search")
    print(f"  neighbour pairs: {pairs.height:,}  ({pairs.height/max(len(cq),1):.1f}/CDR3aa)",
          flush=True)

    qmap = pl.DataFrame({"a_idx": np.arange(len(cq), dtype=np.int64), "cand_aa": cq})
    umap = pl.DataFrame({"b_idx": np.arange(len(universe), dtype=np.int64),
                         S.JUNCTION_AA: universe})
    nb = (pairs.join(qmap, on="a_idx", how="inner").join(umap, on="b_idx", how="inner")
          .select("cand_aa", S.JUNCTION_AA))                      # candidate -> neighbour aa

    # ---- match donors + rows to each candidate -------------------------------------------
    t0 = time.perf_counter()
    if args.scope == "aa":
        j = nb.join(rows, on=S.JUNCTION_AA, how="inner").rename({"cand_aa": "cand"})
        gcols = ["cand"]
    else:
        # V must match EXACTLY: attach the candidate's V, then require rows to carry it.
        j = (nb.join(cand.select(S.JUNCTION_AA, "v_gene")
                     .rename({S.JUNCTION_AA: "cand_aa", "v_gene": "cand_v"}),
                     on="cand_aa", how="inner")
             .join(rows, left_on=[S.JUNCTION_AA, "cand_v"],
                   right_on=[S.JUNCTION_AA, "v_gene"], how="inner")
             .rename({"cand_aa": "cand"}))
        gcols = ["cand", "cand_v"]

    agg = (j.group_by(gcols).agg(
        pl.col("sample_id").n_unique().alias("incidence"),
        pl.col("sample_id").filter(pl.col("_pos")).n_unique().alias("a"),
        pl.col("n_rearr").filter(pl.col("_pos")).sum().alias("k_pos"),
        pl.col("n_rearr").filter(~pl.col("_pos")).sum().alias("k_neg"),
        pl.col("n_rearr").max().alias("max_donor_rows")))
    t(t0, "match donors+rows")

    a = agg["a"].to_numpy().astype(np.int64)
    b = (agg["incidence"].to_numpy().astype(np.int64) - a)
    c, d = n_pos - a, n_neg - b
    k_pos = agg["k_pos"].fill_null(0).to_numpy().astype(np.int64)
    k_neg = agg["k_neg"].fill_null(0).to_numpy().astype(np.int64)

    p_donor = stats.fisher_p(a, b, c, d, alternative="greater")
    p_rows = binom_p(k_pos, k_neg, N_pos, N_neg)
    p_g = gtest_p(k_pos, N_pos - k_pos, k_neg, N_neg - k_neg)

    res = agg.with_columns(
        pl.Series("n_pos_present", a), pl.Series("n_neg_present", b),
        pl.Series("odds_ratio", stats.odds_ratio(a, b, c, d)),
        pl.Series("log2_or", np.log2(stats.odds_ratio(a, b, c, d))),
        # rate ratio on the ROW unit: (k_pos/N_pos) / (k_neg/N_neg)
        pl.Series("rate_ratio", ((k_pos / N_pos) / np.maximum(k_neg / N_neg, 1e-30))),
        pl.Series("p_donor", p_donor), pl.Series("q_donor", stats.fdr_bh(p_donor)),
        pl.Series("p_rows", p_rows), pl.Series("q_rows", stats.fdr_bh(p_rows)),
        pl.Series("p_gtest", p_g),
        # what fraction of k_pos comes from the single biggest donor -> overdispersion flag
        pl.Series("overdispersion",
                  agg["max_donor_rows"].to_numpy() / np.maximum(k_pos, 1)),
    ).sort("p_donor")

    print(f"\nfeatures tested : {res.height:,}", flush=True)
    print("\n--- UNIT = DONOR (Fisher exact on incidence) ---", flush=True)
    for lab, col, thr in (("q<0.01", "q_donor", 0.01), ("q<0.05", "q_donor", 0.05),
                          ("P<1e-4", "p_donor", 1e-4)):
        print(f"  {lab:8s}: {res.filter(pl.col(col) < thr).height:,}", flush=True)
    print("\n--- UNIT = ROWS (conditional binomial on rearrangements) ---", flush=True)
    for lab, col, thr in (("q<0.01", "q_rows", 0.01), ("q<0.05", "q_rows", 0.05),
                          ("P<1e-4", "p_rows", 1e-4)):
        print(f"  {lab:8s}: {res.filter(pl.col(col) < thr).height:,}", flush=True)
    both = res.filter((pl.col("q_donor") < 0.01) & (pl.col("q_rows") < 0.01)).height
    print(f"\n  BOTH units q<0.01: {both:,}", flush=True)
    od = res.filter(pl.col("q_rows") < 0.01)["overdispersion"]
    if od.len():
        print(f"  overdispersion of row-significant hits (frac of k_pos from top donor): "
              f"median {od.median():.3f}  p90 {od.quantile(0.9):.3f}", flush=True)

    show = gcols + ["incidence", "n_pos_present", "n_neg_present", "odds_ratio", "k_pos",
                    "k_neg", "rate_ratio", "p_donor", "q_donor", "p_rows", "q_rows"]
    print("\ntop 15 by donor test:", flush=True)
    print(res.head(15).select(show), flush=True)
    print("\ntop 15 by row test:", flush=True)
    print(res.sort("p_rows").head(15).select(show), flush=True)

    res.write_parquet(out / f"bm_{tag}.parquet")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    for ax, (xc, pc, qc, name) in zip(axes, (
            ("log2_or", "p_donor", "q_donor", "unit = donor (Fisher, incidence)"),
            ("rate_ratio", "p_rows", "q_rows", "unit = rows (binomial, rearrangements)"))):
        x = res[xc].to_numpy()
        if xc == "rate_ratio":
            x = np.log2(np.clip(x, 1e-6, 1e6))
        y = np.clip(-np.log10(np.clip(res[pc].to_numpy(), 1e-320, 1)), 0, 320)
        sig = res[qc].to_numpy() < 0.01
        ax.scatter(x[~sig], y[~sig], s=4, c="#c8ccd4", alpha=.35, lw=0, label="ns")
        ax.scatter(x[sig], y[sig], s=12, c="#D55E00", alpha=.8, lw=0,
                   label=f"q<0.01 (n={int(sig.sum()):,})")
        ax.axhline(4, ls=":", c="#0072B2", lw=.8, label="P<1e-4")
        ax.axvline(0, ls="-", c="k", lw=.6, alpha=.3)
        ax.set(xlabel="log2 ratio  (enriched in COVID →)", ylabel="−log10 p", title=name)
        ax.legend(fontsize=8, frameon=False, loc="upper left")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"COVID biomarkers {args.chain}  scope={args.scope}  "
                 f"({n_pos}+/{n_neg}− donors; {N_pos:,}/{N_neg:,} rows)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / f"bm_volcano_{tag}.png", dpi=140)
    print(f"\n-> {out}/bm_{tag}.parquet + volcano", flush=True)


if __name__ == "__main__":
    main()
