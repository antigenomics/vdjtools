"""Emerson 2017 biomarker benchmark: CMV- and HLA-associated public TCRβs at cohort scale.

Reproduces the core analysis of Emerson et al. (*Nat Genet* 2017, doi:10.1038/ng.3822) with
:func:`vdjtools.biomarker.fisher_association`: an incidence-based Fisher's-exact screen over
the 786-subject ``isalgo/airr_hip`` cohort (the Emerson HIP repertoires, VDJtools format,
with per-subject CMV serostatus and 2-digit HLA-A/B typing). CMV-associated TCRβs are found
by a one-tailed enrichment test among CMV+ subjects; HLA-A*02-associated TCRβs by a two-tailed
test. Hits are then validated against the local VDJdb dump by matching CMV epitope + HLA allele.

Scale is the point: the cohort is streamed into a hive-partitioned Parquet dataset one sample
at a time (``ingest_cohort``), analysed as a single out-of-core ``polars`` LazyFrame
(``scan_cohort``), and the millions of per-feature Fisher tests are vectorised through the
hypergeometric tail — no per-feature Python loop, the cohort never fully in RAM.

Run::

    python examples/emerson_cmv_hla.py                     # full 786-subject cohort (~4 GB)
    python examples/emerson_cmv_hla.py --max-samples 150   # fast subset (progressive scaling)
    python examples/emerson_cmv_hla.py --with-1mm          # also run the metaclonotype (1mm) screen

Needs the ``[overlap]`` extra for ``--with-1mm``/1mm-vdjdb matching, ``huggingface_hub`` for the
data, and ``matplotlib`` for the volcano plots (``pip install 'vdjtools[overlap,examples]'``).
"""
from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S
from vdjtools.biomarker import fisher_association

REPO = "isalgo/airr_hip"
VDJDB = Path("/Users/mikesh/vcs/code/vdjdb-db/database/vdjdb.slim.txt")


def _peak_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1e6 if sys.platform == "darwin" else rss / 1e3


def _t(t0: float) -> str:
    return f"{time.perf_counter() - t0:6.1f}s  RSS {_peak_mb():6.0f} MB"


# ── data ────────────────────────────────────────────────────────────────────────────

def fetch(max_samples: int) -> tuple[Path, pl.DataFrame]:
    """Download metadata + the chosen ``corr/*.txt.gz`` repertoires; return (corr_dir, metadata)."""
    from huggingface_hub import hf_hub_download, snapshot_download

    meta_path = hf_hub_download(REPO, "metadata.txt", repo_type="dataset")
    # NB split on TAB only — the `race` field itself contains commas.
    meta = pl.read_csv(meta_path, separator="\t", infer_schema_length=0)
    known = meta.filter(pl.col("cmv").is_in(["+", "-"]))          # CMV-serotyped subjects
    if max_samples:                                              # balanced subset: both classes
        known = pl.concat([known.filter(pl.col("cmv") == c).head(max_samples // 2)
                           for c in ("+", "-")])
    ids = known["sample_id"].to_list()
    patterns = ["metadata.txt"] + [f"corr/{s}.txt.gz" for s in ids]
    root = Path(snapshot_download(REPO, repo_type="dataset", allow_patterns=patterns))
    return root / "corr", known


def phenotypes(meta: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Build the CMV (+/−) and HLA-A*02 (present/absent) per-subject phenotype tables."""
    cmv = meta.select(
        S_ID := "sample_id",
        pl.when(pl.col("cmv") == "+").then(True)
          .when(pl.col("cmv") == "-").then(False).otherwise(None).alias("cmv_pos"))
    a02 = meta.select(
        "sample_id",
        pl.when(pl.col("hla").is_null() | pl.col("hla").is_in(["", "NA"])).then(None)
          .otherwise(pl.col("hla").str.contains(r"HLA-A\*02")).alias("a02"))
    return cmv, a02


# ── vdjdb validation ────────────────────────────────────────────────────────────────

def load_vdjdb_cmv() -> pl.DataFrame:
    """Human TRB CMV-specific records from the local VDJdb slim dump (cdr3, V/J, epitope, MHC)."""
    v = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    return (v.filter((pl.col("gene") == "TRB") & (pl.col("species") == "HomoSapiens")
                     & pl.col("antigen.species").str.contains("CMV"))
            .select(pl.col("cdr3"),
                    S.strip_allele(pl.col("v.segm")).alias("vdjdb_v"),
                    pl.col("antigen.epitope").alias("epitope"),
                    pl.col("mhc.a").alias("mhc"))
            .unique())


def validate_cmv(hits: pl.DataFrame, vdjdb: pl.DataFrame, fdr: float) -> None:
    """Report overlap of significant CMV-enriched hits with VDJdb CMV entries (exact + 1mm)."""
    sig = hits.filter((pl.col("q_value") < fdr) & (pl.col("direction") == "enriched"))
    print(f"\n=== vdjdb validation — {sig.height} significant CMV-enriched TCRβ (q<{fdr}) ===")
    if sig.height == 0:
        print("  (no significant hits at this scale)")
        return

    # Exact CDR3 match (+ V-family agreement where both call a V).
    ex = sig.join(vdjdb, left_on=S.CDR3_AA, right_on="cdr3", how="inner")
    v_agree = ex.filter(pl.col(S.V_CALL) == pl.col("vdjdb_v")).height if "vdjdb_v" in ex.columns else 0
    a02 = ex.filter(pl.col("mhc").str.starts_with("HLA-A*02")).height
    print(f"  exact CDR3 in vdjdb-CMV : {ex['cdr3_aa'].n_unique():4d} / {sig.height} "
          f"hits   (V-family agrees on {v_agree}; HLA-A*02-restricted in vdjdb: {a02})")
    epi = (ex.group_by("epitope").len().sort("len", descending=True).head(5))
    for r in epi.iter_rows(named=True):
        print(f"      epitope {r['epitope']:<14s} {r['len']:3d} matches")
    for r in ex.sort("q_value").head(6).iter_rows(named=True):
        print(f"      {r['cdr3_aa']:<20s} {r.get('v_call',''):<10s} "
              f"OR={r['odds_ratio']:6.1f} q={r['q_value']:.1e}  ↔ {r['epitope']}/{r['mhc']}")

    # 1-mismatch match to vdjdb (bonus; needs the [overlap] extra).
    try:
        import vdjmatch.cluster as vc
        ours = sig[S.CDR3_AA].unique().to_list()
        ref = vdjdb["cdr3"].unique().to_list()
        pairs = vc.overlap(ours, ref, scope="1,0,0,1")
        n1 = pairs["a_idx"].n_unique() if pairs.height else 0
        print(f"  ≤1-mismatch to vdjdb-CMV: {n1:4d} / {sig.height} hits have a Hamming-1 CMV neighbour")
    except ImportError:
        print("  (skip 1mm-vdjdb match — install the [overlap] extra)")


# ── plotting ────────────────────────────────────────────────────────────────────────

def volcano(hits: pl.DataFrame, vdjdb_cdr3: set[str] | None, title: str, out: Path,
            fdr: float) -> None:
    """Volcano: log2 odds-ratio (x) vs −log10 p (y); significant + vdjdb-validated highlighted."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    h = hits.with_columns(
        (-pl.col("p_value").log10()).alias("y"),
        (pl.col("q_value") < fdr).alias("sig"))
    x = h["log2_or"].to_numpy()
    y = np.clip(h["y"].to_numpy(), 0, 320)
    sig = h["sig"].to_numpy()

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(x[~sig], y[~sig], s=4, c="#c8ccd4", alpha=0.4, linewidths=0, label="ns")
    ax.scatter(x[sig], y[sig], s=8, c="#d1495b", alpha=0.7, linewidths=0, label=f"q<{fdr}")
    if vdjdb_cdr3:
        # Highlight only *significant* hits that are also in vdjdb-CMV — the meaningful
        # validation. (Circling every vdjdb member floods the plot: at low incidence many
        # features share one (OR, p) coordinate, so a vdjdb hit lands on nearly all of them.)
        in_db = h[S.CDR3_AA].is_in(list(vdjdb_cdr3)).to_numpy()
        val = sig & in_db
        ax.scatter(x[val], y[val], s=30, facecolors="none", edgecolors="#00798c",
                   linewidths=1.3, label="sig & in vdjdb-CMV")
    ax.axhline(-np.log10(0.05), ls="--", c="k", lw=0.6, alpha=0.5)
    ax.axvline(0, ls="-", c="k", lw=0.6, alpha=0.3)
    ax.set_xlabel("log2 odds ratio  (enriched →)")
    ax.set_ylabel("−log10 p-value")
    ax.set_title(title)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"  volcano → {out}")


# ── driver ──────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-samples", type=int, default=0, help="0 = full 786-subject cohort")
    ap.add_argument("--min-incidence", type=int, default=2)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--with-1mm", action="store_true", help="also run the metaclonotype (1mm) CMV screen")
    ap.add_argument("--workdir", type=Path, default=Path("examples/.data/emerson"))
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    t0 = time.perf_counter()
    corr_dir, meta = fetch(args.max_samples)
    print(f"[0] fetch {meta.height} subjects ({(meta['cmv']=='+').sum()} CMV+/"
          f"{(meta['cmv']=='-').sum()} CMV−)   {_t(t0)}")
    cmv_ph, a02_ph = phenotypes(meta)

    cohort = args.workdir / "cohort"
    t0 = time.perf_counter()
    if any(cohort.glob("sample_id=*/*.parquet")):
        print(f"[1] ingest_cohort — reusing existing parquet cohort at {cohort}")
    else:
        vio.ingest_cohort(meta.select("sample_id", "cmv", "hla", "age", "sex"), corr_dir, cohort,
                          sample_col="sample_id", file_template="{sample}.txt.gz", fmt="vdjtools")
        print(f"[1] ingest_cohort → hive parquet (one sample in RAM)   {_t(t0)}")
    lf = vio.scan_cohort(cohort, join_metadata=False)

    # CMV: one-tailed enrichment among CMV+ (Emerson's setting). Full V+CDR3+J key.
    t0 = time.perf_counter()
    cmv = fisher_association(lf, cmv_ph, pheno_col="cmv_pos", alternative="greater",
                             min_incidence=args.min_incidence)
    n_sig = cmv.filter(pl.col("q_value") < args.fdr).height
    print(f"[2] CMV Fisher: {cmv.height} public TCRβ tested, {n_sig} sig (q<{args.fdr})   {_t(t0)}")
    for r in cmv.head(6).iter_rows(named=True):
        print(f"      {r['cdr3_aa']:<20s} {r['v_call']:<10s} {r['j_call']:<9s} "
              f"{r['n_pos_present']:3d}+/{r['n_neg_present']:3d}−  OR={r['odds_ratio']:6.1f}  "
              f"p={r['p_value']:.1e} q={r['q_value']:.1e}")

    # HLA-A*02: two-tailed (positive + negative association).
    t0 = time.perf_counter()
    hla = fisher_association(lf, a02_ph, pheno_col="a02", alternative="two-sided",
                             min_incidence=args.min_incidence)
    n_sig_h = hla.filter(pl.col("q_value") < args.fdr).height
    print(f"[3] HLA-A*02 Fisher: {hla.height} tested, {n_sig_h} sig (q<{args.fdr})   {_t(t0)}")

    vdjdb = load_vdjdb_cmv()
    validate_cmv(cmv, vdjdb, args.fdr)

    if args.with_1mm:
        t0 = time.perf_counter()
        # 1mm over public keys only (incidence≥min) — the 89M-unique full set is intractable to
        # cluster whole; public-key restriction is the documented benchmark approximation.
        pub = (lf.group_by([S.CDR3_AA, S.V_CALL, S.J_CALL]).agg(pl.len().alias("_n"))
               .filter(pl.col("_n") >= args.min_incidence).select([S.CDR3_AA, S.V_CALL, S.J_CALL]))
        cmv1 = fisher_association(lf.join(pub, on=[S.CDR3_AA, S.V_CALL, S.J_CALL], how="semi"),
                                  cmv_ph, pheno_col="cmv_pos", alternative="greater",
                                  min_incidence=args.min_incidence, match="1mm")
        n1 = cmv1.filter(pl.col("q_value") < args.fdr).height
        print(f"[4] CMV 1mm metaclonotype screen: {cmv1.height} metas, {n1} sig   {_t(t0)}")

    if not args.no_plots:
        vdjdb_cdr3 = set(vdjdb["cdr3"].to_list())
        volcano(cmv, vdjdb_cdr3, "CMV-associated TCRβ (Emerson HIP)",
                args.workdir / "emerson_cmv_volcano.png", args.fdr)
        volcano(hla, None, "HLA-A*02-associated TCRβ (Emerson HIP)",
                args.workdir / "emerson_hla_volcano.png", args.fdr)

    # runnable check: p-values are valid probabilities and the top CMV hit is enriched.
    assert cmv["p_value"].is_between(0, 1).all()
    if cmv.height:
        assert cmv.row(0, named=True)["direction"] == "enriched"
    print(f"\nOK — Emerson CMV/HLA screen complete   total peak RSS {_peak_mb():.0f} MB")


if __name__ == "__main__":
    main()
