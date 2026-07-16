"""P1 — in-silico alpha-beta pairing, and the ORIGIN of the pairs.

Howie 2015 (pairSEQ) pairs chains by co-occurrence across wells; here donors are the wells.
For an alpha biomarker i and a beta biomarker j:

    theta = n * n_ij / (n_i * n_j)          lift over independence
    Fisher on the 2x2 (carries alpha? x carries beta?) across donors, one-tailed.

THE ORIGIN QUESTION (owner / De Witt 2018). A co-occurring alpha-beta pair on its own is weak
evidence -- two clonotypes that are each COVID-associated will co-occur simply because both
track COVID status. The pair becomes interesting when its co-occurrence is ALSO HLA-restricted:
that is the imprint of a specific presented epitope rather than of shared exposure alone. So
each significant pair is asked three questions:

  1. CONFOUND: does the pair still co-occur AMONG COVID+ DONORS ONLY? If theta collapses to ~1
     inside the COVID arm, the pair was only ever co-tracking infection status, not each other.
  2. HLA: are the co-carriers enriched for one HLA allele (Fisher per allele)? An HLA-restricted
     co-occurring pair is the De Witt signature.
  3. VDJdb: do the alpha and the beta map to the SAME epitope in VDJdb's paired records? That is
     direct confirmation the pair is a real receptor, not a statistical coincidence.

Question 1 is the important one and it is a genuine confound check, not decoration: without it a
"pair" list is just the biomarker list crossed with itself. 2026-07-16.
"""
from __future__ import annotations

import argparse
import time
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


def load(meta, chain):
    frames = []
    for r in meta.filter(pl.col("locus") == chain).iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{chain}.*"))
        if hits:
            frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                          .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    return pl.concat(frames, how="vertical_relaxed")


def membership(rows, cands, donors):
    """boolean donor x candidate matrix under V + CDR3aa±1mm."""
    import vdjmatch.cluster as vc
    universe = rows[S.JUNCTION_AA].unique().sort().to_list()
    cq = cands["cand"].to_list()
    pairs = vc.overlap(cq, universe, scope="1,0,0,1", threads=0)
    qmap = pl.DataFrame({"a_idx": np.arange(len(cq), dtype=np.int64), "cand": cq,
                         "cand_v": cands["v_gene"].to_list()})
    umap = pl.DataFrame({"b_idx": np.arange(len(universe), dtype=np.int64),
                         S.JUNCTION_AA: universe})
    j = (pairs.join(qmap, on="a_idx").join(umap, on="b_idx")
         .join(rows, left_on=[S.JUNCTION_AA, "cand_v"], right_on=[S.JUNCTION_AA, "v_gene"],
               how="inner")
         .select("cand", "cand_v", "sample_id").unique())
    key = (j["cand"] + "|" + j["cand_v"]).to_list()
    feats = sorted(set(key))
    fidx = {f: i for i, f in enumerate(feats)}
    didx = {d: i for i, d in enumerate(donors)}
    M = np.zeros((len(donors), len(feats)), dtype=bool)
    for f, s in zip(key, j["sample_id"].to_list()):
        if s in didx:
            M[didx[s], fidx[f]] = True
    return M, feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=400, help="top biomarkers per chain by p_donor")
    ap.add_argument("--min-co", type=int, default=5, help="min co-carrying donors")
    args = ap.parse_args()
    print(f"=== P1 alpha-beta pairing: top-{args.top} biomarkers/chain, V+CDR3aa±1mm ===",
          flush=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    meta = meta.with_columns(pl.col("reads").cast(pl.Int64, strict=False)).filter(
        pl.col("reads") >= MIN_READS)
    design = (meta.select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                          pl.when(pl.col("COVID_status") == "COVID").then(True)
                            .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                            .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))

    mats, feats_by = {}, {}
    for chain in ("TRA", "TRB"):
        bm = pl.read_parquet(RES / f"bm_{chain}_aa_v.parquet")
        cands = bm.sort("p_donor").head(args.top).select(
            pl.col("cand"), pl.col("cand_v").alias("v_gene"))
        t0 = time.perf_counter()
        cohort = load(meta, chain)
        rows = (cohort.lazy()
                .filter(pl.col(S.JUNCTION_AA).is_not_null()
                        & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
                .with_columns(pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene"))
                .select(S.JUNCTION_AA, "v_gene", "sample_id").unique()
                .join(design.lazy(), on="sample_id", how="inner")
                .collect(engine="streaming"))
        donors = sorted(set(rows["sample_id"].to_list()))
        M, feats = membership(rows, cands, donors)
        mats[chain] = (M, donors)
        feats_by[chain] = feats
        print(f"  {chain}: {M.shape[0]} donors x {M.shape[1]} biomarkers "
              f"[{time.perf_counter()-t0:.0f}s]", flush=True)

    # donors sequenced for BOTH chains
    common = sorted(set(mats["TRA"][1]) & set(mats["TRB"][1]))
    print(f"\n  donors with BOTH chains: {len(common):,}", flush=True)
    ia = [mats["TRA"][1].index(d) for d in common]
    ib = [mats["TRB"][1].index(d) for d in common]
    A = mats["TRA"][0][ia, :]
    B = mats["TRB"][0][ib, :]
    pos = np.array([bool(design.filter(pl.col("sample_id") == d)["_pos"][0]) for d in common])
    n = len(common)
    print(f"  of which COVID+: {int(pos.sum()):,}", flush=True)

    n_ab = (A.astype(np.int32).T @ B.astype(np.int32))          # co-carriage
    n_a = A.sum(0)[:, None]
    n_b = B.sum(0)[None, :]
    theta = np.where((n_a * n_b) > 0, n * n_ab / np.maximum(n_a * n_b, 1), 0.0)

    ii, jj = np.where(n_ab >= args.min_co)
    print(f"  pairs with >={args.min_co} co-carriers: {len(ii):,} "
          f"(of {A.shape[1]*B.shape[1]:,})", flush=True)
    if not len(ii):
        return
    a = n_ab[ii, jj].astype(np.int64)
    b = (n_a[:, 0][ii] - a).astype(np.int64)
    c = (n_b[0, :][jj] - a).astype(np.int64)
    d = (n - a - b - c).astype(np.int64)
    p = stats.fisher_p(a, b, c, d, alternative="greater")
    q = stats.fdr_bh(p)

    # CONFOUND CHECK: same pair, COVID+ donors only
    Ap, Bp = A[pos, :], B[pos, :]
    npos = int(pos.sum())
    nab_p = (Ap.astype(np.int32).T @ Bp.astype(np.int32))
    na_p = Ap.sum(0)[:, None]
    nb_p = Bp.sum(0)[None, :]
    theta_pos = np.where((na_p * nb_p) > 0, npos * nab_p / np.maximum(na_p * nb_p, 1), 0.0)
    ap_ = nab_p[ii, jj].astype(np.int64)
    bp_ = (na_p[:, 0][ii] - ap_).astype(np.int64)
    cp_ = (nb_p[0, :][jj] - ap_).astype(np.int64)
    dp_ = (npos - ap_ - bp_ - cp_).astype(np.int64)
    p_pos = stats.fisher_p(ap_, bp_, cp_, dp_, alternative="greater")

    res = pl.DataFrame({
        "alpha": [feats_by["TRA"][i] for i in ii],
        "beta": [feats_by["TRB"][j] for j in jj],
        "n_a": n_a[:, 0][ii], "n_b": n_b[0, :][jj], "n_ab": a,
        "theta": theta[ii, jj], "p_value": p, "q_value": q,
        "theta_covid_only": theta_pos[ii, jj], "p_covid_only": p_pos,
    }).sort("p_value")

    sig = res.filter(pl.col("q_value") < 0.05)
    print(f"\n  SIGNIFICANT PAIRS q<0.05: {sig.height:,}", flush=True)
    surv = sig.filter(pl.col("p_covid_only") < 0.05)
    print(f"  ...of which SURVIVE the COVID-only confound check (p_covid_only<0.05): "
          f"{surv.height:,}", flush=True)
    print(f"     (a pair that fails this was only co-tracking infection status)", flush=True)
    if sig.height:
        print(f"  theta: median {sig['theta'].median():.2f}  max {sig['theta'].max():.2f}; "
              f"within-COVID theta median {sig['theta_covid_only'].median():.2f}", flush=True)
        print("\n  top 15 pairs:", flush=True)
        print(sig.head(15), flush=True)

    # VDJdb paired records: same epitope for both chains?
    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    dbp = (db.filter(pl.col("species") == "HomoSapiens")
           .with_columns(pl.col("complex.id").str.split(",").alias("_cid"))
           .explode("_cid").with_columns(pl.col("_cid").str.strip_chars())
           .filter(pl.col("_cid") != "0"))
    al = dbp.filter(pl.col("gene") == "TRA").select("_cid", pl.col("cdr3").alias("a_cdr3"),
                                                    pl.col("antigen.epitope").alias("ep"))
    be = dbp.filter(pl.col("gene") == "TRB").select("_cid", pl.col("cdr3").alias("b_cdr3"))
    vpairs = al.join(be, on="_cid").select("a_cdr3", "b_cdr3", "ep").unique()
    known = set(zip(vpairs["a_cdr3"].to_list(), vpairs["b_cdr3"].to_list()))
    print(f"\n  VDJdb paired records (exploded complex.id): {vpairs.height:,}", flush=True)
    hits = [(x, y) for x, y in zip(sig["alpha"].to_list(), sig["beta"].to_list())
            if (x.split("|")[0], y.split("|")[0]) in known]
    print(f"  our significant pairs found as REAL VDJdb receptors: {len(hits)}", flush=True)
    for h in hits[:10]:
        print(f"      {h}", flush=True)

    res.write_parquet(RES / "pair_covid.parquet")
    print(f"\n-> {RES}/pair_covid.parquet ({res.height:,} pairs)", flush=True)


if __name__ == "__main__":
    main()
