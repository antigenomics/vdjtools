"""H1/H2 — COVID association RESTRICTED to HLA carriers (per allele, top-N strata).

A TCR biomarker only works in donors that can present its epitope. Pooling all HLA types
dilutes every HLA-restricted response by the fraction of donors lacking the presenting allele;
restricting to carriers should CONCENTRATE the signal. That is the test.

For each allele: subset donors to its carriers, then run COVID+ vs control inside that stratum
(healthy+precovid = control, per owner). Reported per allele so the strata are comparable.

Both units are shipped as options (owner's call), and the row unit is reported under BOTH nulls:
  p0_rows  = N_pos/(N_pos+N_neg)  -- total sequencing effort
  p0_donor = n_pos/(n_pos+n_neg)  -- correct for clonotypes whose rows do not scale with depth
They differ by ~15-20% on this cohort because controls are ~1.4-1.5x deeper per donor, so any
clonotype carried ~once per donor regardless of depth is biased under the row null.

BIOLOGY NOTE (owner): this cohort is CONVALESCENT -- the responding clones may be largely
contracted by sampling. A weak effect is the prior expectation here, not evidence of a bug.

The expensive step (the 1mm search) is done ONCE over the whole cohort and reused for every
stratum; only the donor subset changes. 2026-07-16.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import binom, fisher_exact

from vdjtools import io as vio
from vdjtools.biomarker import stats
from vdjtools.io import schema as S

ROOT = Path("/projects/biomarkers/raw/covid19")
RES = Path("/projects/biomarkers/results")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")
MIN_READS = 10_000


def oracle(chain, pattern):
    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    return set(db.filter((pl.col("gene") == chain) & (pl.col("species") == "HomoSapiens")
                         & pl.col("antigen.species").str.contains(pattern))["cdr3"]
               .unique().to_list())


def enrich(sig, rest, truth):
    a = len(sig & truth); b = len(sig) - a
    c = len(rest & truth); d = len(rest) - c
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    return a, orr, p


def load(meta, chain):
    frames = []
    for r in meta.filter(pl.col("locus") == chain).iter_rows(named=True):
        hits = sorted(ROOT.glob(f"{r['file_id']}.{chain}.*"))
        if hits:
            frames.append(vio.read(str(hits[0]), fmt="vdjtools")
                          .with_columns(pl.lit(str(r["sample_id"])).alias("sample_id")))
    return pl.concat(frames, how="vertical_relaxed")


def hla_carriers(meta, chain, locus):
    """sample_id -> set of alleles at `locus` (both slots), 4-digit as typed."""
    cols = [c for c in meta.columns if c.startswith(f"HLA-{locus}_")]
    m = meta.filter(pl.col("locus") == chain)
    return (m.select(pl.col("sample_id").cast(pl.String), *[pl.col(c) for c in cols])
            .unpivot(index="sample_id", on=cols, variable_name="slot", value_name="allele")
            .filter(pl.col("allele").is_not_null() & (pl.col("allele") != "")
                    & (pl.col("allele") != "nan"))
            .unique(["sample_id", "allele"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRA")
    ap.add_argument("--locus", default="A", help="HLA locus: A B C DRB1 DQB1 DPB1")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--two-digit", action="store_true",
                    help="collapse to 2-digit groups (A*02 rather than A*02:01)")
    ap.add_argument("--complement", action="store_true",
                    help="run in NON-carriers -- the control that proves restriction. A hit that "
                         "is significant in carriers AND absent here is genuinely restricted; a "
                         "hit significant in both was only ever tracking cohort size/power.")
    args = ap.parse_args()
    print(f"=== COVID x HLA-{args.locus} {args.chain}  top-{args.top}"
          f"{' (2-digit)' if args.two_digit else ''}  key=V+CDR3aa±1mm ===", flush=True)

    meta = pl.read_csv(ROOT / "metadata.tsv", separator="\t", infer_schema_length=0)
    meta = meta.with_columns(pl.col("reads").cast(pl.Int64, strict=False)).filter(
        pl.col("reads") >= MIN_READS)
    design = (meta.filter(pl.col("locus") == args.chain)
              .select(pl.col("sample_id").cast(pl.String).alias("sample_id"),
                      pl.when(pl.col("COVID_status") == "COVID").then(True)
                        .when(pl.col("COVID_status").is_in(["healthy", "precovid"])).then(False)
                        .otherwise(None).alias("_pos"))
              .drop_nulls("_pos").unique(subset="sample_id"))

    hla = hla_carriers(meta, args.chain, args.locus).join(design, on="sample_id", how="inner")
    if args.two_digit:
        hla = hla.with_columns(pl.col("allele").str.replace(r"(\*\d+):.*", r"$1").alias("allele")
                               ).unique(["sample_id", "allele"])
    top = (hla.group_by("allele").agg(pl.len().alias("carriers"),
                                      pl.col("_pos").sum().alias("covid"))
           .with_columns((pl.col("carriers") - pl.col("covid")).alias("control"))
           .filter((pl.col("covid") >= 20) & (pl.col("control") >= 20))
           .sort("carriers", descending=True).head(args.top))
    print("  strata:", flush=True)
    print(top, flush=True)

    t0 = time.perf_counter()
    cohort = load(meta, args.chain)
    vgene = pl.col(S.V_CALL).str.split("*").list.first().alias("v_gene")
    rows = (cohort.lazy()
            .filter(pl.col(S.JUNCTION_AA).is_not_null()
                    & ~pl.col(S.JUNCTION_AA).str.contains(r"[*_]"))
            .with_columns(vgene)
            .group_by([S.JUNCTION_AA, "v_gene", "sample_id"])
            .agg(pl.col(S.JUNCTION_NT).n_unique().alias("n_rearr"))
            .join(design.lazy(), on="sample_id", how="inner")
            .collect(engine="streaming"))
    print(f"[ingest {time.perf_counter()-t0:.0f}s]", flush=True)

    # ONE 1mm search for every stratum
    import vdjmatch.cluster as vc
    universe = rows[S.JUNCTION_AA].unique().sort().to_list()
    t0 = time.perf_counter()
    pairs = vc.overlap(universe, universe, scope="1,0,0,1", threads=0)
    umap = pl.DataFrame({"idx": np.arange(len(universe), dtype=np.int64), "aa": universe})
    nb = (pairs.join(umap.rename({"idx": "a_idx", "aa": "cand"}), on="a_idx")
          .join(umap.rename({"idx": "b_idx", "aa": S.JUNCTION_AA}), on="b_idx")
          .select("cand", S.JUNCTION_AA))
    print(f"[1mm search {time.perf_counter()-t0:.0f}s] pairs {nb.height:,}", flush=True)

    cov = oracle(args.chain, "SARS-CoV-2")
    out_all = []
    all_typed = set(hla["sample_id"].to_list())
    for row in top.iter_rows(named=True):
        allele = row["allele"]
        carriers = set(hla.filter(pl.col("allele") == allele)["sample_id"].to_list())
        # complement is drawn from HLA-TYPED donors only, so the two arms are comparable
        who = (all_typed - carriers) if args.complement else carriers
        allele = f"NOT({allele})" if args.complement else allele
        r = rows.filter(pl.col("sample_id").is_in(who))
        dep = r.group_by(["sample_id", "_pos"]).agg(pl.col("n_rearr").sum().alias("N_i"))
        n_pos = int(dep["_pos"].sum()); n_neg = dep.height - n_pos
        N_pos = int(dep.filter(pl.col("_pos"))["N_i"].sum())
        N_neg = int(dep.filter(~pl.col("_pos"))["N_i"].sum())
        if not (n_pos and n_neg):
            continue
        p0_rows = N_pos / (N_pos + N_neg)
        p0_donor = n_pos / (n_pos + n_neg)

        cv = r.select(pl.col(S.JUNCTION_AA).alias("cand"), "v_gene").unique()
        j = (nb.join(cv, on="cand", how="inner")
             .join(r, on=[S.JUNCTION_AA, "v_gene"], how="inner"))
        agg = (j.group_by(["cand", "v_gene"]).agg(
            pl.col("sample_id").filter(pl.col("_pos")).n_unique().alias("a"),
            pl.col("sample_id").filter(~pl.col("_pos")).n_unique().alias("b"),
            pl.col("n_rearr").filter(pl.col("_pos")).sum().alias("k_pos"),
            pl.col("n_rearr").filter(~pl.col("_pos")).sum().alias("k_neg"))
            .filter(pl.col("a") + pl.col("b") >= 2))
        a = agg["a"].to_numpy().astype(np.int64)
        b = agg["b"].to_numpy().astype(np.int64)
        kp = agg["k_pos"].fill_null(0).to_numpy().astype(np.int64)
        kn = agg["k_neg"].fill_null(0).to_numpy().astype(np.int64)
        p_d = stats.fisher_p(a, b, n_pos - a, n_neg - b, alternative="greater")
        p_r = binom.sf(kp - 1, kp + kn, p0_rows)
        p_rd = binom.sf(kp - 1, kp + kn, p0_donor)
        res = agg.with_columns(
            pl.lit(allele).alias("allele"),
            pl.Series("odds_ratio", stats.odds_ratio(a, b, n_pos - a, n_neg - b)),
            pl.Series("p_donor", p_d), pl.Series("q_donor", stats.fdr_bh(p_d)),
            pl.Series("p_rows", p_r), pl.Series("q_rows", stats.fdr_bh(p_r)),
            pl.Series("p_rows_dn", p_rd), pl.Series("q_rows_dn", stats.fdr_bh(p_rd)))

        nd = res.filter(pl.col("q_donor") < 0.01).height
        nr = res.filter(pl.col("q_rows") < 0.01).height
        nrd = res.filter(pl.col("q_rows_dn") < 0.01).height
        allc = set(res["cand"].to_list())
        sig_d = set(res.filter(pl.col("q_donor") < 0.01)["cand"].to_list())
        cd = enrich(sig_d, allc - sig_d, cov) if sig_d else (0, 0, 1)
        print(f"\n  --- {allele}: {n_pos}+/{n_neg}- carriers; {res.height:,} tested; "
              f"depth {N_pos/n_pos:,.0f} vs {N_neg/n_neg:,.0f}", flush=True)
        print(f"      donor q<0.01 {nd:6,} | rows(rowNull) {nr:6,} | rows(donorNull) {nrd:6,}",
              flush=True)
        print(f"      VDJdb SARS-CoV-2 among donor-sig: {cd[0]} (OR={cd[1]:.1f})", flush=True)
        if nd:
            print(res.filter(pl.col("q_donor") < 0.01).sort("p_donor").head(6)
                  .select("cand", "v_gene", "a", "b", "odds_ratio", "p_donor", "q_donor"),
                  flush=True)
        out_all.append(res.filter((pl.col("q_donor") < 0.05) | (pl.col("q_rows_dn") < 0.05)))

    if out_all:
        allres = pl.concat(out_all, how="vertical_relaxed")
        f = RES / f"hla_{args.chain}_{args.locus}{'_2d' if args.two_digit else ''}.parquet"
        allres.write_parquet(f)
        print(f"\n-> {f}  ({allres.height:,} rows)", flush=True)


if __name__ == "__main__":
    main()
