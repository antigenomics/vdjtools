"""The IMPRINT catalog — HLA-restricted clonotype families and what (if anything) they are.

De Witt 2018: a cohort of repertoires always carries imprints of shared exposure, whether or not
the phenotype you are testing is one of them. So "no COVID-specific hit" is not a null result --
the HLA-restricted families are there regardless, and they are the finding. The question is only
which exposure each one records.

VDJdb can answer that for MHC class I (108,551 human records: CMV 33k, EBV 11k, InfluenzaA 11k,
SARS-CoV-2 8k) but barely for class II (4,076 records total -- 26.6x less). So:

  * class I  family + VDJdb hit  -> a known imprint (EBV / CMV / flu / SARS-CoV-2). Confirms the
    pipeline and dates the exposure.
  * class II family + NO VDJdb hit -> a REAL imprint of a common exposure that VDJdb does not
    cover. Unannotatable by lookup, and therefore exactly what an embedding-based method (mirpy)
    exists to characterise. These are the deliverable, not the leftovers.

For each significant biomarker family we report: the HLA it is restricted BY (measured within
cases, so the phenotype is not the confound), the MHC class implied by that allele, and every
VDJdb species/epitope within Hamming-1 of the CDR3 (1mm, not exact -- VDJdb is sparse and an
exact lookup understates coverage, which would inflate the "novel" call).

Ranked output feeds /Users/mikesh/vcs/projects/2026-mirpy-analysis/plans. 2026-07-17.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import polars as pl

RES = Path("/projects/biomarkers/results")
VDJDB = Path("/projects/immunestatus/vdjdb/vdjdb-2025-07-30/vdjdb.slim.txt")

# an allele name tells you the class; our measured restriction is the ground truth here
CLASS_I = ("A*", "B*", "C*")
CLASS_II = ("DRB", "DQB", "DPB", "DQA", "DPA", "DRA")


def mhc_class(allele: str) -> str:
    if allele.startswith(CLASS_I):
        return "I"
    if allele.startswith(CLASS_II):
        return "II"
    return "?"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", default="TRA", choices=["TRA", "TRB"])
    args = ap.parse_args()
    gene = args.chain
    print(f"=== IMPRINT catalog {gene} ===", flush=True)

    mh = pl.read_parquet(RES / f"motif_hla_{gene}.parquet")
    print(f"  motifs with a measured HLA restriction: {mh.height}", flush=True)

    db = pl.read_csv(VDJDB, separator="\t", infer_schema_length=0)
    v = (db.filter((pl.col("gene") == gene) & (pl.col("species") == "HomoSapiens"))
         .select("cdr3", "antigen.species", "antigen.epitope", "mhc.a", "mhc.class").unique())
    print(f"  VDJdb {gene} human records: {v.height:,} "
          f"(MHCI {v.filter(pl.col('mhc.class')=='MHCI').height:,} / "
          f"MHCII {v.filter(pl.col('mhc.class')=='MHCII').height:,})", flush=True)

    # 1mm match OUR motifs -> VDJdb. Exact lookup understates VDJdb coverage and would
    # over-call "novel"; vdjmatch makes the sensitive version cheap.
    import vdjmatch.cluster as vc
    ours = mh["cand"].unique().sort().to_list()
    vd = v["cdr3"].unique().sort().to_list()
    pairs = vc.overlap(ours, vd, scope="1,0,0,1", threads=0)
    qm = pl.DataFrame({"a_idx": np.arange(len(ours), dtype=np.int64), "cand": ours})
    um = pl.DataFrame({"b_idx": np.arange(len(vd), dtype=np.int64), "cdr3": vd})
    ann = (pairs.join(qm, on="a_idx").join(um, on="b_idx").join(v, on="cdr3", how="inner")
           .select("cand", "antigen.species", "antigen.epitope", "mhc.a", "mhc.class").unique())
    print(f"  1mm VDJdb annotations for our motifs: {ann.height:,} "
          f"(covering {ann['cand'].n_unique()} of {len(ours)} motifs)", flush=True)

    per = (ann.group_by("cand").agg(
        pl.col("antigen.species").unique().alias("vdjdb_species"),
        pl.col("antigen.epitope").unique().alias("vdjdb_epitopes"),
        pl.col("mhc.class").unique().alias("vdjdb_class")))
    cat = (mh.join(per, on="cand", how="left")
           .with_columns(pl.col("allele").map_elements(mhc_class, return_dtype=pl.String)
                         .alias("our_class"),
                         pl.col("vdjdb_species").list.len().fill_null(0).alias("n_species")))

    print("\n  --- families by OUR measured MHC class x VDJdb annotation ---", flush=True)
    summ = (cat.with_columns((pl.col("n_species") > 0).alias("annotated"))
            .group_by(["our_class", "annotated"]).agg(pl.len().alias("n_motifs"))
            .sort(["our_class", "annotated"]))
    print(summ, flush=True)

    print("\n  === KNOWN IMPRINTS (VDJdb-annotated) ===", flush=True)
    known = cat.filter(pl.col("n_species") > 0).sort("p")
    for r in known.head(14).iter_rows(named=True):
        sp = ",".join(sorted(r["vdjdb_species"] or []))
        ep = ",".join(sorted(r["vdjdb_epitopes"] or [])[:2])
        print(f"    {r['cand'][:20]:20s} {r['v'][:10]:10s} -> {r['allele']:11s} "
              f"(class {r['our_class']}) OR={r['or']:7.1f}  VDJdb: {sp:26s} {ep}", flush=True)

    print("\n  === NOVEL IMPRINTS (HLA-restricted, NO VDJdb match at 1mm) ===", flush=True)
    print("      -> these are the mirpy targets: real, reproducible, unannotatable by lookup",
          flush=True)
    novel = cat.filter(pl.col("n_species") == 0).sort("p")
    for cls in ("II", "I"):
        sub = novel.filter(pl.col("our_class") == cls)
        print(f"\n    -- class {cls}: {sub.height} motifs", flush=True)
        for r in sub.head(12).iter_rows(named=True):
            print(f"    {r['cand'][:20]:20s} {r['v'][:10]:10s} -> {r['allele']:11s} "
                  f"OR={r['or']:8.1f}  p={r['p']:9.2e}  nCOVID={r['n']}", flush=True)

    # which alleles carry the novel class II signal -> the exposure axes to chase in mirpy
    nov2 = novel.filter(pl.col("our_class") == "II")
    if nov2.height:
        print("\n    novel class-II signal by allele:", flush=True)
        print(nov2.group_by("allele").agg(pl.len().alias("n_motifs"),
                                          pl.col("or").max().alias("max_or"))
              .sort("n_motifs", descending=True), flush=True)

    cat.write_parquet(RES / f"imprint_{gene}.parquet")
    print(f"\n-> {RES}/imprint_{gene}.parquet", flush=True)


if __name__ == "__main__":
    main()
