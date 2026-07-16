"""Biomarker association / co-occurrence benchmark on large real cohorts.

Runs the full :mod:`vdjtools.biomarker` framework on the three benchmark cohorts and validates
the hits against an external ground truth. Standalone (run via ``python appendix/bench_biomarker.py``,
locally or on Aldan-3 via ``scripts/biomarker_*.sh``); not a pytest test.

Cohorts (``--dataset``):

- ``hip``      — Emerson HIP (786 TCRβ subjects, HF ``isalgo/airr_hip``): CMV association (all five
  tests), per-HLA-allele association, CMV | HLA-A*02 (Cochran–Mantel–Haenszel). Validated vs VDJdb-CMV.
- ``covid19``  — FMBA covid (TRA+TRB; ``--data-dir`` of VDJtools tables + a metadata TSV with
  ``COVID_status`` and HLA): COVID association + **α-β co-occurrence** (in-silico pairing).
  Validated vs the study's published ``covid_associated_clonotypes.csv``.
- ``covid19_vacc`` — FMBA vaccine (TRA+TRB; ``--data-dir`` + metadata with ``timepoint``/``vaccine``):
  before/after-vaccination association.

Every run prints wall time + peak RSS per stage and a hit/validation summary. See ``SOURCES.md``
(Phase 5 / Phase 6) for data provenance and the exact ``aldan3 pull`` staging commands.
"""
from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S
from vdjtools.biomarker import association, condition, cooccurrence

TESTS = ["fisher", "chi2", "bayes_logodds", "bayes_bf", "permutation"]


def _rss_gb() -> float:
    ru = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return ru / 1024**3 if sys.platform == "darwin" else ru / 1024**2   # macOS bytes / Linux KB


class Timer:
    def __init__(self, label):
        self.label = label

    def __enter__(self):
        self.t = time.perf_counter()
        return self

    def __exit__(self, *a):
        print(f"  [{self.label}] {time.perf_counter() - self.t:6.1f}s  peakRSS {_rss_gb():.1f} GB")


# ── loaders ──────────────────────────────────────────────────────────────────────────
def load_cohort(files: "dict[str, str]", fmt: str = "vdjtools") -> pl.LazyFrame:
    """Read a {sample_id: path} map of clonotype tables into one long frame, tagged by sample_id."""
    frames = [vio.read(path, fmt=fmt).with_columns(pl.lit(sid).alias("sample_id"))
              for sid, path in files.items()]
    return pl.concat(frames, how="vertical_relaxed").lazy()


def hip_metadata(meta_txt: Path) -> pl.DataFrame:
    m = pl.read_csv(meta_txt, separator="\t", infer_schema_length=0)   # race has commas → TAB
    return m.rename({"sample_id": "sample_id"}) if "sample_id" in m.columns else m


# ── benchmark bodies ─────────────────────────────────────────────────────────────────
def run_association_suite(cohort, meta, pheno_col, *, key=None, min_incidence=8, label=""):
    key = key or (S.JUNCTION_AA, S.V_CALL, S.J_CALL)
    print(f"\n[{label}] association: {len(TESTS)} tests on '{pheno_col}' "
          f"(key={key}, min_incidence={min_incidence})")
    with Timer("association"):
        res = association(cohort, condition.binary(meta, pheno_col), test=TESTS, key=key,
                          min_incidence=min_incidence, alternative="greater")
    fis = res.filter(pl.col("test") == "fisher").sort("q_value")
    print(f"  features tested: {fis.height}   significant (q<0.05, fisher): "
          f"{fis.filter(pl.col('q_value') < 0.05).height}")
    # cross-test agreement on the top-100 fisher hits
    top = set(fis.head(100)[S.JUNCTION_AA].to_list())
    for t in ("chi2", "permutation", "bayes_bf"):
        sub = res.filter(pl.col("test") == t)
        col = "q_value" if t != "bayes_bf" else "log_bf"
        thr = (pl.col(col) < 0.05) if t != "bayes_bf" else (pl.col(col) > 3)
        hit = set(sub.filter(thr)[S.JUNCTION_AA].to_list())
        print(f"    {t:12s}: {len(top & hit)}/100 of the top fisher hits also flagged")
    return fis


def validate_vs_reference(hits: pl.DataFrame, ref_cdr3: set, name: str, top: int = 200):
    sig = hits.filter(pl.col("q_value") < 0.05)
    ours = set(sig.head(top)[S.JUNCTION_AA].to_list()) if sig.height else set()
    ov = ours & ref_cdr3
    print(f"  vs {name}: {len(ov)}/{len(ours)} of the top significant hits are in the reference "
          f"({len(ref_cdr3)} ref CDR3s)")
    return ov


def run_hip(args):
    from huggingface_hub import snapshot_download
    cache = Path(args.data_dir or "appendix/.data/hip")
    root = Path(snapshot_download("isalgo/airr_hip", repo_type="dataset", local_dir=cache,
                                  allow_patterns=["metadata.txt", "corr/*.txt.gz"]))
    meta = hip_metadata(root / "metadata.txt").filter(pl.col("cmv").is_in(["+", "-"]))
    if args.max_samples:
        meta = meta.head(args.max_samples)
    files = {r["sample_id"]: str(root / "corr" / f"{r['sample_id']}.txt.gz")
             for r in meta.iter_rows(named=True)}
    print(f"hip: {len(files)} subjects")
    with Timer("ingest"):
        cohort = load_cohort(files).collect().lazy()
    cmv_hits = run_association_suite(cohort, meta, "cmv", min_incidence=args.min_incidence, label="CMV")
    # per-HLA-allele + CMH conditioned on HLA-A*02
    hla = meta.with_columns(pl.col("hla").str.contains(r"HLA-A\*02").alias("a02"))
    with Timer("CMH (CMV | HLA-A*02)"):
        cmh = association(cohort, condition.stratified(
            hla.with_columns(pl.col("cmv"), pl.col("a02").cast(pl.Utf8)), "cmv", "a02"),
            stratum_col="_stratum", min_incidence=args.min_incidence)
    print(f"  CMH: {cmh.filter(pl.col('q_value') < 0.05).height} significant (q<0.05)")
    if args.vdjdb:
        ref = _vdjdb_cmv(Path(args.vdjdb))
        validate_vs_reference(cmv_hits, ref, "VDJdb-CMV")


def run_covid(args, pheno_col):
    data = Path(args.data_dir)
    meta = pl.read_csv(args.metadata, separator="\t", infer_schema_length=0)
    files = _glob_tables(data, meta, args.sample_col)
    print(f"{args.dataset}: {len(files)} subjects (TRA+TRB)")
    with Timer("ingest"):
        cohort = load_cohort(files).collect().lazy()
    hits = run_association_suite(cohort, meta, pheno_col, min_incidence=args.min_incidence,
                                 label=pheno_col)
    if args.oracle:
        ref = set(pl.read_csv(args.oracle)["cdr3"].to_list())
        validate_vs_reference(hits, ref, Path(args.oracle).name)
    if args.dataset == "covid19":
        print("\n[α-β co-occurrence]")
        with Timer("cooccurrence"):
            cc = cooccurrence(cohort, chain_a="TRA", chain_b="TRB", min_incidence=args.min_incidence,
                              min_cooccurrence=3, min_incidence_frac=0.03, evalue=True)
        sig = cc.filter(pl.col("q_value") < 0.05)
        print(f"  candidate pairs: {cc.height}   significant α-β pairs (q<0.05): {sig.height}")
        print(sig.select("a_junction_aa", "b_junction_aa", "theta", "q_value", "e_value").head(10))


def _glob_tables(data: Path, meta: pl.DataFrame, sample_col: str) -> "dict[str, str]":
    ids = set(meta[sample_col].to_list())
    files = {}
    for p in sorted(data.glob("*.txt*")) + sorted(data.glob("*.tsv*")):
        sid = p.name.split("_")[0].split(".")[0]
        if sid in ids and sid not in files:
            files[sid] = str(p)
    return files


def _vdjdb_cmv(path: Path) -> set:
    db = pl.read_csv(path, separator="\t", infer_schema_length=0)
    cmv = db.filter((pl.col("gene") == "TRB") & (pl.col("species") == "HomoSapiens")
                    & pl.col("antigen.species").str.contains("CMV"))
    return set(cmv["cdr3"].to_list())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True, choices=["hip", "covid19", "covid19_vacc"])
    ap.add_argument("--data-dir")
    ap.add_argument("--metadata")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--oracle", help="reference clonotype CSV with a 'cdr3' column")
    ap.add_argument("--vdjdb", help="VDJdb slim TSV (hip CMV validation)")
    ap.add_argument("--min-incidence", type=int, default=8)
    ap.add_argument("--max-samples", type=int, default=0)
    args = ap.parse_args()
    print(f"=== biomarker benchmark: {args.dataset} ===")
    t0 = time.perf_counter()
    if args.dataset == "hip":
        run_hip(args)
    elif args.dataset == "covid19":
        run_covid(args, pheno_col="COVID_status")
    else:
        run_covid(args, pheno_col="timepoint")
    print(f"\ntotal {time.perf_counter() - t0:.1f}s  peakRSS {_rss_gb():.1f} GB")


if __name__ == "__main__":
    main()
