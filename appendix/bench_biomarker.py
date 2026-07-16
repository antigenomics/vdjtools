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
def load_cohort(pairs: "list[tuple[str, str]]", fmt: str = "vdjtools") -> pl.LazyFrame:
    """Read (sample_id, path) tables into one long frame tagged by sample_id.

    A list (not a dict) so a subject's TRA **and** TRB tables both load under the same
    ``sample_id`` — required for α-β co-occurrence. Unreadable/malformed tables (a few
    near-empty vaccine repertoires are ragged) are skipped and counted, never silently.
    """
    frames, skipped = [], 0
    for sid, path in pairs:
        try:
            frames.append(vio.read(path, fmt=fmt).with_columns(pl.lit(sid).alias("sample_id")))
        except Exception as e:                         # noqa: BLE001 — report and skip bad tables
            skipped += 1
            if skipped <= 5:
                print(f"  skip unreadable {Path(path).name}: {type(e).__name__}")
    if skipped:
        print(f"  skipped {skipped} unreadable tables (of {len(pairs)})")
    return pl.concat(frames, how="vertical_relaxed").lazy()


def hip_metadata(meta_txt: Path) -> pl.DataFrame:
    m = pl.read_csv(meta_txt, separator="\t", infer_schema_length=0)   # race has commas → TAB
    return m.rename({"sample_id": "sample_id"}) if "sample_id" in m.columns else m


# ── benchmark bodies ─────────────────────────────────────────────────────────────────
def _binary_design(meta, sample_col, col, pos, neg):
    """Design frame with ``_pos``: True if ``col`` ∈ ``pos``, False if ∈ ``neg``, else dropped."""
    e = pl.col(col).cast(pl.String).str.strip_chars()
    return (meta.select(pl.col(sample_col).cast(pl.String).alias("sample_id"),
                        pl.when(e.is_in(pos)).then(True).when(e.is_in(neg)).then(False)
                          .otherwise(None).alias("_pos"))
            .drop_nulls("_pos").unique(subset="sample_id"))


def run_association_suite(cohort, design, *, key=None, min_incidence=8, label=""):
    key = key or (S.JUNCTION_AA, S.V_CALL, S.J_CALL)
    npos = int(design["_pos"].sum())
    print(f"\n[{label}] association: {len(TESTS)} tests, design has {npos} pos / "
          f"{design.height - npos} neg (key={key}, min_incidence={min_incidence})")
    with Timer("association"):
        res = association(cohort, design, test=TESTS, key=key,
                          min_incidence=min_incidence, alternative="greater")
    fis = res.filter(pl.col("test") == "fisher").sort("q_value")
    # The DESIGN size is not the ANALYSED size: association() inner-joins the design to the
    # cohort, so subjects with metadata but no repertoire (or vice versa) silently drop out.
    # n_pos/n_neg on the result are the arms the test actually ran on — report those.
    if fis.height:
        a_pos, a_neg = int(fis["n_pos"][0]), int(fis["n_neg"][0])
        print(f"  ANALYSED arms (design ∩ cohort): {a_pos} pos / {a_neg} neg"
              + (f"   ⚠ imbalanced {a_pos / max(a_neg, 1):.0f}:1" if a_pos > 5 * a_neg or
                 a_neg > 5 * a_pos else "")
              + (f"   ⚠ only {min(a_pos, a_neg)} in the smaller arm"
                 if min(a_pos, a_neg) < 50 else ""))
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
    """Overlap **and** enrichment of significant hits with a VDJdb reference CDR3 set.

    The meaningful validation is a 2×2 over the *unique tested CDR3s* — are reference
    (antigen-specific) CDR3s enriched among the significant hits vs the rest? — plus the
    raw overlap counts (all significant, and the top ``top`` by q).
    """
    from scipy.stats import fisher_exact

    tested = set(hits[S.JUNCTION_AA].to_list())
    sig_df = hits.filter(pl.col("q_value") < 0.05)
    sig = set(sig_df[S.JUNCTION_AA].to_list())
    top_ov = set(sig_df.sort("q_value").head(top)[S.JUNCTION_AA].to_list()) & ref_cdr3
    a = len(sig & ref_cdr3)                       # significant & in reference
    b = len(sig) - a                             # significant & not in reference
    c = len((tested - sig) & ref_cdr3)           # not significant & in reference
    d = len(tested - sig) - c                    # not significant & not in reference
    orr, p = fisher_exact([[a, b], [c, d]], alternative="greater") if (a + c) and (b + d) else (float("nan"), 1.0)
    print(f"  vs {name}: {a}/{len(sig)} significant CDR3s in reference ({len(top_ov)} in top {top}); "
          f"{a + c} of {len(tested)} tested CDR3s are reference members; "
          f"enrichment OR={orr:.2f} p={p:.1e}  ({len(ref_cdr3)} ref CDR3s)")
    return sig & ref_cdr3


def _corr_path(root: Path, sid: str) -> "str | None":
    for ext in (".txt.gz", ".txt"):                    # HF hip is gzipped; the cluster copy is not
        p = root / "corr" / f"{sid}{ext}"
        if p.exists():
            return str(p)
    return None


def run_hip(args):
    root = (Path(args.data_dir) if args.data_dir and (Path(args.data_dir) / "metadata.txt").exists()
            else None)
    if root is None:
        from huggingface_hub import snapshot_download
        root = Path(snapshot_download("isalgo/airr_hip", repo_type="dataset",
                                      local_dir=args.data_dir or "appendix/.data/hip",
                                      allow_patterns=["metadata.txt", "corr/*.txt.gz"]))
    meta = hip_metadata(root / "metadata.txt").filter(pl.col("cmv").is_in(["+", "-"]))
    if args.max_samples:
        meta = meta.head(args.max_samples)
    pairs = [(r["sample_id"], _corr_path(root, r["sample_id"])) for r in meta.iter_rows(named=True)]
    pairs = [(s, p) for s, p in pairs if p]
    print(f"hip: {len(pairs)} subjects")
    with Timer("ingest"):
        cohort = load_cohort(pairs).collect().lazy()
    cmv = _binary_design(meta, "sample_id", "cmv", ["+"], ["-"])
    cmv_hits = run_association_suite(cohort, cmv, min_incidence=args.min_incidence, label="CMV")
    # CMH: CMV association stratified by HLA-A*02 carriage
    hla = meta.with_columns(pl.col("hla").fill_null("").str.contains(r"HLA-A\*02").alias("a02"))
    with Timer("CMH (CMV | HLA-A*02)"):
        cmh = association(cohort, condition.stratified(hla, "cmv", "a02"),
                          stratum_col="_stratum", min_incidence=args.min_incidence)
    print(f"  CMH: {cmh.filter(pl.col('q_value') < 0.05).height} significant (q<0.05)")
    if args.vdjdb:
        validate_vs_reference(cmv_hits, _vdjdb_cmv(Path(args.vdjdb)), "VDJdb-CMV")


def run_covid(args):
    data = Path(args.data_dir)
    meta = pl.read_csv(args.metadata, separator="\t", infer_schema_length=0)
    chains = tuple(args.chains.split(",")) if args.chains else None
    pairs = _tables_from_metadata(data, meta, file_col=args.file_col,
                                  sample_col=args.sample_col, locus_col=args.locus_col,
                                  chains=chains)
    print(f"{args.dataset}: {len({s for s, _ in pairs})} subjects, {len(pairs)} tables "
          f"({args.chains or 'all chains'})")
    with Timer("ingest"):
        cohort = load_cohort(pairs).collect().lazy()
    design = _binary_design(meta, args.sample_col, args.pheno_col,
                            args.pos_values.split(","), args.neg_values.split(","))
    hits = run_association_suite(cohort, design, min_incidence=args.min_incidence, label=args.pheno_col)
    if args.vdjdb:
        # Oracle = VDJdb SARS-CoV-2 CDR3s (both chains); is that antigen-specific set enriched
        # among the significant COVID-associated hits? (parallels the hip VDJdb-CMV validation.)
        validate_vs_reference(hits, _vdjdb_antigen(Path(args.vdjdb), "SARS-CoV-2"), "VDJdb-SARS-CoV-2")
    if args.oracle:
        ref = set(pl.read_csv(args.oracle, infer_schema_length=0)["cdr3"].to_list())
        validate_vs_reference(hits, ref, Path(args.oracle).name)
    if args.dataset == "covid19":
        print("\n[α-β co-occurrence]")
        with Timer("cooccurrence"):
            cc = cooccurrence(cohort, chain_a="TRA", chain_b="TRB", min_incidence=args.min_incidence,
                              min_cooccurrence=3, min_incidence_frac=0.03, evalue=True)
        sig = cc.filter(pl.col("q_value") < 0.05)
        print(f"  candidate pairs: {cc.height}   significant α-β pairs (q<0.05): {sig.height}")
        print(sig.select("a_junction_aa", "b_junction_aa", "theta", "q_value", "e_value").head(10))
        if args.vdjdb and cc.height:
            # Do the in-silico pairs recover VDJdb's true (single-cell) α-β pairs? Antigen-agnostic:
            # pairing is a receptor property, not a condition — any VDJdb complex is a valid pair.
            validate_pairs(cc, _vdjdb_pairs(Path(args.vdjdb)), "VDJdb α-β pairs (any antigen)")


def _tables_from_metadata(data: Path, meta: pl.DataFrame, *, file_col: str = "file_name",
                          sample_col: str = "sample_id",
                          locus_col: "str | None" = "locus",
                          chains: "tuple[str, ...] | None" = None) -> "list[tuple[str, str]]":
    """(sample_id, path) driven by the metadata's own ``file_name -> sample_id`` map.

    NEVER derive the sample id from the filename. The old ``p.name.split("_")[0]`` collapsed every
    ``p17_*_DNA_*`` file to the id ``"p17"`` — merging 73 distinct donors into one pseudo-subject
    (and 32 more into ``"p18"``), i.e. 1,157 ids instead of 1,258, with two "donors" holding 105
    repertoires between them. It only looked correct on the all-12-digit cluster directory.
    The HF in-repo sheet maps file->sample exactly (2,516/2,516), so use it.
    """
    m = meta
    if chains and locus_col and locus_col in m.columns:
        m = m.filter(pl.col(locus_col).is_in(list(chains)))
    pairs, missing = [], 0
    for r in m.iter_rows(named=True):
        fn = str(r[file_col])
        p = data / fn
        if not p.exists():                       # HF gz vs the sheet's legacy .txt name
            stem = fn.split(".")[0]
            hits = sorted(data.glob(f"{stem}.*"))
            if locus_col and locus_col in m.columns:
                hits = [h for h in hits if f".{r[locus_col]}." in h.name]
            p = hits[0] if hits else None
        if p is None or not p.exists():
            missing += 1
            continue
        pairs.append((str(r[sample_col]), str(p)))
    if missing:
        print(f"  {missing} metadata rows have no file on disk", flush=True)
    return pairs


def _vdjdb_antigen(path: Path, pattern: "str | None" = None, gene: "str | None" = None,
                   mhc: "str | None" = None) -> set:
    """Human VDJdb CDR3s, scoped to the question being asked.

    - **antigen-specific association** (infection status): ``pattern="SARS-CoV-2"`` / ``"CMV"``.
    - **HLA association**: ``mhc="HLA-A*02"`` — *any* epitope restricted by that allele; the
      screen tests HLA carriage, not a particular antigen, so the oracle must not be
      antigen-restricted.
    - ``gene=None`` keeps both chains — for the covid TRA+TRB cohort a hit's ``junction_aa``
      may be either chain (VDJdb SARS-CoV-2: 3,796 TRA + 5,333 TRB human records).
    """
    db = pl.read_csv(path, separator="\t", infer_schema_length=0)
    f = pl.col("species") == "HomoSapiens"
    if pattern:
        f = f & pl.col("antigen.species").str.contains(pattern)
    if gene:
        f = f & (pl.col("gene") == gene)
    if mhc:
        f = f & pl.col("mhc.a").str.starts_with(mhc)
    return set(db.filter(f)["cdr3"].to_list())


def _vdjdb_cmv(path: Path) -> set:
    return _vdjdb_antigen(path, "CMV", gene="TRB")


def _vdjdb_pairs(path: Path, pattern: "str | None" = None, mhc: "str | None" = None) -> set:
    """Known ``(cdr3_alpha, cdr3_beta)`` pairs from VDJdb complexes.

    A ``complex.id`` groups the TRA and TRB of one single-cell-resolved clonotype, so these are
    *true* α-β pairs. Scope the oracle to the question being asked:

    - **α-β pairing** (co-occurrence): ``pattern=None`` — **any antigen** (21,679 human complexes).
      Whether two chains pair is a property of the receptor, not of a condition, so requiring the
      pair to be e.g. COVID-associated would needlessly shrink the oracle (3,031) and understate
      the enrichment.
    - **HLA association**: ``mhc="HLA-A*02"`` — any epitope restricted by that allele.
    - **Antigen-specific association** (COVID status): ``pattern="SARS-CoV-2"``.
    """
    db = pl.read_csv(path, separator="\t", infer_schema_length=0)
    f = pl.col("species") == "HomoSapiens"
    if pattern:
        f = f & pl.col("antigen.species").str.contains(pattern)
    if mhc:
        f = f & pl.col("mhc.a").str.starts_with(mhc)
    # `complex.id` in the SLIM dump is a COMMA-SEPARATED LIST — slim collapses duplicate records,
    # so one row can belong to many complexes (7,087 rows do). Joining on the raw string pairs two
    # records only when their entire list is byte-identical, which silently drops most real pairs:
    # measured 20,169 vs 60,798 exploded. Every earlier alpha-beta validation used the short oracle.
    d = (db.filter(f)
         .with_columns(pl.col("complex.id").str.split(",").alias("_cid"))
         .explode("_cid")
         .with_columns(pl.col("_cid").str.strip_chars())
         .filter(pl.col("_cid") != "0"))
    a = d.filter(pl.col("gene") == "TRA").select("_cid", pl.col("cdr3").alias("a"))
    b = d.filter(pl.col("gene") == "TRB").select("_cid", pl.col("cdr3").alias("b"))
    j = a.join(b, on="_cid", how="inner").unique(subset=["a", "b"])
    return set(zip(j["a"].to_list(), j["b"].to_list()))


def validate_pairs(cc: pl.DataFrame, ref_pairs: set, name: str) -> None:
    """Enrichment of known VDJdb α-β pairs among the significant co-occurring pairs.

    2×2 over the *tested* pairs: significant × known. An honest zero is reported as such —
    the tested set is capped by ``max_features`` and skewed to high-incidence public clones,
    which need not be the single-cell-characterised ones in VDJdb.
    """
    from scipy.stats import fisher_exact

    tested = set(zip(cc["a_junction_aa"].to_list(), cc["b_junction_aa"].to_list()))
    sig_df = cc.filter(pl.col("q_value") < 0.05)
    sig = set(zip(sig_df["a_junction_aa"].to_list(), sig_df["b_junction_aa"].to_list()))
    a = len(sig & ref_pairs)
    b = len(sig) - a
    c = len((tested - sig) & ref_pairs)
    d = len(tested - sig) - c
    orr, p = (fisher_exact([[a, b], [c, d]], alternative="greater")
              if (a + c) and (b + d) else (float("nan"), 1.0))
    print(f"  vs {name}: {a}/{len(sig)} significant α-β pairs are known VDJdb pairs; "
          f"{a + c} of {len(tested)} tested pairs are known; enrichment OR={orr:.2f} p={p:.1e} "
          f"({len(ref_pairs)} ref pairs)")


#: Per-dataset covid phenotype defaults (positive / negative metadata values).
_COVID_DEFAULTS = {
    "covid19": ("sample.COVID_status", "current,past", "healthy"),
    "covid19_vacc": ("timepoint", "20d_after_vaccination", "before_vaccination"),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True, choices=["hip", "covid19", "covid19_vacc"])
    ap.add_argument("--data-dir")
    ap.add_argument("--metadata")
    ap.add_argument("--sample-col", default="sample_id")
    ap.add_argument("--pheno-col", help="metadata phenotype column (covid; per-dataset default)")
    ap.add_argument("--pos-values", help="comma list of positive phenotype values (covid)")
    ap.add_argument("--neg-values", help="comma list of negative phenotype values (covid)")
    ap.add_argument("--oracle", help="reference clonotype CSV with a 'cdr3' column")
    ap.add_argument("--vdjdb", help="VDJdb slim TSV (CMV / SARS-CoV-2 validation)")
    ap.add_argument("--file-col", default="file_name", help="metadata column holding the filename")
    ap.add_argument("--locus-col", default="locus", help="metadata column holding the chain")
    ap.add_argument("--chains", default=None, help="comma list, e.g. TRA or TRB or TRA,TRB")
    ap.add_argument("--min-incidence", type=int, default=8)
    ap.add_argument("--max-samples", type=int, default=0)
    args = ap.parse_args()
    print(f"=== biomarker benchmark: {args.dataset} ===")
    t0 = time.perf_counter()
    if args.dataset == "hip":
        run_hip(args)
    else:
        pc, pos, neg = _COVID_DEFAULTS[args.dataset]
        args.pheno_col = args.pheno_col or pc
        args.pos_values = args.pos_values or pos
        args.neg_values = args.neg_values or neg
        run_covid(args)
    print(f"\ntotal {time.perf_counter() - t0:.1f}s  peakRSS {_rss_gb():.1f} GB")


if __name__ == "__main__":
    main()
