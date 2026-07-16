"""How many of the covid19 significant alpha-beta pairs survive the depth + HLA ceilings?

theta_depth = 1+CV^2(N) is the lift depth alone induces; (1+CV^2)/f is the composed ceiling for
an allele carried at frequency f. A pair is only "not explicable by depth/shared HLA" if its
observed theta exceeds these. Reports the theta distribution of the significant set against both,
and repeats with a depth floor (junk samples inflate CV). 2026-07-16.
"""
import sys
import numpy as np
import polars as pl
from pathlib import Path
sys.path.insert(0, "appendix")
from bench_biomarker import load_cohort, _glob_tables          # noqa: E402
from vdjtools.io.cohort import SAMPLE_ID                       # noqa: E402
from vdjtools.biomarker import cooccurrence                    # noqa: E402

FMBA = Path("/projects/fmba_covid")
meta = pl.read_csv(FMBA / "metadata_fmba_full.txt", separator="\t", infer_schema_length=0)
cohort = load_cohort(_glob_tables(FMBA / "data", meta, "id")).collect().lazy()

def report(lf, label, min_depth=0):
    d = lf.group_by(SAMPLE_ID).agg(pl.len().alias("n")).collect()
    if min_depth:
        keep = d.filter(pl.col("n") >= min_depth)[SAMPLE_ID].to_list()
        lf = lf.filter(pl.col(SAMPLE_ID).is_in(keep))
        d = d.filter(pl.col("n") >= min_depth)
    N = d["n"].to_numpy().astype(float)
    cv = N.std() / N.mean()
    td = 1 + cv**2
    cc = cooccurrence(lf, chain_a="TRA", chain_b="TRB", min_incidence=10,
                      min_cooccurrence=3, min_incidence_frac=0.03, evalue=True)
    sig = cc.filter(pl.col("q_value") < 0.05)
    th = sig["theta"].to_numpy()
    print(f"\n=== {label}: subjects={len(N)}  CV={cv:.3f}  theta_depth={td:.3f} ===")
    print(f"  tested={cc.height}  significant={sig.height}")
    if not len(th):
        return
    qs = np.percentile(th, [50, 75, 90, 95, 99, 100])
    print(f"  theta of significant: median={qs[0]:.2f} p75={qs[1]:.2f} p90={qs[2]:.2f} "
          f"p95={qs[3]:.2f} p99={qs[4]:.2f} max={qs[5]:.2f}")
    print(f"  above depth ceiling (theta>{td:.2f}): {(th > td).sum()} / {len(th)} "
          f"({(th > td).mean():.1%})")
    # Composed HLA ceilings (1+CV^2)/f for the real carrier frequencies (466 typed subjects).
    for allele, f in (("A*02:01", 0.489), ("A*03:01", 0.273), ("A*01:01", 0.221),
                      ("A*24:02", 0.212), ("B*07:02", 0.202), ("B*08:01", 0.148),
                      ("A*11:01", 0.114)):
        ceil = td / f
        print(f"    vs {allele:8s} f={f:.3f}  ceiling={ceil:5.2f}  "
              f"pairs above: {(th > ceil).sum():4d} ({(th > ceil).mean():5.1%})")

report(cohort, "as shipped (no depth floor)")
report(cohort, "depth floor >= 1000 clonotypes", min_depth=1000)
