"""Acceptance gate for depth_strata: does CMH on the FULL cohort land near the depth-floored 754?

If depth-stratified CMH on all 552 subjects gives a hit count near the 754 that the >=1000-clonotype
floor produced (rather than near the uncorrected 2802), the simulation transfers to real data and
the floor becomes a diagnostic rather than required preprocessing. If it stays near 2802, reopen.
2026-07-16.
"""
import sys
import numpy as np
import polars as pl
from pathlib import Path
sys.path.insert(0, "appendix")
from bench_biomarker import load_cohort, _glob_tables            # noqa: E402
from vdjtools.biomarker import cooccurrence                      # noqa: E402

FMBA = Path("/projects/fmba_covid")
meta = pl.read_csv(FMBA / "metadata_fmba_full.txt", separator="\t", infer_schema_length=0)
cohort = load_cohort(_glob_tables(FMBA / "data", meta, "id")).collect().lazy()
kw = dict(chain_a="TRA", chain_b="TRB", min_incidence=10, min_cooccurrence=3,
          min_incidence_frac=0.03, evalue=True)

for label, ds in (("pooled Fisher (depth_strata=0, old default)", 0),
                  ("CMH over depth strata (depth_strata=10, new default)", 10)):
    import time; t = time.perf_counter()
    cc = cooccurrence(cohort, depth_strata=ds, **kw)
    sig = cc.filter(pl.col("q_value") < 0.05)
    th = sig["theta"].to_numpy()
    extra = ""
    if "or_mh" in cc.columns and sig.height:
        extra = f"  median or_mh={np.median(sig['or_mh'].to_numpy()):.2f}"
    print(f"{label}\n  tested={cc.height}  significant={sig.height}  "
          f"theta median={np.median(th):.2f} max={th.max():.2f}{extra}  [{time.perf_counter()-t:.0f}s]"
          if sig.height else f"{label}\n  tested={cc.height}  significant=0")
print("\nreference: >=1000-clonotype depth floor + pooled Fisher gave 754 significant (of 733,178)")
