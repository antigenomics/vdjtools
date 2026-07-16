"""Is the covid19 ASSOCIATION confounded by depth? (84% of features significant is a red flag.)

Depth inflates presence of EVERY feature in deep subjects. If depth correlates with COVID status,
every feature looks enriched in the deeper arm. Checks (1) depth vs status, and (2) whether a
depth-stratified CMH association survives. 2026-07-16.
"""
import sys
import numpy as np
import polars as pl
from pathlib import Path
sys.path.insert(0, "appendix")
from bench_biomarker import load_cohort, _glob_tables, _binary_design   # noqa: E402
from vdjtools.io.cohort import SAMPLE_ID                                # noqa: E402
from vdjtools.biomarker import association                              # noqa: E402
from scipy.stats import mannwhitneyu                                    # noqa: E402

FMBA = Path("/projects/fmba_covid")
meta = pl.read_csv(FMBA / "metadata_fmba_full.txt", separator="\t", infer_schema_length=0)
cohort = load_cohort(_glob_tables(FMBA / "data", meta, "id")).collect().lazy()
design = _binary_design(meta, "id", "sample.COVID_status", ["current", "past"], ["healthy"])

d = cohort.group_by(SAMPLE_ID).agg(pl.len().alias("S")).collect()
j = design.join(d, on=SAMPLE_ID, how="inner")
pos = j.filter(pl.col("_pos"))["S"].to_numpy().astype(float)
neg = j.filter(~pl.col("_pos"))["S"].to_numpy().astype(float)
u, p = mannwhitneyu(pos, neg, alternative="two-sided")
print(f"depth by COVID status: pos n={len(pos)} median={np.median(pos):,.0f} | "
      f"neg n={len(neg)} median={np.median(neg):,.0f}")
print(f"  ratio of medians = {np.median(pos)/np.median(neg):.2f}x   Mann-Whitney p={p:.3g}")
print(f"  => depth {'IS' if p < 0.05 else 'is NOT'} associated with the phenotype")

# Depth-stratified association: same design + a depth-decile stratum -> CMH.
q = np.quantile(j["S"].to_numpy().astype(float), np.linspace(0, 1, 11)[1:-1])
strat = j.with_columns(
    pl.col("S").map_elements(lambda s: str(int(np.searchsorted(q, s, side="right"))),
                             return_dtype=pl.String).alias("_stratum")
).select(SAMPLE_ID, "_pos", "_stratum")
print("\nassociation, depth-stratified (CMH over 10 depth deciles):")
r = association(cohort, strat, min_incidence=10, alternative="greater")
sig = r.filter(pl.col("q_value") < 0.05)
print(f"  tested={r.height}  significant={sig.height}  ({sig.height/max(r.height,1):.1%})")
print("  reference: unstratified Fisher gave 44,125 / 52,528 = 84.0%")
