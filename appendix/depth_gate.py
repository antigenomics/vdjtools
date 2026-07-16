"""Depth gate: is the covid19 alpha-beta co-occurrence significance depth-driven?

theta_depth = 1 + CV^2(N) is the co-occurrence lift induced by per-subject repertoire-size
variation ALONE -- independent of HLA, exposure, or physical pairing (De Witt 2018 folds the
same effect in as a per-subject bias factor b_i). If theta_depth is comparable to the observed
theta, the shipped Fisher q-values are not calibrated. 2026-07-16.
"""
import sys
import numpy as np
import polars as pl
sys.path.insert(0, "appendix")
from bench_biomarker import load_cohort, _glob_tables          # noqa: E402
from vdjtools.io import schema as S                            # noqa: E402
from vdjtools.io.cohort import SAMPLE_ID                       # noqa: E402
from vdjtools.biomarker.association import select_candidates   # noqa: E402
from pathlib import Path                                       # noqa: E402

FMBA = Path("/projects/fmba_covid")
meta = pl.read_csv(FMBA / "metadata_fmba_full.txt", separator="\t", infer_schema_length=0)
pairs = _glob_tables(FMBA / "data", meta, "id")
cohort = load_cohort(pairs).collect().lazy()

# 1) Depth per subject = unique clonotypes (the unit incidence is counted in).
depth = cohort.group_by(SAMPLE_ID).agg(pl.len().alias("n")).collect()["n"].to_numpy().astype(float)
cv = depth.std() / depth.mean()
sd_log = np.log(depth).std()
theta_depth_rare = 1.0 + cv**2
print(f"subjects={len(depth)}  median depth={np.median(depth):,.0f}  "
      f"min={depth.min():,.0f}  max={depth.max():,.0f}  spread={depth.max()/depth.min():.1f}x")
print(f"CV(N)={cv:.3f}  sd(log N)={sd_log:.3f}  -> theta_depth (rare-clone bound) = {theta_depth_rare:.3f}")

# 2) The tested candidates are top-incidence (max_features cap), where theta_depth attenuates.
#    Report the incidence distribution of what actually gets tested per chain.
for chain in ("TRA", "TRB"):
    sub = cohort.filter(pl.col(S.V_CALL).str.slice(0, 3) == chain)
    cand = select_candidates(sub, min_incidence=10, min_incidence_frac=0.03).head(2000)
    inc = cand["incidence"].to_numpy()
    n_sub = len(depth)
    print(f"{chain}: candidates={cand.height}  incidence median={np.median(inc):.0f} "
          f"({np.median(inc)/n_sub:.1%} of subjects)  min={inc.min()}  max={inc.max()}")

# 3) Attenuated theta_depth: for a clone with detection prob p_i = 1-exp(-N_i*pi), the induced
#    lift is E[p^2]/E[p]^2 evaluated at the pi matching the observed median incidence.
w = depth / depth.mean()
for label, target in (("median-incidence candidate", 0.15), ("rare (min_incidence=10)", 10/len(depth))):
    lo, hi = 1e-12, 1.0
    for _ in range(200):                     # solve pi s.t. mean detection == target
        mid = (lo + hi) / 2
        if (1 - np.exp(-w * mid * depth.mean())).mean() < target:
            lo = mid
        else:
            hi = mid
    p = 1 - np.exp(-w * mid * depth.mean())
    print(f"theta_depth @ {label} (mean detect={p.mean():.3f}): {(p**2).mean()/p.mean()**2:.3f}")
