"""Measure the false-positive rate of the co-occurrence tests under a depth-only null.

This is the provenance for the `depth_strata=10` default (see SOURCES.md "Co-occurrence
confounding"). Two clonotypes are drawn INDEPENDENTLY given each subject's repertoire size, so
every significant pair is a false positive. Depth is the only thing linking them.

The three statistics compared are the ones actually on the table:
  1. pooled Fisher            — the pre-2.7.0 default
  2. depth-weighted incidence — "normalized incidence" / "found together in x of X rearrangements"
  3. CMH over depth strata    — the 2.7.0 default

A calibrated test rejects at ~the nominal rate. Run: python appendix/cooccurrence_fpr.py
Takes ~1 min. 2026-07-16.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import chi2 as _chi2
from scipy.stats import hypergeom

SEED = 20260716
N_SUBJECTS = 552          # FMBA covid19
SD_LOG_DEPTH = 0.80       # tuned to CV ~= 0.9, the measured cohort value
N_PAIRS = 2000
INCIDENCES = (0.02, 0.03, 0.11)   # rare -> the ~11% regime max_features steers callers into


def fisher_greater(a, b, c, d):
    """One-tailed hypergeometric upper tail P(X >= a) — same convention as stats.fisher_p."""
    n = a + b + c + d
    return hypergeom.sf(a - 1, n, a + b, a + c)


def cmh(A, B, C, D):
    """Two-sided CMH chi2 over strata (rows=pairs, cols=strata) — mirrors stats.cmh."""
    n = A + B + C + D
    ok = n > 0
    nz = np.where(ok, n, 1.0)
    obs = np.where(ok, A, 0.0).sum(1)
    exp = np.where(ok, (A + B) * (A + C) / nz, 0.0).sum(1)
    var = np.where(ok & (n > 1), (A + B) * (C + D) * (A + C) * (B + D) / (nz**2 * (nz - 1.0)),
                   0.0).sum(1)
    stat = np.where(var > 0, (np.abs(obs - exp) - 0.5) ** 2 / np.where(var > 0, var, 1.0), 0.0)
    return _chi2.sf(stat, 1)


def run(rng, pi, S, label):
    n = S.size
    p_det = 1 - np.exp(-S * pi)                       # Poisson detection at depth S
    A = rng.random((N_PAIRS, n)) < p_det              # independent given depth
    B = rng.random((N_PAIRS, n)) < p_det

    n_a, n_b = A.sum(1), B.sum(1)
    n_ab = (A & B).sum(1)
    a = n_ab.astype(float)
    b, c = n_a - n_ab, n_b - n_ab
    d = n - n_a - n_b + n_ab
    keep = (n_a > 0) & (n_b > 0) & (d > 0)            # testable pairs only

    # 1. pooled Fisher
    p_fisher = fisher_greater(a, b, c, d)

    # 2. depth-weighted ("normalized") incidence: subject COUNTS -> sums of subject depths.
    #    NB this is NOT a well-defined test, which is the point. A hypergeometric is a statement
    #    about drawing discrete units; feeding it sums of depths silently asserts sum(S)~7e6
    #    exchangeable pseudo-observations instead of 552 subjects. The resulting "p" is an
    #    artifact of the variant you happen to write: this formulation degenerates conservative
    #    (FPR -> 0), while scaling the weights to keep n fixed degenerates anticonservative
    #    (FPR -> 0.9). Reported to show it is unusable in EITHER direction, not to pin a rate.
    wa = (A * S).sum(1)
    wb = (B * S).sum(1)
    wab = ((A & B) * S).sum(1)
    W = S.sum()
    p_weighted = fisher_greater(wab, wa - wab, wb - wab, W - wa - wb + wab)

    # 3. CMH over 10 equal-count depth strata
    edges = np.quantile(S, np.linspace(0, 1, 11)[1:-1])
    strata = np.searchsorted(edges, S, side="right")
    ns = strata.max() + 1
    SA = np.empty((N_PAIRS, ns))
    SB, SC, SD = (np.empty_like(SA) for _ in range(3))
    for k in range(ns):
        m = strata == k
        ak = (A[:, m] & B[:, m]).sum(1)
        nak, nbk = A[:, m].sum(1), B[:, m].sum(1)
        SA[:, k], SB[:, k], SC[:, k] = ak, nak - ak, nbk - ak
        SD[:, k] = m.sum() - nak - nbk + ak
    p_cmh = cmh(SA, SB, SC, SD)

    row = {"pooled Fisher": p_fisher, "depth-weighted": p_weighted, "CMH/depth strata": p_cmh}
    print(f"\n  {label}: mean incidence {p_det.mean():.3f}, testable {keep.sum()}/{N_PAIRS}")
    for name, p in row.items():
        pk = p[keep]
        print(f"    {name:20s} FPR@0.05 = {np.mean(pk < 0.05):.4f}   "
              f"FPR@1e-6 = {np.mean(pk < 1e-6):.5f}")


def main():
    rng = np.random.default_rng(SEED)
    print(f"seed={SEED}  n_subjects={N_SUBJECTS}  pairs/config={N_PAIRS}")
    S = np.exp(rng.normal(np.log(14060), SD_LOG_DEPTH, N_SUBJECTS))
    print(f"depth: median={np.median(S):,.0f} CV={S.std()/S.mean():.3f} "
          f"-> theta_depth = 1+CV^2 = {1 + (S.std()/S.mean())**2:.3f}")
    print("A calibrated test gives FPR@0.05 ~ 0.05 and FPR@1e-6 ~ 1e-6.")
    for target in INCIDENCES:
        pi = -np.log(1 - target) / S.mean()           # detection ~= target at the mean depth
        run(rng, pi, S, f"incidence ~{target:.0%}")


if __name__ == "__main__":
    main()
