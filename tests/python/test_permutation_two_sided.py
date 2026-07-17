"""association(test='permutation', alternative='two-sided') must be genuinely two-sided.

The bug: it silently substituted the one-sided upper tail, so a feature DEPLETED in cases got a
p ~ 1.0 (upper tail) instead of its real (small) two-sided p.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from vdjtools.biomarker import association, condition


def _clono(sid, v, j, aa):
    return dict(sample_id=sid, v_call=v, j_call=j, junction_aa=aa,
                junction_nt="ACG", duplicate_count=1, frequency=0.0)


def _cohort(seed=1):
    """A feature strongly DEPLETED in the condition-positive arm."""
    rng = np.random.default_rng(seed)
    rows, meta = [], []
    for i in range(60):
        pos = i < 30
        sid = f"S{i:02d}"
        meta.append(dict(sample_id=sid, cmv="+" if pos else "-"))
        rows.append(_clono(sid, "TRBV1", "TRBJ1", "CASSBG"))                       # ubiquitous
        if (not pos and rng.random() < 0.85) or (pos and rng.random() < 0.1):
            rows.append(_clono(sid, "TRBV4", "TRBJ4", "CASSZF"))                   # depleted in +
    return pl.DataFrame(rows), pl.DataFrame(meta)


def test_two_sided_permutation_is_not_the_upper_tail():
    cohort, meta = _cohort()
    ph = condition.binary(meta, "cmv")
    up = association(cohort, ph, test="permutation", alternative="greater", n_perm=2000, seed=0)
    two = association(cohort, ph, test="permutation", alternative="two-sided", n_perm=2000, seed=0)

    z_up = up.filter(pl.col("junction_aa") == "CASSZF")["p_value"][0]
    z_two = two.filter(pl.col("junction_aa") == "CASSZF")["p_value"][0]
    # Depleted-in-cases: upper tail large (~1), two-sided small and clearly different.
    assert z_up > 0.5, f"upper-tail p for a depleted feature should be large, got {z_up}"
    assert z_two < 0.1, f"two-sided p should catch the depletion, got {z_two}"
    assert z_two < z_up


def test_two_sided_permutation_p_is_bounded():
    """The doubling convention must stay in [0, 1]."""
    cohort, meta = _cohort()
    ph = condition.binary(meta, "cmv")
    r = association(cohort, ph, test="permutation", alternative="two-sided", n_perm=500, seed=0)
    p = r["p_value"].drop_nulls().to_numpy()
    assert ((p >= 0) & (p <= 1)).all()
