"""Explicit post-selection mass scaling; this is not a fitted selection model."""
import numpy as np

# ALICE's published average on the OLGA TRB scale; not established for TRA.
ALICE_TRB_SELECTION_Q = 9.41


def selection_mass(pgen, selection_q=1.0):
    """Scale generation mass by an explicit positive selection factor.

    Returns Q * Pgen without clipping or renormalization. The result is the
    mass used in a generative expectation, not a calibrated binding probability.
    Empirical control frequencies already reflect selection; do not rescale them.
    """
    mass = np.asarray(pgen, dtype=float)
    q = float(selection_q)
    if not np.isfinite(q) or q <= 0 or not np.all(np.isfinite(mass) & (mass >= 0)):
        raise ValueError("finite nonnegative generation mass and positive selection_q required")
    with np.errstate(over="ignore"):
        scaled = mass * q
    if not np.all(np.isfinite(scaled)):
        raise ValueError("selection-scaled mass overflow")
    return scaled.item() if scaled.ndim == 0 else scaled
