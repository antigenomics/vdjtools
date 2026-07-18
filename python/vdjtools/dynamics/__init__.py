"""Longitudinal clonotype dynamics — did a clonotype change *within* a donor?

Where :mod:`vdjtools.biomarker` tests **incidence across subjects** (is this clonotype
carried by more cases than controls?), this subpackage tests **frequency across timepoints**
within one subject (did this clonotype's frequency change between pre and post?). The two ask
different questions of different designs, and a paired design answered with an incidence test
is a known way to manufacture hits.

Three complementary tools:

* :func:`test_pair` — the per-clonotype paired test. The engine is the Ayestaran (Cambridge,
  2024) effective-sample-size method: sequencing is a *two-step* sampling process (molecules,
  then reads), so the sample size that drives the noise is the harmonic sum
  ``1/N_eff = 1/N_S1 + 1/N_seq`` — dominated by whichever step is smaller. ``N_eff`` is estimated
  **per pair** from that pair's own mean–variance scaling and is not computable from a single
  sample, so there is deliberately no way to spell "normalise the whole cohort to a common
  N_eff" with this API.
* :func:`test_metaclonotypes` — the same test on **1-Hamming / 1-Levenshtein CDR3 groups**, so a
  convergent (near-variant) expansion is tested as one feature.
* :func:`capture_rates` / :func:`capture_test` — the VDJtrack **recapture model** (Pavlova,
  Zvyagin & Shugay 2024): bin clonotypes by size (singleton/doubleton/tripleton/large), measure
  the Poisson-driven recapture rate per group, and test whether an annotated group persists more
  than that size baseline.
"""
from __future__ import annotations

from .capture import (
    SIZE_CLASSES,
    capture_paired_test,
    capture_rates,
    capture_test,
    poisson_capture,
    size_class,
)
from .expansion import expansion_test
from .groups import test_metaclonotypes
from .paired import DEFAULT_KEY, estimate_neff, test_pair

__all__ = [
    "DEFAULT_KEY",
    "estimate_neff",
    "test_pair",
    "test_metaclonotypes",
    "expansion_test",
    "SIZE_CLASSES",
    "size_class",
    "poisson_capture",
    "capture_rates",
    "capture_test",
    "capture_paired_test",
]
