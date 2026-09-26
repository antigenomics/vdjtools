"""Every knob that moves a number must move only the columns it is documented to move.

This file exists for a class of bug rather than for one bug. `DEFAULT_CSTAR = 0.20` produced a
fully populated, entirely plausible `vsig:div` block at a coverage level nobody established, and
nothing in the emitted vector said so -- so a matrix built on the fallback and one built on
measured constants looked identical and could be joined. The fix was a reported fraction; the
guard is this: turn each knob, take the set of columns that moved, and check it against the set
the knob is allowed to touch.

Seven loci on purpose. A TRB-only fixture cannot see that a reference covers TRA and TRB and
leaves five loci of diversity as holes, which is the case that used to raise `KeyError`.
"""
from __future__ import annotations

import warnings

import numpy as np
import polars as pl
import pytest

from vdjtools.signature import assemble, layout as L

_AA = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
_LOCI = {"TRB": ("TRBV{}*01", "TRBJ{}*01", 28, 13), "IGH": ("IGHV{}*01", "IGHJ{}*01", 7, 6),
         "IGK": ("IGKV{}*01", "IGKJ{}*01", 6, 5), "TRA": ("TRAV{}*01", "TRAJ{}*01", 40, 55),
         "IGL": ("IGLV{}*01", "IGLJ{}*01", 5, 5), "TRG": ("TRGV{}*01", "TRGJ{}*01", 6, 3),
         "TRD": ("TRDV{}*01", "TRDJ{}*01", 3, 4)}


@pytest.fixture(scope="module")
def sample():
    r = np.random.default_rng(7)
    out = {}
    for loc, (vf, jf, nv, nj) in _LOCI.items():
        n = int(r.integers(120, 400))
        ln = r.integers(11, 18, n)
        out[loc] = pl.DataFrame({
            "junction_aa": ["C" + "".join(_AA[r.integers(0, 20, k)]) + "F" for k in ln],
            "v_call": [vf.format(v) for v in r.integers(1, nv, n)],
            "j_call": [jf.format(v) for v in r.integers(1, nj, n)],
            "duplicate_count": np.ceil(r.zipf(1.4, n).clip(1, 2000)).astype(int).tolist()})
    return out


def _vsig(sample, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return assemble.vsig(sample, **kw)


def _moved(a: dict, b: dict) -> set[str]:
    return {k for k in a
            if not (a[k] == b[k] or (np.isnan(a[k]) and np.isnan(b[k])))}


def test_cstar_moves_only_diversity_and_its_own_provenance(sample):
    """The knob that used to be invisible. Its blast radius is three column families.

    `vsig:div:*` is the quantity, `mask:*:estimable` says whether it could be produced at all,
    and `qc:-:cstar_fallback_frac` says where the level came from. Anything else moving means a
    coverage level has leaked into a block that does not standardise to one.
    """
    flat = _vsig(sample, tier="standard")
    measured = _vsig(sample, tier="standard",
                     cstar={loc: 0.11 for loc in _LOCI})
    moved = _moved(flat, measured)
    assert moved, "0.20 and 0.11 must give different Hill numbers on this sample"
    allowed = {c for c in moved
               if c.startswith("vsig:div:") or c.endswith(":estimable")
               or c == "vsig:qc:-:cstar_fallback_frac"}
    assert moved == allowed, f"cstar leaked into {sorted(moved - allowed)}"
    assert flat["vsig:qc:-:cstar_fallback_frac"] == 1.0
    assert measured["vsig:qc:-:cstar_fallback_frac"] == 0.0


def test_a_reference_covering_two_loci_holes_the_other_five_and_says_so(sample):
    """What the shipped amplicon reference actually is. This used to be a `KeyError`."""
    with pytest.warns(UserWarning, match="no coverage level was supplied"):
        out = assemble.vsig(sample, tier="standard", cstar={"TRA": 0.545, "TRB": 0.408})
    for loc in ("TRA", "TRB"):
        assert np.isfinite(out[f"vsig:div:{loc}:1D_c"]) and out[f"vsig:mask:{loc}:estimable"] == 1.0
    for loc in ("IGH", "IGK", "IGL", "TRG", "TRD"):
        assert np.isnan(out[f"vsig:div:{loc}:1D_c"]), loc
        assert out[f"vsig:mask:{loc}:estimable"] == 0.0, loc
        assert np.isfinite(out[f"vsig:depth:{loc}:reads"]), f"{loc} depth needs no coverage level"
    # a level was established for every locus that got one, so nothing was borrowed
    assert out["vsig:qc:-:cstar_fallback_frac"] == 0.0


def test_cstar_none_and_a_covering_dict_agree_everywhere_outside_diversity(sample):
    """Refusing the level must not perturb a single count, length or pairing statistic."""
    holed = _vsig(sample, tier="standard", cstar=None)
    measured = _vsig(sample, tier="standard", cstar={loc: 0.11 for loc in _LOCI})
    moved = _moved(holed, measured)
    assert moved, "declining the level must change the diversity block"
    assert all(c.startswith("vsig:div:") or c.endswith(":estimable") for c in moved), \
        sorted(c for c in moved if not c.startswith("vsig:div:"))


def test_pgen_n_max_moves_only_the_pgen_block(sample):
    """It is a documented accuracy/cost knob, and the accuracy it costs is confined."""
    few = _vsig(sample, tier="standard", pgen_n_max=25)
    more = _vsig(sample, tier="standard", pgen_n_max=120)
    moved = _moved(few, more)
    assert moved, "a different subsample must give different Pgen moments"
    assert all(c.startswith("vsig:pgen:") for c in moved), sorted(moved)[:5]


def test_the_kernel_thread_count_moves_nothing(sample):
    """`threads` buys cores, never a different answer -- on all seven loci."""
    assert not _moved(_vsig(sample, tier="standard", threads=1),
                      _vsig(sample, tier="standard", threads=0))


def test_declining_a_block_moves_nothing_else_on_seven_loci(sample):
    """The 3.16.0 claim, re-checked at the cohort's real locus breadth."""
    keep = [c for c in L.columns("standard", "vsig") if not c.startswith("vsig:pgen:")]
    full = _vsig(sample, tier="standard")
    cut = _vsig(sample, tier="standard", columns=keep)
    assert set(cut) == set(keep)
    assert not _moved(cut, {k: full[k] for k in cut})
