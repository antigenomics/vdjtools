"""Phase 1f — native C++ Pgen matches the Python reference (and OLGA) exactly.

The native path (``vdjtools.model.native``) reads the same PackedModel and must agree with the
pure-Python ``pgen`` to numerical tolerance — it is only faster (≈90x for VDJ).
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from vdjtools.model import from_olga, native
from vdjtools.model.generate import generate
from vdjtools.model.pgen import pgen_aa as py_pgen_aa
from vdjtools.model.pgen import pgen_nt as py_pgen_nt
from vdjtools.model.pgen import prepare

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        "/Users/mikesh/vcs/code/mirpy/mir/resources/olga/default_models",
    )
)
pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
pytestmark = pytest.mark.skipif(not OLGA_MODELS.exists(), reason=f"OLGA models not at {OLGA_MODELS}")


def test_native_nt_beta_oracle():
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    p = native.pgen_nt(m, "TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC", "TRBV30*01", "TRBJ1-2*01")
    assert np.isclose(p, 2.3986503758867323e-12, rtol=1e-9)


@pytest.mark.parametrize("sub,locus", [("human_T_alpha", "TRA"), ("human_T_beta", "TRB")])
def test_native_matches_python(sub, locus):
    m = from_olga(OLGA_MODELS / sub, locus=locus)
    prep = prepare(m)
    df = generate(m, 8, seed=3, productive_only=True)
    for r in df.to_dicts():
        nt, v, j = r["cdr3_nt"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_nt(m, nt, v, j), py_pgen_nt(prep, nt, v, j), rtol=1e-9)
    # unrestricted (sum over all genes) must also agree
    nt = df["cdr3_nt"][0]
    assert np.isclose(native.pgen_nt(m, nt), py_pgen_nt(prep, nt), rtol=1e-9)


def test_native_aa_beta_oracle():
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    p = native.pgen_aa(m, "CAWSVAPDRGGYTF", "TRBV30*01", "TRBJ1-2*01")
    assert np.isclose(p, 1.203646865765782e-10, rtol=1e-6)


def test_native_aa_matches_python_vj():
    m = from_olga(OLGA_MODELS / "human_T_alpha", locus="TRA")
    prep = prepare(m)
    for r in generate(m, 6, seed=3, productive_only=True).to_dicts():
        aa, v, j = r["cdr3_aa"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_aa(m, aa, v, j), py_pgen_aa(prep, aa, v, j), rtol=1e-6)


@pytest.mark.slow
def test_native_aa_matches_python_vdj():
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    prep = prepare(m)
    # short beta CDR3s (the VDJ aa enumeration is O(D x deletions x positions))
    for r in sorted(generate(m, 60, seed=3, productive_only=True).to_dicts(), key=lambda r: len(r["cdr3_aa"]))[:3]:
        aa, v, j = r["cdr3_aa"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_aa(m, aa, v, j), py_pgen_aa(prep, aa, v, j), rtol=1e-6)
