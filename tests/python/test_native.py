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


def _olga_vdj(sub):
    """Build OLGA's own GenerationProbabilityVDJ for the given model subdir."""
    from olga import load_model as ol
    from olga import generation_probability as gp

    d = OLGA_MODELS / sub
    gen = ol.GenerativeModelVDJ()
    gen.load_and_process_igor_model(str(d / "model_marginals.txt"))
    gd = ol.GenomicDataVDJ()
    gd.load_igor_genomic_data(
        str(d / "model_params.txt"),
        str(d / "V_gene_CDR3_anchors.csv"),
        str(d / "J_gene_CDR3_anchors.csv"),
    )
    return gp.GenerationProbabilityVDJ(gen, gd)


def test_native_aa_vdj_matches_olga():
    """Native transfer-matrix VDJ aa-Pgen == OLGA's compute_aa_CDR3_pgen (both fast, exact).

    This is the primary VDJ aa-Pgen check: it compares the native ``_core`` transfer matrix
    directly against OLGA on unrestricted (sum over all V/J/D) and V/J-restricted queries. The
    pure-Python enumeration is too slow to be a routine oracle here (``test_native_aa_matches_
    python_vdj`` covers it under ``-m slow``); OLGA is the authoritative fast oracle.
    """
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    olga = _olga_vdj("human_T_beta")
    seqs = sorted(
        generate(m, 40, seed=7, productive_only=True).to_dicts(), key=lambda r: len(r["cdr3_aa"])
    )
    for r in seqs[:12]:
        aa, v, j = r["cdr3_aa"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_aa(m, aa), olga.compute_aa_CDR3_pgen(aa), rtol=1e-9)  # all V/J/D
        assert np.isclose(native.pgen_aa(m, aa, v, j), olga.compute_aa_CDR3_pgen(aa, v, j), rtol=1e-9)


def test_native_aa_hamming1_matches_olga():
    """Native 1-mismatch aa-Pgen (``mismatches=1``) == OLGA's compute_hamming_dist_1_pgen.

    The native path evaluates the inclusion-exclusion identity with the fast transfer matrix and
    a single wildcard mask per position (no 19x per-neighbour enumeration), so it is several times
    faster than OLGA while matching it to machine precision — unrestricted and V/J-restricted.
    """
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    olga = _olga_vdj("human_T_beta")
    seqs = sorted(
        generate(m, 40, seed=7, productive_only=True).to_dicts(), key=lambda r: len(r["cdr3_aa"])
    )
    for r in seqs[:8]:
        aa, v, j = r["cdr3_aa"], r["v_call"], r["j_call"]
        assert np.isclose(
            native.pgen_aa(m, aa, mismatches=1), olga.compute_hamming_dist_1_pgen(aa), rtol=1e-9
        )
        assert np.isclose(
            native.pgen_aa(m, aa, v, j, mismatches=1),
            olga.compute_hamming_dist_1_pgen(aa, v, j),
            rtol=1e-9,
        )


def test_native_aa_matches_python_vj():
    m = from_olga(OLGA_MODELS / "human_T_alpha", locus="TRA")
    prep = prepare(m)
    for r in generate(m, 6, seed=3, productive_only=True).to_dicts():
        aa, v, j = r["cdr3_aa"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_aa(m, aa, v, j), py_pgen_aa(prep, aa, v, j), rtol=1e-6)


def test_native_estep_matches_python():
    """The native EM E-step's soft counts match the pure-Python E-step to machine epsilon."""
    from collections import defaultdict

    from vdjtools._core import estep_batch, make_counts
    from vdjtools.model.infer import _estep_seq
    from vdjtools.model.native import _encode, pack

    m = from_olga(OLGA_MODELS / "human_T_alpha", locus="TRA")
    prep = prepare(m)
    seqs = generate(m, 8, seed=3, productive_only=True)["cdr3_nt"].to_list()
    pm, vi, ji = pack(m)
    counts = make_counts(pm)
    estep_batch(pm, [_encode(s) for s in seqs], [], [], [], counts)

    pyc = {k: defaultdict(float) for k in m.tables if k != "n_d"}
    for s in seqs:
        _estep_seq(prep, s.upper(), pyc)
    nJ = len(ji)
    nat = np.array(counts.v_choice)
    py = np.zeros(len(nat))
    for (v,), w in pyc["v_choice"].items():
        py[vi[v]] = w
    assert np.allclose(nat, py, atol=1e-12)
    assert np.isclose(nat.sum(), len(seqs))  # each read's posteriors sum to 1

    nat_dn = np.array(counts.dinucl_vj)
    py_dn = np.zeros(16)
    for (fr, to), w in pyc["vj_dinucl"].items():
        py_dn[to * 4 + fr] = w
    assert np.allclose(nat_dn, py_dn, atol=1e-12)


@pytest.mark.slow
def test_native_em_recovers_tra():
    """Native EM recovers the source model's marginals (same as pure-Python EM, ~100x faster)."""
    import collections

    from vdjtools.model.infer import infer_native

    m = from_olga(OLGA_MODELS / "human_T_alpha", locus="TRA")
    seqs = generate(m, 1500, seed=11)["cdr3_nt"].to_list()
    fit, rep = infer_native(m, seqs, max_iter=10)

    def corr(a, b):
        k = set(a) | set(b)
        return float(np.corrcoef([a.get(x, 0) for x in k], [b.get(x, 0) for x in k])[0, 1])

    def td(model, ev, keys):
        return {tuple(r[:-1]): r[-1] for r in model.tables[ev].select([*keys, "p"]).iter_rows()}

    grp = {a: s for a, s in zip(m.genomic["genes_v"]["v_allele"], m.genomic["genes_v"]["cut_segment"])}

    def gv(model):
        d = collections.defaultdict(float)
        for v, p in zip(model.tables["v_choice"]["v_allele"], model.tables["v_choice"]["p"]):
            d[grp.get(v, v)] += p
        return d

    assert corr(gv(fit), gv(m)) > 0.9
    assert corr(td(fit, "vj_ins", ["length"]), td(m, "vj_ins", ["length"])) > 0.95
    assert corr(td(fit, "vj_dinucl", ["from_nt", "to_nt"]), td(m, "vj_dinucl", ["from_nt", "to_nt"])) > 0.98


@pytest.mark.slow
def test_native_aa_matches_python_vdj():
    m = from_olga(OLGA_MODELS / "human_T_beta", locus="TRB")
    prep = prepare(m)
    # short beta CDR3s (the VDJ aa enumeration is O(D x deletions x positions))
    for r in sorted(generate(m, 60, seed=3, productive_only=True).to_dicts(), key=lambda r: len(r["cdr3_aa"]))[:3]:
        aa, v, j = r["cdr3_aa"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_aa(m, aa, v, j), py_pgen_aa(prep, aa, v, j), rtol=1e-6)
