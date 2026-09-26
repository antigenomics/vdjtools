"""The D-segment walk shares its prefixes across 3' trims — and nothing it returns moved.

``pgen_aa_vdj`` used to re-thread the surviving D germline from scratch for every ``(idx5, idx3)``
trim, though the emission at ``idx3`` is a *prefix* of the one at ``idx3 - 1``. One forward pass per
``(D, idx5)`` now covers every 3' trim of it. The identity is exact, so this file pins two things
a tolerance-based oracle cannot:

* **every 3' trim still contributes** — a prefix loop that stops one nucleotide short drops the
  longest D of each ``(D, idx5)`` chain, which on a real model is a few percent of the mass and
  looks like nothing at ``rtol=1e-6``;
* **the values are bit-for-bit what 3.16.0 returned** — ``vsig:pgen:*:frac_atypical`` is compared
  against a frozen ``pgen_q05``, so a last-bit move is a silent data bug, not a rounding detail.

Why a tolerance is not enough here, measured 2026-09-26: stopping the shared walk one nucleotide
short (so each ``(D, 5' cut)`` chain loses its longest D) moves IGH Pgen by a **median relative
error of 7.9e-7** over 12 junctions off the bundled model -- under the ``rtol=1e-6`` every
OLGA-comparison test in this repo uses, even though the worst junction moves 4.3%. The frozen
reference below catches that mutation on every value; a tolerance catches it only if the draw
happens to include a long junction.

IGH is where the cost was (9,212 live ``(D, ndel5, ndel3)`` states against TRB's 297) and it had no
OLGA comparison anywhere in the suite until this file.
"""
from __future__ import annotations

import itertools
import json
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from vdjtools.model import Event, EventKind, Manifest, Model, from_olga, load_bundled, native
from vdjtools.model.reference import _CODON_TABLE as CODON_TABLE
from vdjtools.model.generate import generate

OLGA_MODELS = Path(
    os.environ.get(
        "VDJTOOLS_OLGA_MODELS",
        str(Path(__file__).resolve().parent / "fixtures" / "olga" / "default_models"),
    )
)
GOLDEN = Path(__file__).resolve().parent / "fixtures" / "pgen_golden.json"


# ---- every 3' trim is covered, on a model small enough to enumerate by hand ------------------

def _tiny_vdj(d_del_rows) -> Model:
    """A minimal single-D VDJ model: V "TGT", D "CAG", J "TTT", insertions of length 0 or 1.

    ``d_del_rows`` is ``[(ndel5, ndel3, p), ...]``. D is 3 nt, so a 9-nt (3-aa) CDR3 forces
    ``insVD + surviving D + insDJ == 3`` and several 3' trims of the same 5' cut are reachable at
    once — which is exactly the case the shared prefix walk has to cover.
    """
    def genes(seg, allele, cut):
        return pl.DataFrame({f"{seg}_allele": [allele], "cut_segment": [cut], "functional": [True]})

    ins = pl.DataFrame({"length": pl.Series([0, 1, 2], dtype=pl.Int16), "p": [0.5, 0.3, 0.2]})
    dinucl = pl.DataFrame({
        "from_nt": pl.Series([f for f in range(4) for _ in range(4)], dtype=pl.UInt8),
        "to_nt": pl.Series([t for _ in range(4) for t in range(4)], dtype=pl.UInt8),
        "p": [0.25] * 16,
    })
    ddel = pl.DataFrame({
        "d_allele": ["D"] * len(d_del_rows),
        "ndel5": pl.Series([r[0] for r in d_del_rows], dtype=pl.Int16),
        "ndel3": pl.Series([r[1] for r in d_del_rows], dtype=pl.Int16),
        "p": [r[2] for r in d_del_rows],
    })
    tables = {
        "v_choice": pl.DataFrame({"v_allele": ["V"], "p": [1.0]}),
        "j_choice": pl.DataFrame({"j_allele": ["J"], "p": [1.0]}),
        "d_gene": pl.DataFrame({"j_allele": ["J"], "d_allele": ["D"], "p": [1.0]}),
        "n_d": pl.DataFrame({"n_d": pl.Series([1], dtype=pl.UInt8), "p": [1.0]}),
        "v_3_del": pl.DataFrame({"v_allele": ["V"], "ndel": pl.Series([0], dtype=pl.Int16), "p": [1.0]}),
        "j_5_del": pl.DataFrame({"j_allele": ["J"], "ndel": pl.Series([0], dtype=pl.Int16), "p": [1.0]}),
        "d_del": ddel,
        "vd_ins": ins, "dj_ins": ins,
        "vd_dinucl": dinucl, "dj_dinucl": dinucl,
    }
    events = {
        "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
        "j_choice": Event("j_choice", EventKind.GENE_CHOICE),
        "d_gene": Event("d_gene", EventKind.GENE_CHOICE, ("j_choice",)),
        "n_d": Event("n_d", EventKind.N_D),
        "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
        "j_5_del": Event("j_5_del", EventKind.DELETION, ("j_choice",)),
        "d_del": Event("d_del", EventKind.DELETION_2D, ("d_gene",)),
        "vd_ins": Event("vd_ins", EventKind.INS_LENGTH),
        "dj_ins": Event("dj_ins", EventKind.INS_LENGTH),
        "vd_dinucl": Event("vd_dinucl", EventKind.DINUCLEOTIDE),
        "dj_dinucl": Event("dj_dinucl", EventKind.DINUCLEOTIDE),
    }
    manifest = Manifest(locus="TEST", organism="synthetic", chain_type="VDJ", events=events,
                        palindrome_max={"v_3": 0, "j_5": 0, "d_5": 0, "d_3": 0}, source="tiny")
    genomic = {"genes_v": genes("v", "V", "TGT"), "genes_j": genes("j", "J", "TTT"),
               "genes_d": genes("d", "D", "CAG")}
    return Model(manifest=manifest, tables=tables, genomic=genomic)


#: Five live D trims. ``(0, 0)``, ``(0, 1)`` and ``(0, 2)`` share one 5' cut and differ only in how
#: far the shared walk runs, so a prefix loop that stops early loses whichever it stops before.
TRIMS = [(0, 0, 0.3), (0, 1, 0.25), (0, 2, 0.2), (1, 0, 0.15), (1, 1, 0.1)]


#: ``CQF`` reaches all five: the D's own "CAG" untrimmed, then its 2-nt and 1-nt prefixes with the
#: insertions making up the difference, then the same from the 5'-trimmed "AG". That is one whole
#: prefix chain plus a second one, which is the structure the shared walk has to reproduce.
@pytest.mark.parametrize("aa,n_live", [("CQF", 5), ("CSF", 4), ("CPF", 3)])
def test_every_3p_trim_contributes_its_own_share(aa, n_live):
    """Pgen over five trims == the weighted sum of the five one-trim Pgens.

    Pgen is linear in ``P(delD)``, so this decomposition is exact — and it fails loudly and
    specifically if the shared walk covers only some prefixes of a ``(D, 5' cut)`` chain. An
    oracle comparison would only say "a bit low"; this says which trim went missing.
    """
    full = native.pgen_aa(_tiny_vdj(TRIMS), aa)
    parts = {(n5, n3): p * native.pgen_aa(_tiny_vdj([(n5, n3, 1.0)]), aa) for n5, n3, p in TRIMS}
    assert sum(v > 0 for v in parts.values()) == n_live, f"model changed shape: {parts}"
    assert full == pytest.approx(sum(parts.values()), rel=1e-12), parts


def test_the_shared_walk_agrees_with_the_untouched_nt_path():
    """aa Pgen == sum of nt Pgen over synonymous codons, on a model with several live trims.

    ``pgen_nt`` reaches the D through ``d_middle``, a separate enumeration this change never
    touched, so it is an independent check on the aa transfer matrix rather than a restatement
    of it.
    """
    m = _tiny_vdj(TRIMS)
    syn = defaultdict(list)
    for cod, a in CODON_TABLE.items():
        syn[a].append(cod)
    for aa in ("CQF", "CSF", "CPF"):
        brute = sum(native.pgen_nt(m, "".join(c)) for c in itertools.product(*[syn[a] for a in aa]))
        assert brute > 0
        assert native.pgen_aa(m, aa) == pytest.approx(brute, rel=1e-12)


# ---- IGH against OLGA: the locus this change is for, and it had no oracle test ---------------

def _olga_vdj(sub):
    from olga import generation_probability as gp
    from olga import load_model as ol

    d = OLGA_MODELS / sub
    gen = ol.GenerativeModelVDJ()
    gen.load_and_process_igor_model(str(d / "model_marginals.txt"))
    gd = ol.GenomicDataVDJ()
    gd.load_igor_genomic_data(str(d / "model_params.txt"), str(d / "V_gene_CDR3_anchors.csv"),
                              str(d / "J_gene_CDR3_anchors.csv"))
    return gp.GenerationProbabilityVDJ(gen, gd)


@pytest.fixture(scope="module")
def igh():
    """IGH, its OLGA oracle, and a length-spread draw.

    Spread, not the three shortest: a D walk that drops mass is a *relative* error that grows with
    how much D there is to drop, so a short-junction draw is the least sensitive one available.
    """
    pytest.importorskip("olga.load_model", reason="olga (the [oracle] extra) not installed")
    if not OLGA_MODELS.exists():
        pytest.skip(f"OLGA models not at {OLGA_MODELS}")
    m = from_olga(OLGA_MODELS / "human_B_heavy", locus="IGH")
    rows = sorted(generate(m, 60, seed=7, productive_only=True).to_dicts(),
                  key=lambda r: len(r["junction_aa"]))
    return m, _olga_vdj("human_B_heavy"), [rows[0], rows[len(rows) // 2], rows[-1]]


def test_igh_nt_pgen_matches_olga(igh):
    """IGH is the 9,212-D-state locus the walk was rewritten for; nothing compared it to OLGA."""
    m, olga, rows = igh
    checked = 0
    for r in rows:
        nt, v, j = r["junction_nt"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_nt(m, nt, v, j), olga.compute_nt_CDR3_pgen(nt, v, j), rtol=1e-9)
        checked += 1
    assert checked == 3


def test_igh_aa_pgen_matches_olga(igh):
    """The aa transfer matrix on IGH, V/J-restricted and marginalised, against OLGA."""
    m, olga, rows = igh
    checked = 0
    for r in rows:
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        assert np.isclose(native.pgen_aa(m, aa, v, j), olga.compute_aa_CDR3_pgen(aa, v, j), rtol=1e-9)
        checked += 1
    assert checked == 3


@pytest.mark.slow
def test_igh_aa_pgen_matches_olga_marginalised(igh):
    """Same, summing over every V/J/D — the form ``pgen_block`` actually calls. OLGA's own aa path
    is the slow one here (it is why the concordance harness draws 12 aa against 150 nt)."""
    m, olga, rows = igh
    for r in rows[:2]:
        aa = r["junction_aa"]
        assert np.isclose(native.pgen_aa(m, aa), olga.compute_aa_CDR3_pgen(aa), rtol=1e-9)


# ---- the frozen reference -------------------------------------------------------------------

def test_pgen_matches_the_frozen_reference():
    """Bit-for-bit equality with the values 3.16.0 returned, on the three D-bearing loci.

    Hex float64, so the fixture round-trips exactly and the comparison is ``==`` rather than a
    tolerance. This is the only stored Pgen reference in the repo: the signature block's
    ``frac_atypical`` is measured against a ``pgen_q05`` frozen outside it, so a shift of 1e-6
    lands in a plausible range and no other test would see it.

    Regenerate ONLY when a germline or a bundled model changes — never to make a red test green.
    """
    want = json.loads(GOLDEN.read_text())
    for locus, block in sorted(want["values"].items()):
        m = load_bundled(locus, "olga")
        aas = [s[1] for s in block["seqs"]]
        vs = [s[2] for s in block["seqs"]]
        js = [s[3] for s in block["seqs"]]
        got = {
            "aa_free": native.pgen_aa_batch(m, aas, threads=1),
            "aa_vj": native.pgen_aa_batch(m, aas, vs, js, threads=1),
            "mm1_free": native.pgen_aa_batch(m, aas[:4], mismatches=1, threads=1),
            "nt_vj": [native.pgen_nt(m, s[0], s[2], s[3]) for s in block["seqs"][:8]],
        }
        for key, vals in sorted(got.items()):
            exp = [float.fromhex(h) for h in block[key]]
            assert len(vals) == len(exp)
            bad = [(i, e, g) for i, (e, g) in enumerate(zip(exp, vals)) if e != g]
            assert not bad, f"{locus}/{key}: {len(bad)} of {len(exp)} moved, first {bad[0]}"
