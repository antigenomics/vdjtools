"""Degenerate (motif) aa-Pgen: the masked DP exposed as ``native.pgen_aa_degenerate``.

The masked transfer-matrix DP has always been the engine behind ``pgen_aa`` (every position pinned
to one residue) and ``pgen_aa_hamming1`` (one wildcard position at a time), but only those two
special cases were reachable from Python. These tests pin the general per-position residue-set form
against both, and against the two identities that make it usable as a motif Pgen: the sum over the
sequences a motif matches, and the all-wildcard length marginal.
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest

from vdjtools.model import load_bundled, native
from vdjtools.model.generate import generate

AA20 = "ACDEFGHIKLMNPQRSTVWY"


@pytest.fixture(scope="module", params=["TRB", "TRG"])
def locus_model(request):
    """A bundled model (VDJ TRB and VJ TRG — the masked DP has a separate branch each) + reads."""
    m = load_bundled(request.param, "olga")
    reads = generate(m, 6, seed=5, productive_only=True).to_dicts()
    return m, reads


# --- pinned == pgen_aa -------------------------------------------------------------------------

def test_pinned_positions_equal_pgen_aa(locus_model):
    """One residue per position is exactly the ``pgen_aa`` query — bitwise, not just close."""
    m, reads = locus_model
    for r in reads:
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        assert native.pgen_aa_degenerate(m, list(aa)) == native.pgen_aa(m, aa)
        assert native.pgen_aa_degenerate(m, list(aa), v, j) == native.pgen_aa(m, aa, v, j)


def test_residue_sets_sum_the_matching_sequences(locus_model):
    """A motif's Pgen is the total Pgen of every sequence it matches (the defining identity)."""
    m, reads = locus_model
    r = reads[0]
    aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
    allowed = list(aa)
    for k in (2, 4, 6):  # widen three positions to two residues each -> 8 matching sequences
        allowed[k] = aa[k] + ("A" if aa[k] != "A" else "G")
    members = ["".join(c) for c in itertools.product(*allowed)]
    assert len(members) == 8
    got = native.pgen_aa_degenerate(m, allowed, v, j)
    want = sum(native.pgen_aa(m, s, v, j) for s in members)
    assert got > 0.0
    assert np.isclose(got, want, rtol=1e-12)


def test_all_wildcard_is_the_length_marginal(toy_model):
    """Every position wildcarded == P(L, V, J): the summed Pgen of all 20**L aa sequences.

    Checked by brute force on the toy VJ locus, where L=3 (8000 sequences) is enumerable.
    """
    p = native.pgen_aa_degenerate(toy_model, ["X"] * 3)
    every = native.pgen_aa_batch(toy_model, ["".join(c) for c in itertools.product(AA20, repeat=3)])
    assert p > 0.0
    assert np.isclose(p, sum(every), rtol=1e-12)


# --- monotonicity, wildcards, the Hamming-1 identity --------------------------------------------

def test_widening_a_position_never_decreases_pgen(locus_model):
    """Set inclusion is monotone: a larger residue set at any position can only add mass."""
    m, reads = locus_model
    r = reads[0]
    aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
    for k in range(len(aa)):
        pinned = list(aa)
        pair = list(aa)
        pair[k] = aa[k] + ("W" if aa[k] != "W" else "G")
        wild = list(aa)
        wild[k] = "X"
        p1 = native.pgen_aa_degenerate(m, pinned, v, j)
        p2 = native.pgen_aa_degenerate(m, pair, v, j)
        p3 = native.pgen_aa_degenerate(m, wild, v, j)
        assert p2 >= p1 * (1 - 1e-12)
        assert p3 >= p2 * (1 - 1e-12)


def test_empty_and_X_mean_wildcard(locus_model):
    """``""`` and ``"X"`` are both wildcards — and are NOT what ``pgen_aa`` does with an X.

    ``mask_for_aa`` matches by exact character against the genetic code, so an ``X`` in a plain
    ``pgen_aa`` query is an empty codon set and scores a silent 0.0. The degenerate wrapper maps it
    to the full 20-residue mask instead.
    """
    m, reads = locus_model
    aa = reads[0]["junction_aa"]
    k = len(aa) // 2
    empty, ex, spelled = list(aa), list(aa), list(aa)
    empty[k], ex[k], spelled[k] = "", "X", AA20
    p_empty = native.pgen_aa_degenerate(m, empty)
    assert p_empty > 0.0
    assert p_empty == native.pgen_aa_degenerate(m, ex)
    assert p_empty == native.pgen_aa_degenerate(m, spelled)
    assert native.pgen_aa(m, aa[:k] + "X" + aa[k + 1:]) == 0.0  # the trap this closes


def test_hamming1_identity_from_the_public_api(locus_model):
    """sum_k Pgen(position k wildcarded) - (L-1) Pgen(centre) reproduces ``mismatches=1``."""
    m, reads = locus_model
    for r in reads[:3]:
        aa, v, j = r["junction_aa"], r["v_call"], r["j_call"]
        centre = native.pgen_aa(m, aa, v, j)
        ball = 0.0
        for k in range(len(aa)):
            wild = list(aa)
            wild[k] = "X"
            ball += native.pgen_aa_degenerate(m, wild, v, j)
        ball -= (len(aa) - 1) * centre
        assert np.isclose(ball, native.pgen_aa(m, aa, v, j, mismatches=1), rtol=1e-12)


# --- V/J conditioning and argument validation ---------------------------------------------------

def test_v_j_conditioning_matches_pgen_aa():
    """v/j behave exactly as in ``pgen_aa``: joint when given, marginal when ``None``."""
    m = load_bundled("TRB", "olga")
    allowed = list("CASSVGLYSTDTQYF")
    allowed[4] = "VIL"
    agnostic = native.pgen_aa_degenerate(m, allowed)
    joint = native.pgen_aa_degenerate(m, allowed, "TRBV9*01", "TRBJ2-3*01")
    assert 0.0 < joint < agnostic
    # the marginal is the sum of the joints over all V, at a pinned J
    at_j = native.pgen_aa_degenerate(m, allowed, None, "TRBJ2-3*01")
    v_alleles = m.genomic["genes_v"]["v_allele"].to_list()
    total = sum(native.pgen_aa_degenerate(m, allowed, v, "TRBJ2-3*01") for v in v_alleles)
    assert np.isclose(at_j, total, rtol=1e-9)


def test_gene_level_name_raises_not_silently_marginalizes():
    """A gene-level V/J must raise here exactly as it does for ``pgen_aa`` (the -1 trap)."""
    m = load_bundled("TRB", "olga")
    allowed = list("CASSVGLYSTDTQYF")
    with pytest.raises(KeyError, match="gene name"):
        native.pgen_aa_degenerate(m, allowed, "TRBV9", "TRBJ2-3*01")
    with pytest.raises(KeyError, match="gene name"):
        native.pgen_aa_degenerate_batch(m, [allowed], ["TRBV9"], ["TRBJ2-3*01"])
    with pytest.raises(KeyError, match="not in the model"):
        native.pgen_aa_degenerate(m, allowed, "TRBV999*01")


@pytest.mark.parametrize("bad", ["B", "AZ", "1"])
def test_unknown_residue_raises(bad):
    """An unrecognised residue has no codons — it must raise, not score a silent 0.0."""
    m = load_bundled("TRG", "olga")
    allowed = list("CATWDKQLGKKIKVF")
    allowed[3] = bad
    with pytest.raises(ValueError, match="amino acid"):
        native.pgen_aa_degenerate(m, allowed)


# --- batch --------------------------------------------------------------------------------------

def test_batch_equals_serial_and_is_thread_invariant():
    """The batch form is bitwise-identical to per-sequence, and to ``pgen_aa_batch`` when pinned."""
    m = load_bundled("TRB", "olga")
    g = generate(m, 128, seed=5, productive_only=True)  # >64 so the threaded path actually threads
    seqs = g["junction_aa"].to_list()
    vs, js = g["v_call"].to_list(), g["j_call"].to_list()
    allowed = [list(s) for s in seqs]
    serial = [native.pgen_aa_degenerate(m, a, v, j) for a, v, j in zip(allowed, vs, js)]
    assert native.pgen_aa_degenerate_batch(m, allowed, vs, js, threads=1) == serial
    assert native.pgen_aa_degenerate_batch(m, allowed, vs, js, threads=8) == serial
    assert native.pgen_aa_degenerate_batch(m, allowed, vs, js, threads=0) == serial
    assert serial == native.pgen_aa_batch(m, seqs, vs, js)
    assert native.pgen_aa_degenerate_batch(m, allowed, threads=8) == native.pgen_aa_batch(m, seqs)


def test_batch_validates_args():
    m = load_bundled("TRG", "olga")
    allowed = [list("CATWDKQLGKKIKVF"), list("CATWDRGWDTTGWFKIF")]
    with pytest.raises(ValueError, match="same length"):
        native.pgen_aa_degenerate_batch(m, allowed, v=["TRGV9*01"])
    with pytest.raises(ValueError, match="same length"):
        native.pgen_aa_degenerate_batch(m, allowed, j=["TRGJ1*01"])
