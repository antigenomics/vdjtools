"""Partial-codon pruning must preserve independent sums over nucleotide encodings."""
from itertools import product
import math

import pytest

from vdjtools.model import native


@pytest.mark.parametrize(("fixture", "allowed"), [
    ("toy_model_vdj", list("CGGF")), ("toy_model_vdj", list("CGDGF")),
    ("toy_model_vdj", list("CGDRGF")), ("toy_model", list("CAQF")),
    ("toy_model", list("CAGYTF")), ("toy_model", ["C", "AG", "Q", "F"]),
])
def test_partial_codons_match_nt_enumeration(request, fixture, allowed):
    # Different D lengths/trims cross all three codon-boundary phases. The NT
    # scorer enumerates rearrangements independently of the AA transfer matrix.
    codons = {
        "C": ("TGT", "TGC"), "G": ("GGT", "GGC", "GGA", "GGG"),
        "D": ("GAT", "GAC"), "F": ("TTT", "TTC"),
        "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
        "A": ("GCT", "GCC", "GCA", "GCG"), "Q": ("CAA", "CAG"),
        "T": ("ACT", "ACC", "ACA", "ACG"), "Y": ("TAT", "TAC"),
    }
    model = request.getfixturevalue(fixture)
    options = [tuple(c for aa in residues for c in codons[aa]) for residues in allowed]
    nts = ["".join(parts) for parts in product(*options)]
    expected = math.fsum(native.pgen_nt_batch(model, nts, threads=1))
    assert expected > 0
    assert native.pgen_aa_degenerate(model, allowed) == pytest.approx(expected, rel=1e-12, abs=0)
