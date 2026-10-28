"""The fixed-D AA walk must equal an independent sum over nucleotide encodings."""
from itertools import product
import math

import pytest

from vdjtools.model import native


@pytest.mark.parametrize("aa", ["CGGF", "CGDGF", "CGDRGF"])
def test_fixed_d_prefixes_match_nt_enumeration(toy_model_vdj, aa):
    # Different D lengths/trims cross all three codon-boundary phases. The NT
    # scorer enumerates rearrangements independently of the AA transfer matrix.
    codons = {
        "C": ("TGT", "TGC"), "G": ("GGT", "GGC", "GGA", "GGG"),
        "D": ("GAT", "GAC"), "F": ("TTT", "TTC"),
        "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    }
    nts = ["".join(parts) for parts in product(*(codons[a] for a in aa))]
    expected = math.fsum(native.pgen_nt_batch(toy_model_vdj, nts, threads=1))
    assert expected > 0
    assert native.pgen_aa(toy_model_vdj, aa) == pytest.approx(expected, rel=1e-12, abs=0)
