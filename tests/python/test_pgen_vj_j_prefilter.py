"""Exact VJ J-prefilter oracles against the unchanged nucleotide enumeration."""
import itertools

import pytest

from vdjtools import _core
from vdjtools.model.reference import _CODON_TABLE


def _encode(nt):
    return ["ACGT".index(base) for base in nt]


def _model(v, j, del_j):
    model = _core.PackedModel()
    model.vdj = False
    model.cut_v, model.cut_j = [_encode(v)], [_encode(j)]
    model.func_v, model.func_j = [0], [0]
    model.pv, model.pjv = [1.0], [1.0]
    model.nbins_v, model.nbins_j = len(v) + 1, len(j) + 1
    model.del_v = [1.0] + [0.0] * len(v)
    model.del_j = del_j
    model.ins_vj = [1.0]
    model.R_vj, model.bias_vj = [0.25] * 16, [0.25] * 4
    return model


@pytest.mark.parametrize("v,j,deletions,aa,expected", [
    ("TGT", "GCTTTT", {0: 1.0}, "CAF", 1.0),
    ("TGT", "GCTTTT", {0: 1.0}, "CAY", 0.0),  # incompatible full terminal codon
    ("TGT", "GCTTTT", {}, "CAF", 0.0),  # no nonzero J deletion
    ("TGT", "GCTTTT", {6: 1.0}, "CAF", 0.0),  # fully deleted J cannot contribute
    ("TGTGCT", "ATGTTT", {0: 0.5, 3: 0.5}, "CAF", 0.5),  # shortest live trim wins
    ("TG", "TATGTTT", {0: 1.0}, "CMF", 1.0),  # one J nt crosses the first codon
    ("T", "GTATGTTT", {0: 1.0}, "CMF", 1.0),  # two J nt cross the first codon
    ("TGTATGTT", "T", {0: 1.0}, "CMF", 1.0),  # no full J codon: cannot reject
    ("TGTATGT", "TT", {0: 1.0}, "CMF", 1.0),
])
def test_vj_j_prefilter_matches_unoptimized_nt(v, j, deletions, aa, expected):
    row = [deletions.get(index, 0.0) for index in range(len(j) + 1)]
    model = _model(v, j, row)
    synonymous = [[codon for codon, residue in _CODON_TABLE.items() if residue == a] for a in aa]
    # VJ pgen_nt uses the existing direct enumeration, never pgen_aa_vj.
    oracle = sum(_core.pgen_nt(model, _encode("".join(codons)))
                 for codons in itertools.product(*synonymous))
    assert oracle == expected
    assert _core.pgen_aa(model, aa) == oracle
    assert _core.pgen_aa(model, aa, 0, 0) == oracle


def test_vj_j_prefilter_uses_allowed_codon_masks():
    model = _model("TGT", "GCTTTT", [1.0] + [0.0] * 6)
    assert _core.pgen_aa_degenerate(model, ["C", "A", "FY"]) == 1.0
