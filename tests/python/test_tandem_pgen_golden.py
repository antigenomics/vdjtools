"""Complete pre-226 Pgen on fixed native-generated bundled-model junctions.

Reference probabilities were captured from commit 528baacc, before contraction
reordering. Both contraction orders, codon phases, singleton NT masks and V/J
marginalization must retain all tandem mass (rtol 1e-12, no absolute tolerance).
"""
import json
from pathlib import Path

import numpy as np
import pytest

from vdjtools.model import load_bundled, native


@pytest.mark.parametrize('locus', ['TRB', 'TRD', 'IGH'])
def test_complete_tandem_pgen_golden(locus):
    golden = json.loads((Path(__file__).parent / 'fixtures' / 'tandem_pgen_golden.json').read_text())[locus]
    model = load_bundled(locus)
    assert native.pgen_aa(model, '') == 0
    assert native.pgen_aa_degenerate(model, []) == 0
    data = golden['data']
    for kind in ('aa', 'nt'):
        for masked in (False, True):
            expected = golden['pgen'][f'{kind}_{masked}']
            assert all(p > 0 for p in expected)
            kw = dict(v=data['v_call'], j=data['j_call']) if masked else {}
            fn = getattr(native, f'pgen_{kind}_batch')
            actual = fn(model, data[f'junction_{kind}'], threads=1, **kw)
            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
            assert fn(model, data[f'junction_{kind}'], threads=2, **kw) == actual
