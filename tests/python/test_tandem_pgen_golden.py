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


@pytest.mark.parametrize('locus', ['TRB', 'TRD', 'IGH'])
def test_tandem_worker_scratch_mixed_lengths_and_models(locus):
    golden = json.loads((Path(__file__).parent / 'fixtures' / 'tandem_pgen_golden.json').read_text())[locus]
    model = load_bundled(locus)
    # Long, short and empty queries exercise every reused worker buffer.
    order = sorted(range(len(golden['data']['junction_aa'])),
                   key=lambda i: -len(golden['data']['junction_aa'][i]))
    sequences = [golden['data']['junction_aa'][i] for i in order]
    expected = [golden['pgen']['aa_False'][i] for i in order]
    sequences = (sequences + ['', 'C', 'WWWWWW']) * 4
    expected = (expected + [0.0, native.pgen_aa(model, 'C'), native.pgen_aa(model, 'WWWWWW')]) * 4
    for threads in (1, 2):
        actual = native.pgen_aa_batch(model, sequences, threads=threads)
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
    other = load_bundled('TRA')
    assert native.pgen_aa_batch(other, [''], threads=2) == [0.0]
    np.testing.assert_allclose(native.pgen_aa_batch(model, sequences, threads=2),
                               expected, rtol=1e-12, atol=0)
