"""Exact V/J/length reference masses for lower-dimensional selection profiles."""
from collections import defaultdict

import pytest

from vdjtools.model.reference import translate
from vdjtools.model import selection_reference
from test_joint_history import exact_outcomes, history_model


@pytest.mark.parametrize('duplicates', [False, True])
def test_selection_reference_equals_complete_history_enumeration(duplicates):
    model = history_model(duplicate_alleles=duplicates)
    expected = defaultdict(float)
    for nt, probability in exact_outcomes(model).items():
        if len(nt) % 3:
            continue
        aa = translate(nt)
        if aa.startswith('C') and aa[-1] in 'FW' and '*' not in aa:
            j = 'TRBJ1-1' if nt.endswith('TTT') else 'TRBJ2-1'
            expected[('TRBV1', j, len(aa))] += probability
    result = selection_reference(model, threads=1)
    assert set(result.columns) == {'v_call', 'j_call', 'junction_aa_length', 'p_pre'}
    observed = {(v, j, length): p for v, j, length, p in result.iter_rows()}
    assert observed == pytest.approx(dict(expected), rel=1e-12, abs=0)
    assert 0 < sum(observed.values()) <= 1
    assert selection_reference(model, threads=2).equals(result)
