"""Stationary insertion biases stay probabilities before either E-step consumes them."""
from collections import defaultdict

import numpy as np
import polars as pl
import pytest

from vdjtools.model import Model, native, pgen
from vdjtools.model.infer import _estep_seq
from test_dd import _tiny_dd_model


# Finite, nonnegative, column-stochastic; LAPACK can return a negative near-zero
# component for its eigenvalue-one vector when transitions span many orders.
_TRANSITIONS = np.array([
    [.9993778834713533, .9986084122662894, 5.474480505410291e-197, 9.532948939222661e-169],
    [1.119769178185351e-33, 4.945032347457621e-63, 7.362428235195025e-45, 3.754886845524904e-25],
    [2.4306160565892733e-33, .001391587733710556, 1., 8.335636005866498e-56],
    [.0006221165286467268, 1.629894610102337e-94, 4.5842275442411615e-65, 1.],
])


def test_stationary_bias_is_a_probability_for_ill_conditioned_chain():
    assert np.isfinite(_TRANSITIONS).all() and (_TRANSITIONS >= 0).all()
    np.testing.assert_allclose(_TRANSITIONS.sum(axis=0), 1., rtol=0, atol=2e-16)
    bias = pgen._steady_state(_TRANSITIONS)
    assert np.isfinite(bias).all() and (bias >= 0).all()
    assert bias.sum() == pytest.approx(1., abs=2e-16)
    np.testing.assert_allclose(_TRANSITIONS @ bias, bias, rtol=0, atol=4e-16)


def test_valid_stationary_bias_remains_exactly_unchanged():
    matrix = np.array([[.1, .2, .3, .4], [.2, .3, .4, .1],
                       [.3, .4, .1, .2], [.4, .1, .2, .3]])
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    old = np.real(eigenvectors[:, np.argmin(np.abs(eigenvalues - 1.))])
    old /= old.sum()
    assert (old > 0).all()
    np.testing.assert_array_equal(pgen._steady_state(matrix), old)


def test_positive_model_yields_nonnegative_native_and_reference_counts():
    from vdjtools import _core
    model = _tiny_dd_model(p_nd2=1., pdel_full_only=False, dd_ins1=.5)
    tables = dict(model.tables)
    tables['d_del'] = tables['d_del'].with_columns(p=pl.Series([1., 0., 0.]))
    tables['dd_dinucl'] = tables['dd_dinucl'].with_columns(
        p=pl.Series([_TRANSITIONS[to, frm] for frm in range(4) for to in range(4)]))
    genomic = dict(model.genomic)
    genomic['genes_d'] = genomic['genes_d'].with_columns(cut_segment=pl.lit('GT'))
    model = Model(manifest=model.manifest, tables=tables, genomic=genomic)
    assert all(table['p'].is_finite().all() and (table['p'] >= 0).all() for table in tables.values())
    sequence = 'GGGTGTCC'
    reference = defaultdict(lambda: defaultdict(float))
    probability = _estep_seq(pgen.prepare(model), sequence, reference)
    packed, _, _ = native.pack(model)
    counts = _core.make_counts(packed)
    likelihood = _core.estep_batch(packed, [native._encode(sequence)], [], [], [], counts, 1)
    assert likelihood == pytest.approx(np.log(probability))
    assert all(value >= 0 for event in reference.values() for value in event.values())
    for name in ('d2_gene', 'd2_del', 'ins_dd', 'dinucl_dd'):
        values = np.asarray(getattr(counts, name))
        assert np.isfinite(values).all() and (values >= 0).all()
    expected = np.array([reference['dd_ins'].get((length,), 0.) for length in range(2)])
    np.testing.assert_allclose(counts.ins_dd, expected, rtol=1e-12, atol=0)
