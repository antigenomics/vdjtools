"""Inference-only pruning matches independent conditional models and reports omitted mass."""
from math import exp, log
from copy import deepcopy

import numpy as np
import polars as pl
import pytest

from test_insertion_ambiguity import V, J, _ambiguous_model
from vdjtools import _core
from vdjtools.model import Model, native
from vdjtools.model.dd import to_dd, to_history
from vdjtools.model.io import from_germline

_WEAK = 1e-8
_TELEMETRY = {'pruned_mass_fraction_sum', 'pruned_mass_fraction_max', 'n_pruned_observations'}


def _fixture(kind, weak):
    segments = [('V', 'V0*01', 'GG'), ('J', 'J0*01', 'CC')]
    segments += [('V', 'V1*01', 'GG')] if kind == 'vj' else [('D', 'D0*01', 'AT'), ('D', 'D1*01', 'AT')]
    gl = pl.DataFrame([dict(segment=segment, allele=allele, gene=allele.split('*')[0],
                            sequence=sequence, functional=True, cdr3_anchor=-1, full_germline='')
                       for segment, allele, sequence in segments])
    base = from_germline(gl, locus='TOY', ins_max=0,
                         palindrome_max={'v_3': 0, 'j_5': 0, 'd_5': 0, 'd_3': 0})
    tables = dict(base.tables)
    for name in ('v_3_del', 'j_5_del'):
        tables[name] = tables[name].with_columns(p=(pl.col('ndel') == 0).cast(pl.Float64))
    branch, allele = ('v_choice', 'v_allele') if kind == 'vj' else ('d_gene', 'd_allele')
    tables[branch] = tables[branch].with_columns(
        p=pl.when(pl.col(allele).str.contains('1\\*')).then(weak).otherwise(1 - weak))
    if kind != 'vj':
        tables['d_del'] = tables['d_del'].with_columns(
            p=((pl.col('ndel5') == 0) & (pl.col('ndel3') == 0)).cast(pl.Float64))
    base = Model(manifest=base.manifest, tables=tables, genomic=base.genomic).validate()
    if kind in ('legacy_dd', 'joint_dd'):
        base = to_dd(base, p_nd2=0.4)
        tables = dict(base.tables)
        tables['d2_gene'] = tables['d2_gene'].with_columns(p=(pl.col('d2_allele') == 'D0*01').cast(pl.Float64))
        base = Model(manifest=base.manifest, tables=tables, genomic=base.genomic).validate()
        if kind == 'joint_dd':
            base = to_history(base, p_nd2=0.4)
    return base


def _physical_fields(counts):
    return [name for name in dir(counts) if not name.startswith('_') and name not in _TELEMETRY]


def _run(model, sequence, rate, ratio=None, threads=1, repeats=1, tandem=True, duplicates=False):
    packed, _, _ = native.pack(model)
    counts = _core.make_counts(packed)
    # Uniform duplicate multiplicity leaves the weak posterior fraction unchanged.
    vm = [0, 1] if model.chain_type == 'VJ' else [0]
    dm = [] if model.chain_type == 'VJ' else [0, 1]
    d2 = [] if model.chain_type == 'VJ' else [0]
    if duplicates:
        vm, dm, d2 = [v for v in vm for _ in (0, 1)], [d for d in dm for _ in (0, 1)], d2 * 2
    kw = {} if ratio is None else {'prune_ratio': ratio}
    ll = _core.estep_batch(packed, [native._encode(sequence)] * repeats,
                          [vm] * repeats, [[0]] * repeats, [dm] * repeats, counts,
                          threads=threads, dd_allowed=[int(tandem)] * repeats,
                          d2masks=[d2] * repeats, error_rate=rate, **kw)
    return ll, counts


@pytest.mark.parametrize('kind,tandem', [('vj', True), ('legacy_vdj', False),
                                        ('legacy_dd', True), ('legacy_dd', False),
                                        ('joint_dd', True), ('joint_dd', False)])
@pytest.mark.parametrize('rate', [0.0, 0.07])
@pytest.mark.parametrize('duplicates', [False, True])
def test_pruned_branch_matches_independent_zero_weight_conditional_model(kind, tandem, rate, duplicates):
    model, retained = _fixture(kind, _WEAK), _fixture(kind, 0.0)
    emitted = 'GGCC' if kind == 'vj' else ('GGATATCC' if kind.endswith('dd') and tandem else 'GGATCC')
    sequence = 'A' + emitted[1:] if rate else emitted
    before = native.pgen_nt(model, emitted)
    snapshots = {name: frame.clone() for name, frame in model.tables.items()}
    germlines = {name: frame.clone() for name, frame in model.genomic.items()}
    manifest = deepcopy(model.manifest)
    reference_ll, reference = _run(retained, sequence, rate, 0, tandem=tandem, duplicates=duplicates)
    exact_ll, exact = _run(model, sequence, rate, 0, tandem=tandem, duplicates=duplicates)
    default_ll, default = _run(model, sequence, rate, tandem=tandem, duplicates=duplicates)
    assert default_ll == exact_ll
    for name in _physical_fields(exact):
        assert np.array_equal(np.asarray(getattr(default, name)), np.asarray(getattr(exact, name))), name
    assert all(getattr(default, name) == getattr(exact, name) == 0 for name in _TELEMETRY)
    for threads, repeats in ((1, 1), (2, 64)):
        ll, actual = _run(model, sequence, rate, 1e-5, threads, repeats, tandem, duplicates)
        assert ll / repeats == pytest.approx(reference_ll + log(1 - _WEAK), abs=1e-12)
        for name in _physical_fields(reference):
            assert np.asarray(getattr(actual, name)) / repeats == pytest.approx(
                np.asarray(getattr(reference, name)), abs=1e-12, rel=1e-11), name
        omitted = 1 - exp(ll / repeats - exact_ll)
        assert omitted == pytest.approx(_WEAK, abs=1e-14)
        assert actual.pruned_mass_fraction_sum / repeats == pytest.approx(omitted, abs=1e-14)
        assert actual.pruned_mass_fraction_max == pytest.approx(omitted, abs=1e-14)
        assert actual.n_pruned_observations == repeats
    assert native.pgen_nt(model, emitted) == before
    assert all(model.tables[name].equals(frame) for name, frame in snapshots.items())
    assert all(model.genomic[name].equals(frame) for name, frame in germlines.items())
    assert model.manifest == manifest


def test_low_ratio_keeps_weak_branch_and_coincident_n_insertion_history():
    model = _fixture('vj', _WEAK)
    exact_ll, exact = _run(model, 'GGCC', 0, 0)
    ll, actual = _run(model, 'GGCC', 0, 1e-12)
    assert ll == pytest.approx(exact_ll, abs=1e-14)
    for name in _physical_fields(exact):
        assert np.asarray(getattr(actual, name)) == pytest.approx(np.asarray(getattr(exact, name)))
    assert actual.pruned_mass_fraction_sum == actual.pruned_mass_fraction_max == actual.n_pruned_observations == 0
    packed, _, _ = native.pack(_ambiguous_model())
    counts = _core.make_counts(packed)
    _core.estep_batch(packed, [native._encode(V + J)], [], [], [], counts, prune_ratio=1e-5)
    assert counts.ins_vj == pytest.approx([16 / 17, 0, 1 / 17])
    assert counts.v_3_del[2] == pytest.approx(1 / 17)
    assert counts.n_pruned_observations == 0


def test_weak_component_inside_a_retained_vdj_cell_is_not_individually_removed():
    base = _fixture('legacy_vdj', 0)
    tables, genomic = dict(base.tables), dict(base.genomic)
    genomic['genes_v'] = pl.concat([genomic['genes_v'], genomic['genes_v'].with_columns(
        v_allele=pl.lit('V1*01'), gene=pl.lit('V1'))])
    tables['v_choice'] = pl.DataFrame({'v_allele': ['V0*01', 'V1*01'], 'p': [1 - _WEAK, _WEAK]})
    tables['v_3_del'] = pl.concat([tables['v_3_del'], tables['v_3_del'].with_columns(v_allele=pl.lit('V1*01'))])
    model = Model(manifest=base.manifest, tables=tables, genomic=genomic).validate()
    # V alternatives share one D-realization cell, so its strong component retains both.
    packed, _, _ = native.pack(model)
    exact, actual = _core.make_counts(packed), _core.make_counts(packed)
    sequence = [native._encode('GGATCC')]
    exact_ll = _core.estep_batch(packed, sequence, [], [], [], exact)
    ll = _core.estep_batch(packed, sequence, [], [], [], actual, prune_ratio=1e-5)
    assert ll == pytest.approx(exact_ll, abs=1e-14)
    for name in _physical_fields(exact):
        assert np.asarray(getattr(actual, name)) == pytest.approx(np.asarray(getattr(exact, name)))
    assert actual.v_choice[1] == pytest.approx(_WEAK)
    assert actual.n_pruned_observations == actual.pruned_mass_fraction_sum == 0


@pytest.mark.parametrize('ratio', [-1e-5, float('nan'), float('inf'), 1.0])
def test_invalid_pruning_ratio_is_rejected(ratio):
    with pytest.raises(ValueError, match='prune_ratio'):
        _run(_fixture('vj', _WEAK), 'GGCC', 0, ratio)


def test_high_level_default_records_actual_pruning_in_reports_and_training_frames():
    from vdjtools.model.infer import InferenceReport, infer_native, training_frame
    observations = ['GGCC', 'GGCC', 'AAAA']  # the unsupported read must not dilute mean omitted mass
    fitted, report = infer_native(_fixture('vj', _WEAK), observations,
                                   init='template', max_iter=1, threads=1)
    assert report.prune_ratio == 1e-5
    assert report.n_sequences == 3 and report.n_scoreable == [2]
    assert report.loglik == pytest.approx([log(1 - _WEAK)], abs=1e-14)
    assert report.pruned_mass_mean == pytest.approx([_WEAK], abs=1e-14)
    assert report.pruned_mass_max == pytest.approx([_WEAK], abs=1e-14)
    assert report.n_pruned_observations == [2]
    frame = report.to_frame()
    assert InferenceReport.from_dict(report.to_dict()).to_frame().equals(frame)
    assert training_frame(report).drop('run').equals(frame)
    assert training_frame(fitted).drop('run').equals(frame)
    assert frame['pruned_mass_mean'].to_list() == report.pruned_mass_mean
    assert frame['pruned_mass_max'].to_list() == report.pruned_mass_max
    assert frame['n_pruned_observations'].to_list() == report.n_pruned_observations


@pytest.mark.parametrize('ratio', [0.0, 1e-12, 1e-5])
def test_resume_restores_saved_pruning_policy_and_preserves_run_history(tmp_path, ratio):
    from vdjtools.model.infer import infer_native, resume, training_frame
    observations = ['GGCC']
    fitted, first = infer_native(_fixture('vj', _WEAK), observations, init='template',
                                 max_iter=1, threads=1, prune_ratio=ratio)
    checkpoint = tmp_path / 'synthetic-checkpoint'
    fitted.save(checkpoint)
    resumed, second = resume(checkpoint, observations, max_iter=1, threads=1, checkpoint=None)
    assert second.prune_ratio == ratio
    assert [run['prune_ratio'] for run in resumed.training['runs']] == [ratio, ratio]
    combined = training_frame(resumed)
    assert combined['run'].to_list() == [0, 1]
    if ratio:
        assert combined['pruned_mass_mean'].to_list() == first.pruned_mass_mean + second.pruned_mass_mean
        assert combined['pruned_mass_max'].to_list() == first.pruned_mass_max + second.pruned_mass_max
        assert combined['n_pruned_observations'].to_list() == first.n_pruned_observations + second.n_pruned_observations
    else:
        assert 'pruned_mass_mean' not in combined.columns
    exact, explicit = resume(fitted, observations, max_iter=1, threads=1, prune_ratio=0)
    assert explicit.prune_ratio == 0
    assert exact.training['runs'][-1]['prune_ratio'] == 0
    if ratio:
        # Exact continuation has no pruning columns; diagonal history keeps the earlier diagnostics.
        assert training_frame(exact)['pruned_mass_mean'].to_list() == first.pruned_mass_mean + [None]


def test_resume_of_a_historical_run_without_pruning_metadata_stays_exact():
    from vdjtools.model.infer import infer_native, resume
    fitted, _ = infer_native(_fixture('vj', _WEAK), ['GGCC'], init='template',
                             max_iter=1, threads=1, prune_ratio=0)
    history = deepcopy(fitted.training)
    for key in ('prune_ratio', 'pruned_mass_mean', 'pruned_mass_max', 'n_pruned_observations'):
        history['runs'][-1].pop(key)
    historical = Model(manifest=fitted.manifest, tables=fitted.tables,
                       genomic=fitted.genomic, training=history)
    _, report = resume(historical, ['GGCC'], max_iter=1, threads=1)
    assert report.prune_ratio == 0
    assert report.pruned_mass_mean == report.pruned_mass_max == [0]
    assert report.n_pruned_observations == [0]
    assert report.loglik == pytest.approx([0], abs=1e-14)


def test_python_frame_inference_allows_exact_policy_and_rejects_pruning():
    from vdjtools.model.infer import infer_frame
    frame = pl.DataFrame({'junction': ['GGCC']})
    _, report = infer_frame(_fixture('vj', _WEAK), frame, native=False, use_calls=False,
                            init='template', max_iter=1, prune_ratio=0)
    assert report.prune_ratio == 0 and not report.native
    with pytest.raises(ValueError, match='prun'):
        infer_frame(_fixture('vj', _WEAK), frame, native=False, use_calls=False,
                    init='template', max_iter=1, prune_ratio=1e-5)
