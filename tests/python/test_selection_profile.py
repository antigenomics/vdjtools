"""Lower-dimensional normalized selection shapes on complete-history models."""
import numpy as np
import polars as pl
import pytest

from vdjtools.model import selection_reference
from vdjtools.model.selection import fit_selection_profile, selection_profile_factors
from test_joint_history import history_model


@pytest.mark.parametrize('kind', ['length', 'vj-length'])
def test_selection_profile_normalizes_and_recovers_coarse_distribution(kind):
    model = history_model()
    reference = selection_reference(model, threads=1)
    observations = reference.with_columns(
        (pl.col('p_pre') * 20000 *
         pl.when(pl.col('j_call') == 'TRBJ2-1').then(2).otherwise(1) *
         pl.when(pl.col('junction_aa_length') == 3).then(3).otherwise(1))
        .round().cast(pl.Int64).alias('events')).drop('p_pre')
    profile = fit_selection_profile(model, observations, kind=kind, smoothing=0, threads=1)
    factors = selection_profile_factors(profile, reference['junction_aa_length'].to_numpy(),
        reference['v_call'].to_list(), reference['j_call'].to_list())
    pre = reference['p_pre'].to_numpy() / reference['p_pre'].sum()
    assert np.dot(pre, factors) == pytest.approx(1, rel=1e-12)
    predicted = reference.with_columns(pl.Series('p_post', pre * factors))
    if kind == 'length':
        expected = observations.group_by('junction_aa_length').agg(pl.col('events').sum()).sort('junction_aa_length')
        got = predicted.group_by('junction_aa_length').agg(pl.col('p_post').sum()).sort('junction_aa_length')
    else:
        expected = observations
        got = predicted
    np.testing.assert_allclose(got['p_post'], expected['events'] / expected['events'].sum(), atol=2e-3)
    assert profile['audit']['supported_events'] == observations['events'].sum()
    assert profile['audit']['unsupported_events'] == 0


def test_selection_profile_retains_missing_and_ambiguous_gene_observations():
    model = history_model()
    obs = pl.DataFrame({'v_call': ['TRBV1'] * 4,
        'j_call': ['TRBJ1-1,TRBJ2-1', None, 'UNKNOWN', 'TRBJ2-1'],
        'junction_aa_length': [2, 2, 2, 100], 'events': [5, 3, 2, 1]})
    profile = fit_selection_profile(model, obs, smoothing=1, threads=1)
    assert profile['audit']['input_events'] == 11
    assert profile['audit']['supported_events'] == 10
    assert profile['audit']['unsupported_events'] == 1
    assert profile['audit']['unresolved_gene_events'] == 5
    with pytest.raises(ValueError, match='length'):
        selection_profile_factors(profile, [100], ['TRBV1'], ['TRBJ2-1'])


@pytest.mark.parametrize('kind', ['length', 'vj-length'])
def test_selected_marginal_pgen_equals_explicit_joint_sum_and_keeps_model(kind):
    from vdjtools.model import native
    from vdjtools.model.selection import selection_pgen_aa_batch
    model = history_model()
    reference = selection_reference(model, threads=1)
    obs = reference.with_columns((pl.col('p_pre') * 10000 *
        pl.when(pl.col('j_call') == 'TRBJ2-1').then(2).otherwise(1) *
        pl.when(pl.col('junction_aa_length') == 3).then(3).otherwise(1))
        .round().cast(pl.Int64).alias('events')).drop('p_pre')
    profile = fit_selection_profile(model, obs, kind=kind, threads=1)
    sequences = ['CF', 'CNF', 'CTF']
    before = native.pgen_aa_batch(model, sequences, threads=1)
    expected = []
    for seq in sequences:
        value = 0.
        for j in ('TRBJ1-1', 'TRBJ2-1'):
            factor = selection_profile_factors(profile, [len(seq)], ['TRBV1'], [j])[0]
            value += factor * native.pgen_aa(model, seq, 'TRBV1', j)
        expected.append(value)
    actual = selection_pgen_aa_batch(model, sequences, profile, threads=1)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
    assert selection_pgen_aa_batch(model, sequences, profile, threads=2) == actual
    assert native.pgen_aa_batch(model, sequences, threads=1) == before
    assert not np.allclose(actual, before, rtol=1e-6, atol=0)
    with pytest.raises(ValueError, match='model'):
        selection_pgen_aa_batch(history_model(p_nd2=.4), sequences, profile, threads=1)


def test_profile_report_alice_scales_eligible_ball_and_checks_model(monkeypatch):
    from vdjtools.model.selection import load_selection_profile
    from vdjtools.overlap.alice import _score_locus, alice
    model = history_model()
    reference = selection_reference(model, threads=1)
    profile = fit_selection_profile(model, reference.with_columns(
        (pl.col('p_pre') * 10000).round().cast(pl.Int64).alias('events')), threads=1)
    report = {'selection_profile': profile, 'selection_q': 2.}
    assert load_selection_profile(report, model)[1] == 2.
    with pytest.raises(ValueError, match='model'):
        load_selection_profile(report, history_model(p_nd2=.4))
    # Test the shared consumer with a real native constrained ball.
    import importlib
    module = importlib.import_module('vdjtools.overlap.alice')
    monkeypatch.setattr(module, '_degrees', lambda *args: np.array([3]))
    sample = pl.DataFrame({'junction_aa': ['CF'], 'v_call': ['TRBV1'],
        'j_call': ['TRBJ1-1'], 'duplicate_count': [3]})
    with pytest.raises(ValueError, match='scope'):
        alice(sample, model=model, selection_profile=report, scope='0,0,0,0')
    with pytest.raises(ValueError, match='anchored'):
        alice(sample.with_columns(pl.lit('AF').alias('junction_aa')), model=model, selection_profile=report)
    scored = _score_locus(sample, model, 'TRB', '1,0,0,1', 2., 3, 1, profile=profile)
    from vdjtools.model.native import pgen_aa_gene_batch
    ball = pgen_aa_gene_batch(model, ['CF'], ['TRBV1'], ['TRBJ1-1'],
        mismatches=1, anchored_only=True, threads=1)[0]
    q = 2 * selection_profile_factors(profile, [2], ['TRBV1'], ['TRBJ1-1'])[0]
    assert scored['selection_q'][0] == pytest.approx(q)
    assert scored['E'][0] == pytest.approx(q * ball)


def test_selection_profile_single_length_support(monkeypatch):
    import vdjtools.model.selection as module
    reference = pl.DataFrame({'v_call':['TRBV1']*2, 'j_call':['TRBJ1-1','TRBJ2-1'],
        'junction_aa_length':[2,2], 'p_pre':[.3,.2]})
    monkeypatch.setattr(module,'selection_reference',lambda *args,**kw:reference)
    obs=reference.with_columns(pl.Series('events',[30,20]))
    shape=fit_selection_profile(history_model(),obs,threads=1)
    assert len(shape['log_length'])==1
    assert shape['audit']['supported_events']==50


@pytest.mark.parametrize('species,locus', [('human','TRA'),('human','TRB'),('mouse','TRA'),('mouse','TRB')])
@pytest.mark.parametrize('kind', ['length','vj-length'])
def test_bundled_profile_matches_model_and_normalized_gene_mass(species,locus,kind):
    from vdjtools.model import load_bundled,load_selection_profile
    shape,q=load_selection_profile(kind,load_bundled(locus,organism=species,collapse=False))
    assert shape['kind']==kind and np.isfinite(q) and q>0
    assert sum(row['p'] for row in shape['selected_vj'])==pytest.approx(1,abs=1e-12)
    assert shape['audit']['input_events']==shape['audit']['supported_events']+shape['audit']['unsupported_events']
    assert np.all(np.isfinite(shape['log_length']))


def test_alice_cli_shipped_profile_records_complete_multiplier(tmp_path):
    from typer.testing import CliRunner
    from vdjtools.cli import app
    from vdjtools.model import load_bundled,load_selection_profile
    sample=tmp_path/'sample.tsv'
    pl.DataFrame({'junction_aa':['CASSF','CATSF','CATTF'], 'v_call':['TRBV19']*3,
        'j_call':['TRBJ2-1']*3,'duplicate_count':[5]*3}).write_csv(sample,separator='\t')
    out=tmp_path/'out.tsv'
    result=CliRunner().invoke(app,['tcrnet',str(sample),'--mode','alice','--locus','TRB',
        '--selection-profile','vj-length','--threads','1','--out',str(out)])
    assert result.exit_code==0,result.output
    scored=pl.read_csv(out,separator='\t')
    shape,q=load_selection_profile('vj-length',load_bundled('TRB',collapse=False))
    expected=q*selection_profile_factors(shape,scored['junction_aa'].str.len_chars().to_numpy(),
        scored['v_call'].to_list(),scored['j_call'].to_list())
    np.testing.assert_allclose(scored['selection_q'],expected,rtol=1e-12)
    np.testing.assert_allclose(scored['E'],scored['n_exposure']*scored['pgen_ball']*expected,rtol=1e-12)


def test_native_selected_views_keep_allele_weights_inside_joint_sum():
    from vdjtools.model import native
    from vdjtools._core import selection_pgen_aa_batch as selected, pgen_aa
    model = history_model(duplicate_alleles=True)
    packed, vi, ji = native.pack(model)
    weights = [(.25 + v) * (1 + j * .75) if (v + j) % 3 else 0.
               for v in range(len(vi)) for j in range(len(ji))]
    original = list(packed.pv)
    sequences = ['CF', 'CNF', 'CTF', '', 'C'] * 3
    expected = [sum(weights[v * len(ji) + j] * pgen_aa(packed, seq, v, j)
                    for v in range(len(vi)) for j in range(len(ji))) for seq in sequences]
    for threads in (1, 2):
        np.testing.assert_allclose(selected(packed, sequences, weights, threads), expected,
                                   rtol=1e-12, atol=0)
    assert packed.pv == original


def test_selected_workspace_counts_model_lengths_and_workers():
    from vdjtools.model import native
    from vdjtools._core import _selection_workspace_bytes as estimate
    packed, _, _ = native.pack(history_model())
    baseline = estimate(packed, 10, 100, 1)
    assert estimate(packed, 20, 100, 1) > baseline
    assert estimate(packed, 10, 100, 2) > baseline
    old = list(packed.q_tandem)
    try:
        packed.q_tandem = old + [0.] * 10000
        assert estimate(packed, 10, 100, 1) >= baseline + 80000
    finally:
        packed.q_tandem = old
    with pytest.raises(ValueError):
        estimate(packed, 10**9, 100, 1)
    with pytest.raises(ValueError, match='indexing'):
        estimate(packed, (2**31 - 1) // 100 // 3 + 1, 1, 1)


def test_selected_vj_view_matches_explicit_allele_sum():
    from vdjtools.model import load_bundled, native
    from vdjtools.model.generate import generate
    from vdjtools._core import selection_pgen_aa_batch as selected, pgen_aa
    model = load_bundled('TRA')
    rows = generate(model, 2, seed=226, productive_only=True, engine='native', threads=1)
    packed, vi, ji = native.pack(model)
    v = vi[rows['v_call'][0]]; j = ji[rows['j_call'][0]]
    vs = {v, (v + 1) % len(vi)}; js = {j, (j + 1) % len(ji)}
    weights = [(.5 + a) / (1 + b) if a in vs and b in js else 0.
               for a in range(len(vi)) for b in range(len(ji))]
    seqs = rows['junction_aa'].to_list()
    expected = [sum(weights[a * len(ji) + b] * pgen_aa(packed, seq, a, b)
                    for a in vs for b in js) for seq in seqs]
    assert expected[0] > 0
    np.testing.assert_allclose(selected(packed, seqs, weights, 2), expected, rtol=1e-12, atol=0)


@pytest.mark.parametrize('kind', ['length', 'vj-length'])
def test_selected_memory_guard_rejects_before_native_scoring(monkeypatch, kind):
    from vdjtools.model.selection import selection_pgen_aa_batch
    import vdjtools._core as core
    model = history_model()
    reference = selection_reference(model, threads=1)
    profile = fit_selection_profile(model, reference.with_columns(
        (pl.col('p_pre') * 10000).round().cast(pl.Int64).alias('events')),
        kind=kind, threads=1)
    monkeypatch.setattr('vdjtools.cores.available_memory', lambda: 0)
    def unexpected(*args):
        raise AssertionError('native scorer started before memory admission')
    monkeypatch.setattr(core, 'selection_pgen_aa_batch', unexpected)
    with pytest.raises(MemoryError, match='selection-weighted native Pgen'):
        selection_pgen_aa_batch(model, ['CF'], profile, threads=2)
