"""One projection implementation for statistical, geometry and RFF feature tables."""
import numpy as np
import polars as pl
import pytest
from vdjtools.signature.corpus import TablePCA


@pytest.mark.parametrize('prefix',['vsig','rsig','phi'])
def test_feature_table_pca_roundtrip_and_channel_passthrough(tmp_path,prefix):
    cols = [f'{prefix}:feature:{i}' for i in range(4)]
    table = pl.from_numpy(np.random.default_rng(4).normal(size=(12,4)),schema=cols)
    table = table.with_columns(pl.Series('sample_id',[f'{i:03}' for i in range(12)]),
                               pl.lit(prefix+':fixed').alias('embedding_space'),
                               pl.Series('channel',list(range(12))))
    reducer = TablePCA.fit(table,columns=cols,n_components=2)
    reducer.save(tmp_path/'rotation.npz')
    expected = reducer.transform(table)
    got = TablePCA.load(tmp_path/'rotation.npz').transform(table)
    assert got.equals(expected)
    assert got['sample_id'].to_list() == table['sample_id'].to_list()
    assert got['channel'].to_list() == table['channel'].to_list()
    with pytest.raises(ValueError,match='space'):
        reducer.transform(table.with_columns(pl.lit('other').alias('embedding_space')))
    with pytest.raises(ValueError,match='columns'):
        reducer.transform(table.drop(cols[0]))


def test_raw_features_need_no_corpus_and_match_existing_featuriser():
    from vdjtools.signature.features import gene_vocab, raw_and_channels
    from vdjtools.signature.signature import raw_features_cohort
    df = pl.DataFrame({'v_call':['TRBV10-3*01','TRBV20-1*01'],
                       'j_call':['TRBJ2-7*01','TRBJ1-2*01'],
                       'junction_aa':['CASSIRSSYEQYF','CSARVSGYYGYTF'],
                       'duplicate_count':[12,4]})
    vocab = {'TRB':gene_vocab('TRB')}
    expected,channels = raw_and_channels({'TRB':df},vocab,cstar_target=None)
    expected.update(channels)
    got = raw_features_cohort([('001',lambda:{'TRB':df})],vocab=vocab,cstar_target=None,n_jobs=1)
    assert got['sample_id'].to_list() == ['001']
    assert not any(':pc:' in c for c in got.columns)
    for col,value in expected.items():
        assert col in got.columns
        np.testing.assert_equal(got[col][0],value)


@pytest.mark.parametrize('workers',[1,2])
def test_raw_feature_workers_release_input_frames(workers):
    import weakref
    from threading import Lock
    from vdjtools.cores import available_cores, kernel_threads
    from vdjtools.signature.signature import raw_features_cohort
    from vdjtools.signature.features import gene_vocab
    if workers > available_cores():
        pytest.skip('needs two allocated cores')
    state = {'live':0,'peak':0}
    lock = Lock()
    def release():
        with lock:
            state['live'] -= 1
    def read():
        assert kernel_threads(0) == 1
        df = pl.DataFrame({'v_call':['TRBV10-3*01','TRBV20-1*01'],
                           'j_call':['TRBJ2-7*01','TRBJ1-2*01'],
                           'junction_aa':['CASSIRSSYEQYF','CSARVSGYYGYTF'],
                           'duplicate_count':[12,4]})
        with lock:
            state['live'] += 1
            state['peak'] = max(state['peak'],state['live'])
        weakref.finalize(df,release)
        return {'TRB':df}
    out = raw_features_cohort([(str(i),read) for i in range(16)],
                              vocab={'TRB':gene_vocab('TRB')},n_jobs=workers,cstar_target=None)
    assert out.height == 16
    assert state['live'] == 0 and state['peak'] <= workers


@pytest.mark.benchmark
@pytest.mark.skipif(not __import__('os').getenv('RUN_BENCHMARK'),reason='set RUN_BENCHMARK=1')
def test_large_raw_feature_cohort_memory():
    import os
    import resource
    import sys
    import time
    from threadpoolctl import threadpool_info
    from vdjtools.cores import kernel_threads
    from vdjtools.signature.signature import raw_features_cohort
    from vdjtools.signature.features import gene_vocab
    count = int(os.getenv('VDJ_BENCH_SAMPLES','10000'))
    workers = int(os.getenv('VDJ_BENCH_WORKERS','1'))
    def read():
        assert kernel_threads(0) == 1
        assert all(p['num_threads'] == 1 for p in threadpool_info())
        return {'TRB':pl.DataFrame({'v_call':['TRBV10-3*01','TRBV20-1*01'],
                                   'j_call':['TRBJ2-7*01','TRBJ1-2*01'],
                                   'junction_aa':['CASSIRSSYEQYF','CSARVSGYYGYTF'],
                                   'duplicate_count':[12,4]})}
    start = time.perf_counter()
    out = raw_features_cohort([(str(i),read) for i in range(count)],
                              vocab={'TRB':gene_vocab('TRB')},n_jobs=workers,cstar_target=None)
    assert out.height == count
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    gib = rss/(1024**3 if sys.platform == 'darwin' else 1024**2)
    print(f'Raw feature cohort: samples={count}, columns={out.width}, workers={workers}, '
          f'wall_s={time.perf_counter()-start:.3f}, peak_RSS_GiB={gib:.3f}',flush=True)
    assert gib < 16


def test_table_pca_saves_variance_and_loads_legacy(tmp_path):
    from vdjtools.signature.corpus import fit_pca
    cols = [f'phi:{i}' for i in range(6)]
    values = np.random.default_rng(17).normal(size=(20,6))
    table = pl.from_numpy(values,schema=cols).with_columns(pl.lit('test').alias('embedding_space'))
    expected = fit_pca(values.copy(),4).steps[1][1].explained_variance_ratio_
    rotation = TablePCA.fit(table,n_components=4)
    np.testing.assert_array_equal(rotation.explained_variance_ratio,expected)
    path = tmp_path/'rotation.npz'
    rotation.save(path)
    restored = TablePCA.load(path)
    np.testing.assert_array_equal(restored.explained_variance_ratio,expected)
    with np.load(path,allow_pickle=False) as z:
        legacy = {key:z[key] for key in z.files if key != 'explained_variance_ratio'}
    np.savez(tmp_path/'legacy.npz',**legacy)
    assert TablePCA.load(tmp_path/'legacy.npz').explained_variance_ratio is None
    legacy['explained_variance_ratio'] = np.array([0.5,-0.1,0.4,0.2])
    np.savez(tmp_path/'invalid.npz',**legacy)
    with pytest.raises(ValueError,match='variance'):
        TablePCA.load(tmp_path/'invalid.npz')


def test_pca_phase_budget_restores_blas_and_rejects_nested_workers(monkeypatch):
    from threadpoolctl import threadpool_info, threadpool_limits
    from vdjtools.cores import available_cores, numerical_kernels, serial_kernels
    from vdjtools.signature import corpus as C
    budget = min(2, available_cores())
    table = pl.from_numpy(np.random.default_rng(7).normal(size=(12, 4)),
                          schema=[f'phi:{i}' for i in range(4)])
    table = table.with_columns(pl.lit('fixed').alias('embedding_space'))
    original = C.fit_pca
    seen = []
    def checked(*args, **kwargs):
        seen.append([p['num_threads'] for p in threadpool_info() if p['user_api'] == 'blas'])
        # Accelerate on macOS exposes no threadpoolctl-managed BLAS pool.
        assert all(t == budget for t in seen[-1])
        return original(*args, **kwargs)
    with threadpool_limits(limits=1, user_api='blas'):
        reference = TablePCA.fit(table, n_components=2)
        monkeypatch.setattr(C, 'fit_pca', checked)
        parallel = TablePCA.fit(table, n_components=2, threads=budget)
        assert all(p['num_threads'] == 1 for p in threadpool_info() if p['user_api'] == 'blas')
        with pytest.raises(RuntimeError, match='phase failed'):
            with numerical_kernels(budget):
                raise RuntimeError('phase failed')
        assert all(p['num_threads'] == 1 for p in threadpool_info() if p['user_api'] == 'blas')
    np.testing.assert_allclose(reference.components, parallel.components, rtol=0, atol=1e-12)
    if budget > 1:
        with serial_kernels(), pytest.raises(ValueError, match='between'):
            TablePCA.fit(table, threads=budget)
    with pytest.raises(ValueError, match='between'):
        TablePCA.fit(table, threads=available_cores() + 1)


def test_corpus_map_finishes_before_blas_fit(monkeypatch):
    from threading import Lock
    from threadpoolctl import threadpool_info
    from vdjtools.cores import available_cores, kernel_threads
    from vdjtools.signature import corpus as C, feature_definition
    if available_cores() < 2:
        pytest.skip('needs two allocated cores')
    frame = pl.DataFrame({'v_call': ['TRBV10-3*01', 'TRBV20-1*01'],
                          'j_call': ['TRBJ2-7*01', 'TRBJ1-2*01'],
                          'junction_aa': ['CASSIRSSYEQYF', 'CSARVSGYYGYTF'],
                          'duplicate_count': [12, 4]})
    completed = []
    lock = Lock()
    def read(i):
        assert kernel_threads(0) == 1
        with lock:
            completed.append(i)
        return {'TRB': frame.with_columns(pl.Series('duplicate_count', [12 + 3*i, 4 + i*i]))}
    original = C.fit_locus
    def checked(*args, **kwargs):
        assert sorted(completed) == list(range(4))
        assert kernel_threads(2) == 2
        assert all(p['num_threads'] == 2 for p in threadpool_info() if p['user_api'] == 'blas')
        return original(*args, **kwargs)
    monkeypatch.setattr(C, 'fit_locus', checked)
    result = C.fit_cohort([(str(i), lambda i=i: read(i)) for i in range(4)],
                          loci=('TRB',), feature_definition=feature_definition('basic'),
                          cstar_target=None, n_jobs=2, threads=2, n_components=2)
    assert result.fits
