"""Worker budgets must be effective before numerical runtimes import."""
import os

import pytest

from vdjtools.cores import _KERNEL_ENV, available_cores
from vdjtools.signature.cohort import parallel_rows


def _probe(i):
    import numpy as np
    import polars as pl
    try:
        from threadpoolctl import threadpool_info
    except ImportError:
        threadpool_info = list

    np.ones((4, 4)) @ np.ones((4, 4))
    return {"i": i, "pid": os.getpid(), "polars": pl.thread_pool_size(),
            "limits": {k: os.environ.get(k) for k in _KERNEL_ENV},
            "blas": [p['num_threads'] for p in threadpool_info()]}


@pytest.mark.skipif(available_cores() < 2, reason="requires two allocated cores")
def test_workers_limit_kernels_before_import_and_restore_parent(monkeypatch):
    for key in _KERNEL_ENV:
        monkeypatch.setenv(key, "4")
    rows = parallel_rows(range(4), _probe, 2)
    assert [r['i'] for r in rows] == list(range(4))
    assert all(r['pid'] == os.getpid() for r in rows)
    assert all(r['polars'] == 1 for r in rows)
    assert all(set(r['limits'].values()) == {'4'} for r in rows)
    assert all(all(n == 1 for n in r['blas']) for r in rows)
    assert all(os.environ[k] == '4' for k in _KERNEL_ENV)


@pytest.mark.skipif(available_cores() < 2, reason="requires two allocated cores")
def test_pool_failure_restores_environment(monkeypatch):
    monkeypatch.delenv('OMP_THREAD_LIMIT', raising=False)
    monkeypatch.setenv('POLARS_MAX_THREADS', '7')
    def broken(**kwargs):
        raise RuntimeError('cannot start')
    monkeypatch.setattr('concurrent.futures.ThreadPoolExecutor', broken)
    with pytest.raises(RuntimeError, match='cannot start'):
        parallel_rows([1, 2], abs, 2)
    assert os.environ['POLARS_MAX_THREADS'] == '7'
    assert 'OMP_THREAD_LIMIT' not in os.environ


def test_negative_workers_raise():
    with pytest.raises(ValueError, match='non-negative'):
        parallel_rows([1, 2], abs, -1)


def test_cli_sets_limits_before_polars_import():
    import json
    import subprocess
    import sys

    code = '''
import sys, json
from vdjtools._cli import main
assert 'polars' not in sys.modules
sys.argv = ['vdjtools', 'signature', '--help']
try:
    main()
except SystemExit as e:
    assert e.code == 0
import polars
print(json.dumps({'threads': polars.thread_pool_size()}))
'''
    env = dict(os.environ, POLARS_MAX_THREADS='4')
    result = subprocess.run([sys.executable, '-c', code], env=env, text=True,
                            capture_output=True, check=True)
    assert json.loads(result.stdout.splitlines()[-1])['threads'] == 1


def _invalid_sample(_):
    raise RuntimeError('invalid sample calculation')


@pytest.mark.skipif(available_cores() < 2, reason="requires two allocated cores")
def test_sample_runtime_errors_are_not_reported_as_pool_startup_errors():
    with pytest.raises(RuntimeError, match='^invalid sample calculation$'):
        parallel_rows([1, 2], _invalid_sample, 2)


@pytest.mark.skipif(available_cores() < 2, reason="requires two allocated cores")
def test_threads_share_resources_and_use_distinct_workers():
    from threading import Barrier, get_ident
    from vdjtools.cores import kernel_threads

    barrier = Barrier(2, timeout=10)
    resource = object()

    def work(i):
        barrier.wait()
        return i, get_ident(), resource, kernel_threads(0)

    rows = parallel_rows(range(2), work, 2)
    assert [r[0] for r in rows] == [0, 1]
    assert len({r[1] for r in rows}) == 2
    assert all(r[2] is resource and r[3] == 1 for r in rows)
    assert kernel_threads(0) >= 2


@pytest.mark.skipif(available_cores() < 2, reason="requires two allocated cores")
def test_independent_concurrent_pools_keep_their_own_callbacks():
    from concurrent.futures import ThreadPoolExecutor
    from threading import Barrier

    barrier = Barrier(4, timeout=10)
    def run(label):
        def work(i):
            barrier.wait()
            return label, i
        return parallel_rows(range(2), work, 2)

    with ThreadPoolExecutor(max_workers=2) as ex:
        rows = list(ex.map(run, ('a', 'b')))
    assert rows == [[('a', 0), ('a', 1)], [('b', 0), ('b', 1)]]


def test_one_worker_is_serial_without_starting_a_pool(monkeypatch):
    from vdjtools.cores import kernel_threads
    def forbidden(**kw):
        raise AssertionError('serial mode started a pool')
    monkeypatch.setattr('concurrent.futures.ThreadPoolExecutor', forbidden)
    assert parallel_rows([0, 1], lambda _: kernel_threads(0), 1) == [1, 1]
