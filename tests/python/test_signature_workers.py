"""Worker budgets must be effective before numerical runtimes import."""
import os

import pytest

from vdjtools.signature.cohort import _KERNEL_ENV, WORKER_ENV, parallel_rows


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


def test_workers_limit_kernels_before_import_and_restore_parent(monkeypatch):
    for key in _KERNEL_ENV:
        monkeypatch.setenv(key, "4")
    monkeypatch.delenv(WORKER_ENV, raising=False)
    rows = parallel_rows(range(4), _probe, 2)
    assert [r['i'] for r in rows] == list(range(4))
    assert all(r['pid'] != os.getpid() for r in rows)
    assert all(r['polars'] == 1 for r in rows)
    assert all(set(r['limits'].values()) == {'1'} for r in rows)
    assert all(all(n == 1 for n in r['blas']) for r in rows)
    assert all(os.environ[k] == '4' for k in _KERNEL_ENV)
    assert WORKER_ENV not in os.environ


def test_pool_failure_restores_environment(monkeypatch):
    monkeypatch.delenv('OMP_THREAD_LIMIT', raising=False)
    monkeypatch.setenv('POLARS_MAX_THREADS', '7')
    def broken(**kwargs):
        raise RuntimeError('cannot start')
    monkeypatch.setattr('concurrent.futures.ProcessPoolExecutor', broken)
    with pytest.raises(RuntimeError, match='NOT falling back'):
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


def test_sample_runtime_errors_are_not_reported_as_pool_startup_errors():
    with pytest.raises(RuntimeError, match='^invalid sample calculation$'):
        parallel_rows([1, 2], _invalid_sample, 2)
