"""Apple's immutable ceiling and thread-local policies, with no scientific fits."""
import os
import subprocess
import sys
from types import SimpleNamespace

import pytest

from vdjtools import cores


def test_apple_policy_restores_on_failure_and_rejects_other_ceilings(monkeypatch):
    policy = [1]

    def set_policy(value):
        policy[0] = value
        return 0

    library = SimpleNamespace(BLASGetThreading=lambda: policy[0], BLASSetThreading=set_policy)
    monkeypatch.setattr(cores, '_accelerate', lambda: library)
    monkeypatch.setattr(cores, '_uses_accelerate', lambda: True)
    monkeypatch.setattr(cores, '_ACCELERATE_CAP', 2)
    with pytest.raises(RuntimeError, match='work failed'):
        with cores._accelerate_threads(2):
            assert policy[0] == 0
            raise RuntimeError('work failed')
    assert policy[0] == 1
    with pytest.raises(RuntimeError, match='fresh process'):
        with cores._accelerate_threads(3):
            pytest.fail('unbounded ceiling accepted')
    assert policy[0] == 1


@pytest.mark.parametrize('flags', [[], ['--threads=2'], ['-t', '2'], ['-t2']])
def test_cli_selects_ceiling_before_command_imports(monkeypatch, flags):
    seen = []
    monkeypatch.setattr(cores, 'configure_runtime', lambda **kw: seen.append(kw['blas_threads']))
    cores.configure_cli_runtime(['pca', *flags, '--output', 'example.tsv'])
    assert seen == ([None] if not flags else [2])


@pytest.mark.skipif(cores.available_cores() < 2, reason='requires two allocated CPUs')
def test_cli_pca_configures_its_selected_ceiling(tmp_path):
    import numpy as np
    import polars as pl
    source, output, rotation = (tmp_path / name for name in ('input.parquet', 'out.parquet', 'pca.npz'))
    frame = pl.from_numpy(np.random.default_rng(3).normal(size=(12, 4)),
                          schema=[f'phi:{i}' for i in range(4)])
    frame.with_columns(pl.lit('fixed').alias('embedding_space')).write_parquet(source)
    subprocess.run([sys.executable, '-m', 'vdjtools.cli', 'pca', str(source),
                    '--threads', '2', '--components', '2', '--fit', str(rotation),
                    '-o', str(output)], check=True, timeout=30)
    assert pl.read_parquet(output).height == 12
    assert rotation.is_file()


@pytest.mark.skipif(cores._accelerate() is None, reason='requires the macOS 15+ BLAS policy API')
@pytest.mark.parametrize('mode', ['library', 'cli'])
def test_real_accelerate_ceiling_and_worker_policies(mode):
    code = '''
from vdjtools import cores
if MODE == 'library':
    cores.configure_runtime(blas_threads=2)
else:
    cores.configure_cli_runtime(['pca', '--threads=2'])
import numpy as np
import polars as pl
from vdjtools.signature.corpus import TablePCA
from vdjtools.signature.cohort import parallel_rows
apple = cores._accelerate()
assert apple.BLASGetThreading() == 1
assert cores._ACCELERATE_CAP == 2
x = np.arange(80, dtype=float).reshape(20, 4)
table = pl.from_numpy(x, schema=['phi:0','phi:1','phi:2','phi:3'])
table = table.with_columns(pl.lit('fixed').alias('embedding_space'))
def work(i):
    assert apple.BLASGetThreading() == 1
    with cores.numerical_kernels(0):
        assert apple.BLASGetThreading() == 1
    return i
serial = TablePCA.fit(table, n_components=2, threads=1)
with cores.numerical_kernels(2):
    assert apple.BLASGetThreading() == 0
    cores.configure_runtime()  # A downstream import must preserve its caller's phase.
    assert apple.BLASGetThreading() == 0
    parallel = TablePCA.fit(table, n_components=2, threads=2)
    assert parallel_rows([0,1], work, 2) == [0,1]
    assert apple.BLASGetThreading() == 0
assert apple.BLASGetThreading() == 1
np.testing.assert_allclose(parallel.components, serial.components, rtol=0, atol=1e-12)
assert parallel_rows([0,1], work, 2) == [0,1]
assert apple.BLASGetThreading() == 1
if cores.available_cores() > 2:
    try:
        with cores.numerical_kernels(0):
            raise AssertionError('changed ceiling accepted')
    except RuntimeError as exc:
        assert 'fresh process' in str(exc)
'''
    env = {**os.environ, 'VECLIB_MAXIMUM_THREADS': '8', 'POLARS_MAX_THREADS': '1'}
    subprocess.run([sys.executable, '-c', 'MODE=' + repr(mode) + '\n' + code],
                   check=True, timeout=30, env=env)


@pytest.mark.skipif(cores._accelerate() is None, reason='requires the macOS 15+ BLAS policy API')
def test_preinitialized_accelerate_remains_serial_and_fails_parallel_explicitly():
    code = '''
import numpy as np
np.ones((32,32)) @ np.ones((32,32))
from vdjtools import cores
cores.configure_runtime(blas_threads=1)
assert cores._accelerate().BLASGetThreading() == 1
with cores.numerical_kernels(1):
    assert cores._accelerate().BLASGetThreading() == 1
try:
    with cores.numerical_kernels(2):
        raise AssertionError('unverified initialized ceiling accepted')
except RuntimeError as exc:
    assert 'fresh process' in str(exc)
'''
    env = {**os.environ, 'VECLIB_MAXIMUM_THREADS': '1'}
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30, env=env)
