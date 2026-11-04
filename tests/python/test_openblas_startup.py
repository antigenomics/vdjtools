"""OpenBLAS reserves parallel capacity lazily without expanding serial workers."""
import os
import subprocess
import sys

import pytest


def test_lazy_parent_capacity_and_scoped_child_environment():
    code = '''
import os, sys, subprocess
from vdjtools import cores
assert 'numpy' not in sys.modules and 'scipy' not in sys.modules
assert int(os.environ['OPENBLAS_NUM_THREADS']) == cores.kernel_threads(0)
cores.configure_cli_runtime(['pca', '--threads=1'])
assert os.environ['OPENBLAS_NUM_THREADS'] == '1'
if cores.available_cores() >= 2:
    cores.configure_cli_runtime(['pca', '--threads=2'])
    assert os.environ['OPENBLAS_NUM_THREADS'] == '2'
cores.configure_runtime()
import threadpoolctl
previous = dict(os.environ)
with cores.single_threaded_children():
    assert os.environ['OPENBLAS_NUM_THREADS'] == '1'
    subprocess.run([sys.executable, '-c', "import os; from vdjtools import cores; assert os.environ['OPENBLAS_NUM_THREADS'] == '1'; assert cores.kernel_threads(0) == 1; assert cores.numerical_threads(0) == 1"], check=True, timeout=20)
    cores.configure_runtime()
    assert os.environ['OPENBLAS_NUM_THREADS'] == '1'
    with cores.serial_kernels():
        assert cores.numerical_threads(0) == 1
assert dict(os.environ) == previous
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30,
                   env={**os.environ, 'OPENBLAS_NUM_THREADS': '1'})


def test_openblas_capacity_survives_serial_handoffs(monkeypatch):
    from vdjtools import cores
    import threadpoolctl

    monkeypatch.setattr(cores, '_OPENBLAS_CAPS', {})
    pool = {'internal_api': 'openblas', 'filepath': 'test-openblas', 'num_threads': 4}
    monkeypatch.setattr(threadpoolctl, 'threadpool_info', lambda: [pool])
    cores._check_openblas_capacity(4)
    pool['num_threads'] = 1
    cores._check_openblas_capacity(4)
    with pytest.raises(RuntimeError, match='capacity 4, requested 8'):
        cores._check_openblas_capacity(8)


def test_pca_alternates_serial_and_automatic_without_pool_growth_crash():
    code = '''
from vdjtools.cores import numerical_kernels
import numpy as np
from sklearn.decomposition import PCA
for budget in (0, 1, 0, 1):
    for shape, components in (((200, 700), 6), ((45, 80), 2)):
        x = np.random.RandomState(46).normal(size=shape)
        with numerical_kernels(budget):
            assert np.isfinite(PCA(components, random_state=23).fit_transform(x)).all()
'''
    subprocess.run([sys.executable, '-c', code], check=True, timeout=90,
                   env={**os.environ, 'OPENBLAS_NUM_THREADS': '1'})
