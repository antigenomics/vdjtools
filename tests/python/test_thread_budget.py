"""Allocation and import-time pool regressions, using only synthetic work."""
import os
from unittest.mock import patch

import numpy as np
import polars as pl
import pytest

from vdjtools.signature.cohort import parallel_rows


def _pool_probe(i):
    from vdjtools.cores import kernel_threads
    pl.DataFrame({'x': [i]}).select(pl.col('x').sum())
    np.ones((32, 32)) @ np.ones((32, 32))
    return {'i': i, 'polars': pl.thread_pool_size(), 'native': kernel_threads(0)}


def test_thread_workers_use_serial_kernels_without_environment_mutation():
    with patch.dict(os.environ, {'POLARS_MAX_THREADS': '16', 'OPENBLAS_NUM_THREADS': '16'}):
        rows = parallel_rows(range(4), _pool_probe, 2)
        assert [r['i'] for r in rows] == list(range(4))
        assert all(r['polars'] == r['native'] == 1 for r in rows)
        assert os.environ['POLARS_MAX_THREADS'] == '16'
        assert os.environ['OPENBLAS_NUM_THREADS'] == '16'


def test_counts_respect_allocation_and_reject_invalid_requests():
    from vdjtools.cores import kernel_threads
    with patch('vdjtools.cores.available_cores', return_value=2):
        assert kernel_threads(0) == 2
        assert kernel_threads(1) == 1
        assert kernel_threads(2) == 2
        for bad in (-1, 3, True, 1.5):
            with pytest.raises((ValueError, TypeError)):
                kernel_threads(bad)


def test_preimported_blas_is_limited():
    import subprocess
    import sys
    code = '''import numpy as np
np.ones((64,64)) @ np.ones((64,64))
import vdjtools
from threadpoolctl import threadpool_info
assert all(p['num_threads'] == 1 for p in threadpool_info()), threadpool_info()
'''
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '8', 'OMP_NUM_THREADS': '8'}
    subprocess.run([sys.executable, '-c', code], env=env, check=True, timeout=30)


def test_nested_thread_workers_use_serial_kernels():
    from vdjtools.cores import kernel_threads, serial_kernels
    with patch('vdjtools.cores.available_cores', return_value=4):
        with serial_kernels():
            assert kernel_threads(0) == 1
            with pytest.raises(ValueError):
                kernel_threads(2)
        assert kernel_threads(0) == 4


def test_pool_rejects_explicit_overallocation_even_for_small_inputs():
    with patch('vdjtools.cores.available_cores', return_value=2):
        with pytest.raises(ValueError):
            parallel_rows([], _pool_probe, 3)


@pytest.mark.parametrize('v2', [True, False])
def test_nested_cgroup_uses_tightest_ancestor(tmp_path, v2):
    from vdjtools.cores import _cgroup_quota
    parent, child = tmp_path, tmp_path / 'pod'
    child.mkdir()
    if v2:
        (parent / 'cpu.max').write_text('200000 100000')
        (child / 'cpu.max').write_text('max 100000')
    else:
        for path, quota in [(parent, '200000'), (child, '-1')]:
            (path / 'cpu.cfs_quota_us').write_text(quota)
            (path / 'cpu.cfs_period_us').write_text('100000')
    with patch('vdjtools.cores._cgroup_paths', return_value=[(child, v2), (parent, v2)]):
        assert _cgroup_quota() == 2


def test_preinitialized_polars_fails_before_outer_threads_start():
    import subprocess
    import sys
    code = '''import polars as pl
pl.DataFrame({'x':[1]}).select(pl.col('x').sum())
assert pl.thread_pool_size() == 4
import vdjtools
from vdjtools.cores import require_serial_pools
try:
    require_serial_pools()
except RuntimeError as e:
    assert 'fresh process' in str(e)
else:
    raise AssertionError('unbounded immutable pool accepted')
'''
    env = {**os.environ, 'POLARS_MAX_THREADS': '4'}
    subprocess.run([sys.executable, '-c', code], env=env, check=True, timeout=30)


def test_serial_sample_map_accepts_preinitialized_polars():
    import subprocess
    import sys
    from vdjtools.cores import available_cores
    pool = min(4, available_cores())
    code = """import os
import polars as pl
pl.DataFrame({'x':[1]}).select(pl.col('x').sum())
assert pl.thread_pool_size() == int(os.environ['POLARS_MAX_THREADS'])
import vdjtools
from vdjtools.signature.cohort import parallel_rows
assert parallel_rows([1,2], abs, 1) == [1,2]
from vdjtools.cores import available_cores
if available_cores() >= 2:
    try:
        parallel_rows([1,2], abs, 2)
    except RuntimeError:
        pass
    else:
        raise AssertionError('outer pool accepted multithreaded Polars')
"""
    env = {**os.environ, 'POLARS_MAX_THREADS': str(pool)}
    subprocess.run([sys.executable, '-c', code], env=env, check=True, timeout=30)


@pytest.mark.parametrize('v2', [True, False])
def test_cgroup_mount_mapping_includes_visible_ancestors(v2):
    from io import StringIO
    from pathlib import Path
    from vdjtools.cores import _cgroup_paths
    membership = '0::/parent/pod' if v2 else '2:cpu,cpuacct:/parent/pod'
    mount = ('29 24 0:26 /parent /sys/fs/cgroup rw - cgroup2 cgroup rw' if v2 else
             '29 24 0:26 /parent /sys/fs/cgroup/cpu rw - cgroup cgroup rw,cpu,cpuacct')
    def read(path):
        return StringIO(membership if str(path).endswith('/cgroup') else mount)
    with patch('builtins.open', side_effect=read):
        paths = _cgroup_paths()
    base = Path('/sys/fs/cgroup' if v2 else '/sys/fs/cgroup/cpu')
    assert (base / 'pod', v2) in paths
    assert (base, v2) in paths
    assert (base / 'parent' / 'pod', v2) not in paths


def test_preinitialized_polars_cannot_exceed_allocation_in_serial_map():
    import subprocess
    import sys
    code = '''import polars as pl
assert pl.thread_pool_size() == 4
import vdjtools
from vdjtools.signature.cohort import parallel_rows
try:
    parallel_rows([1], abs, 1)
except RuntimeError as exc:
    assert 'Polars' in str(exc) and 'allocation' in str(exc) and 'fresh process' in str(exc)
else:
    raise AssertionError('four-thread Polars pool accepted under a two-CPU allocation')
'''
    env = {**os.environ, 'POLARS_MAX_THREADS': '4', 'SLURM_CPUS_PER_TASK': '2'}
    subprocess.run([sys.executable, '-c', code], env=env, check=True, timeout=30)


def test_serial_map_rechecks_changed_allocation():
    import subprocess
    import sys
    from vdjtools.cores import available_cores
    if available_cores() < 2:
        pytest.skip('requires a two-thread pool before reducing the allocation')
    code = '''import os
import polars as pl
assert pl.thread_pool_size() == 2
import vdjtools
from vdjtools.signature.cohort import parallel_rows
os.environ['SLURM_CPUS_PER_TASK'] = '1'
try:
    parallel_rows([1], abs, 1)
except RuntimeError as exc:
    assert 'allocation' in str(exc)
else:
    raise AssertionError('single-worker map ignored the reduced allocation')
'''
    env = {**os.environ, 'POLARS_MAX_THREADS': '2', 'SLURM_CPUS_PER_TASK': '2'}
    subprocess.run([sys.executable, '-c', code], env=env, check=True, timeout=30)


def test_process_map_can_leave_oversized_parent_polars_pool_idle(tmp_path):
    import subprocess
    import sys
    from vdjtools.cores import available_cores
    if available_cores() < 2:
        pytest.skip('requires two allocated process workers')
    sample = tmp_path / 'sample.tsv'
    sample.write_text('v_call\tj_call\tjunction_aa\tduplicate_count\nTRBV1\tTRBJ1\tCASSLF\t3\n')
    script = tmp_path / 'process_map.py'
    script.write_text('''import sys
import polars as pl
if __name__ == '__main__':
    assert pl.thread_pool_size() == 4
    import vdjtools
    from vdjtools.io import map_samples
    assert map_samples(len, [('a', sys.argv[1]), ('b', sys.argv[1])],
                       workers=2, backend='process') == [('a', 1), ('b', 1)]
''')
    env = {**os.environ, 'POLARS_MAX_THREADS': '4', 'SLURM_CPUS_PER_TASK': '2'}
    subprocess.run([sys.executable, str(script), str(sample)], env=env, check=True, timeout=30)
