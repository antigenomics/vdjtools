"""``os.cpu_count()`` is the machine's core count, not this process's allowance.

Measured on an Aldan-3 ``medium`` node under ``srun -c 8``: ``os.cpu_count()`` returns **40** and
``len(os.sched_getaffinity(0))`` returns **8**. Sizing a process pool off the former starts 40
interpreters to share 8 cores, each paying a fresh import and a fresh load of the frozen
artifacts. On a 32 GB box that is several GB of overhead for cores that do not exist.

Every pool-sizing site in this package goes through :func:`vdjtools.cores.available_cores`, and
these tests pin the reasons.
"""
from __future__ import annotations

import os
from unittest import mock

import pytest

from vdjtools.cores import _cgroup_quota, available_cores


def test_it_is_positive_and_no_larger_than_the_machine():
    n = available_cores()
    assert n >= 1
    assert n <= (os.cpu_count() or 1)


def test_the_affinity_mask_wins_over_the_machine_count():
    """The SLURM case, which is the one that was actually measured to differ by 5x."""
    with mock.patch.object(os, "cpu_count", return_value=40), \
         mock.patch.object(os, "process_cpu_count", return_value=8, create=True):
        assert available_cores() == 8


def test_a_cgroup_quota_wins_even_when_the_affinity_mask_is_wide():
    """``docker run --cpus=2`` leaves the affinity mask full -- it is a bandwidth quota.

    This is the case no other API reports, and it is how Kubernetes CPU limits work, so a pool
    sized off affinity alone still over-allocates inside a container.
    """
    with mock.patch.object(os, "cpu_count", return_value=40), \
         mock.patch.object(os, "process_cpu_count", return_value=40, create=True), \
         mock.patch("vdjtools.cores._cgroup_quota", return_value=2):
        assert available_cores() == 2


def test_the_smallest_constraint_wins_when_several_apply():
    """Affinity and quota constrain independently, so a box can carry both."""
    with mock.patch.object(os, "cpu_count", return_value=40), \
         mock.patch.object(os, "process_cpu_count", return_value=8, create=True), \
         mock.patch("vdjtools.cores._cgroup_quota", return_value=3):
        assert available_cores() == 3


def test_it_never_returns_zero_even_when_everything_is_unavailable(monkeypatch):
    """A zero would make `range(workers)` empty and the cohort silently produce nothing."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    with mock.patch.object(os, "cpu_count", return_value=None), \
         mock.patch.object(os, "process_cpu_count", return_value=None, create=True), \
         mock.patch("vdjtools.cores._cgroup_quota", return_value=None):
        assert available_cores() == 1
        assert available_cores(default=0) == 1


@pytest.mark.parametrize("text,expected", [
    ("max 100000", None),          # cgroup v2, no limit
    ("200000 100000", 2),          # --cpus=2
    ("250000 100000", 2),          # --cpus=2.5 is honestly two whole workers
    ("50000 100000", 1),           # --cpus=0.5 still has to run something
])
def test_cgroup_v2_quota_parsing(tmp_path, text, expected):
    f = tmp_path / "cpu.max"
    f.write_text(text)
    real_open = open

    def fake_open(path, *a, **k):
        if str(path) == "/sys/fs/cgroup/cpu.max":
            return real_open(f, *a, **k)
        raise OSError("no such file")

    with mock.patch("builtins.open", fake_open):
        assert _cgroup_quota() == expected


def test_the_pool_sizing_sites_use_it():
    """A site that reaches for os.cpu_count() directly re-creates the bug."""
    import pathlib
    root = pathlib.Path(__file__).resolve().parents[2] / "python" / "vdjtools"
    offenders = [
        f"{p.relative_to(root)}:{i}"
        for p in root.rglob("*.py") if p.name != "cores.py"
        for i, line in enumerate(p.read_text().splitlines(), 1)
        if "os.cpu_count()" in line and not line.lstrip().startswith("#")
        and "``" not in line          # docstrings may name it to explain why it is wrong
    ]
    assert not offenders, f"use vdjtools.cores.available_cores() instead: {offenders}"


def test_kernel_auto_budget_respects_allocation_and_pool(monkeypatch):
    from vdjtools.cores import kernel_threads, serial_kernels
    monkeypatch.setattr("vdjtools.cores.available_cores", lambda: 3)
    assert kernel_threads(0) == 3
    with serial_kernels():
        assert kernel_threads(0) == 1
        with pytest.raises(ValueError, match="between 0 and 1"):
            kernel_threads(2)
        with pytest.raises(ValueError, match="threads"):
            kernel_threads(-1)
    assert kernel_threads(0) == 3


def test_memory_guard_uses_remaining_cgroup_bytes(tmp_path, monkeypatch):
    import vdjtools.cores as cores
    (tmp_path/'memory.max').write_text(str(16 * 1024**3))
    (tmp_path/'memory.current').write_text(str(12 * 1024**3))
    monkeypatch.setattr(cores, '_cgroup_paths', lambda *a: [(tmp_path, True)])
    assert cores.available_memory() <= 4 * 1024**3
    cores.require_memory(1024**2, 'small operation')
    with pytest.raises(MemoryError, match='additional workspace'):
        cores.require_memory(5 * 1024**3, 'PCA fit')


@pytest.mark.parametrize("hard_limit", [None, "4"])
def test_import_preserves_caller_openmp_hard_limit(hard_limit):
    import subprocess
    import sys

    env = dict(os.environ)
    if hard_limit is None:
        env.pop("OMP_THREAD_LIMIT", None)
    else:
        env["OMP_THREAD_LIMIT"] = hard_limit
    probe = "import os, vdjtools; print(os.environ.get('OMP_THREAD_LIMIT')); print(os.environ['OMP_NUM_THREADS'])"
    result = subprocess.run([sys.executable, "-c", probe], env=env, timeout=30,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [str(hard_limit), "1"]


def test_auto_cap_preserves_larger_explicit_and_nested_budgets(monkeypatch):
    import vdjtools.cores as cores
    monkeypatch.setattr(cores, 'available_cores', lambda: 32)
    monkeypatch.setattr(cores, 'physical_cores', lambda: 20)
    assert cores.available_cores() == 32
    assert cores.kernel_threads() == cores.numerical_threads() == 20
    assert cores.kernel_threads(32) == cores.numerical_threads(32) == 32
    monkeypatch.setattr(cores._LOCAL, 'numerical', 32, raising=False)
    assert cores.numerical_threads() == 32
    with cores.serial_kernels():
        assert cores.kernel_threads() == cores.numerical_threads() == 1


def test_pool_validation_does_not_mutate_owned_blas_phase(monkeypatch):
    import polars as pl
    import threadpoolctl
    import vdjtools.cores as cores

    monkeypatch.setattr(pl, 'thread_pool_size', lambda: 1)
    with cores.numerical_kernels(1):
        monkeypatch.setattr(threadpoolctl, 'threadpool_limits',
                            lambda *a, **k: pytest.fail('validation changed global BLAS limits'))
        cores.require_serial_pools()
        assert cores.numerical_threads(0) == 1


@pytest.mark.parametrize('timeout', [None, '4', '26'])
def test_import_preserves_openblas_idle_timeout(timeout):
    import subprocess
    import sys

    env = dict(os.environ)
    env.pop('OPENBLAS_THREAD_TIMEOUT', None)
    if timeout is not None:
        env['OPENBLAS_THREAD_TIMEOUT'] = timeout
    result = subprocess.run([sys.executable, '-c',
                             "import vdjtools, os; print(os.environ.get('OPENBLAS_THREAD_TIMEOUT'))"],
                            env=env, capture_output=True, text=True, check=True, timeout=30)
    assert result.stdout.strip() == str(timeout)


def test_physical_topology_counts_smt_once(monkeypatch):
    import vdjtools.cores as cores
    from pathlib import Path
    monkeypatch.setattr(cores, "available_cores", lambda: 40)
    monkeypatch.setattr(cores.os, "sched_getaffinity", lambda pid: set(range(40)), raising=False)
    def topology(path):
        cpu = int(str(path).split("/cpu/cpu")[1].split("/")[0])
        return "0" if path.name == "physical_package_id" else str(cpu % 20)
    monkeypatch.setattr(Path, "read_text", topology)
    assert cores.physical_cores() == 20
    monkeypatch.setattr(cores, "available_cores", lambda: 8)
    assert cores.physical_cores() == 8
