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


def test_it_never_returns_zero_even_when_everything_is_unavailable():
    """A zero would make `range(workers)` empty and the cohort silently produce nothing."""
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
