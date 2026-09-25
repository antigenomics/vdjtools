"""How many cores this process may actually use — which is not ``os.cpu_count()``.

``os.cpu_count()`` reports the machine's cores, not the ones this process is allowed to run on.
Everywhere repertoire analysis actually runs, those differ:

* **SLURM.** Measured on an Aldan-3 ``medium`` node allocated with ``srun -c 8``:
  ``os.cpu_count()`` returns **40** and ``len(os.sched_getaffinity(0))`` returns **8**.
  Sizing a pool off the former starts 40 workers on 8 cores.
* **Containers.** ``docker run --cpuset-cpus`` shows up in the affinity mask, but ``--cpus=N``
  does not -- it is a CFS bandwidth quota, invisible to every API except the cgroup file itself.
  Kubernetes CPU limits are the same mechanism.
* **taskset / numactl**, and any scheduler that pins a job.

Over-allocating is not merely wasteful. Each worker process pays a fresh interpreter and a fresh
load of the frozen artifacts -- on the order of 1.3 s and a few hundred MB for a signature run --
so 40 workers on an 8-core, 32 GB box is several GB of pure overhead competing for cores it does
not have, and it can take the box out of memory for a cohort that would have fitted five times
over.

The order below is deliberate: affinity and cgroup quota are *both* consulted and the smaller
wins, because they constrain independently and a box can carry both.
"""
from __future__ import annotations

import os

__all__ = ["available_cores"]


def _cgroup_quota() -> int | None:
    """Cores implied by a CFS bandwidth quota, or ``None`` if there is no finite one.

    cgroup v2 writes ``"<quota> <period>"`` (or ``"max <period>"``) to ``cpu.max``; v1 splits the
    same pair across two files. ``--cpus=2.5`` is honestly 2 whole workers, so this rounds down,
    and never below 1.
    """
    try:
        with open("/sys/fs/cgroup/cpu.max") as fh:              # cgroup v2
            quota, period = fh.read().split()
        if quota != "max":
            return max(1, int(int(quota) // int(period)))
        return None
    except (OSError, ValueError):
        pass
    try:                                                         # cgroup v1
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as fh:
            quota = int(fh.read())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as fh:
            period = int(fh.read())
        return max(1, quota // period) if quota > 0 and period > 0 else None
    except (OSError, ValueError):
        return None


def available_cores(default: int = 1) -> int:
    """The number of cores this process may actually use, never less than 1.

    Consults, and takes the smallest of: the CPU affinity mask (``os.process_cpu_count()`` on
    3.13+, otherwise ``os.sched_getaffinity``), the cgroup CFS quota, and ``os.cpu_count()``.
    Falls back to ``default`` only if every one of them is unavailable -- which is macOS with a
    broken ``cpu_count``, and not a case anybody is running a cohort on.

    Args:
        default: Returned when nothing at all could be determined.

    Returns:
        A positive core count suitable for sizing a process pool.

    Example:
        >>> from vdjtools.cores import available_cores
        >>> available_cores() >= 1
        True
    """
    counts = []
    proc_count = getattr(os, "process_cpu_count", None)          # 3.13+, affinity-aware
    if proc_count is not None:
        counts.append(proc_count())
    elif hasattr(os, "sched_getaffinity"):                       # Linux, every version
        try:
            counts.append(len(os.sched_getaffinity(0)))
        except OSError:
            pass
    counts.append(_cgroup_quota())
    counts.append(os.cpu_count())
    usable = [c for c in counts if c]
    return max(1, min(usable)) if usable else max(1, default)
