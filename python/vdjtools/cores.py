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
import operator
import threading
from contextlib import contextmanager

__all__ = ["available_cores"]


def _cgroup_paths(controller_name="cpu"):
    """Visible controller directories, including ancestors that can impose tighter quotas."""
    from pathlib import Path

    paths = [(Path("/sys/fs/cgroup"), True), (Path("/sys/fs/cgroup") / controller_name, False),
             (Path("/sys/fs/cgroup/cpu,cpuacct"), False)]
    try:
        with open("/proc/self/cgroup") as fh:
            memberships = [line.strip().split(":", 2) for line in fh]
        with open("/proc/self/mountinfo") as fh:
            mounts = [line.strip().split(" - ", 1) for line in fh]
    except OSError:
        return paths
    for left, right in mounts:
        fields, controller = left.split(), right.split()
        v2 = controller[0] == "cgroup2"
        if not v2 and (controller[0] != "cgroup" or controller_name not in controller[2].split(",")):
            continue
        root, mount = (Path(x.replace(r"\040", " ")) for x in fields[3:5])
        for _, ctrls, name in memberships:
            if (v2 and ctrls) or (not v2 and controller_name not in ctrls.split(",")):
                continue
            group = Path(name)
            if group.is_relative_to(root):
                current = mount / group.relative_to(root)
                while current.is_relative_to(mount):
                    paths.append((current, v2))
                    if current == mount:
                        break
                    current = current.parent
    return list(dict.fromkeys(paths))


def _cgroup_quota() -> int | None:
    """Smallest finite v1/v2 quota in the current cgroup and its visible ancestors.

    Fractional cores round down, with one worker for allocations below one CPU.
    CPU requests are scheduling weights, not hard allocations, and cannot set this ceiling.
    """
    limits = []
    for directory, v2 in _cgroup_paths():
        try:
            if v2:
                with open(directory / "cpu.max") as fh:
                    quota, period = fh.read().split()
                if quota == "max":
                    continue
            else:
                with open(directory / "cpu.cfs_quota_us") as fh:
                    quota = fh.read().strip()
                with open(directory / "cpu.cfs_period_us") as fh:
                    period = fh.read().strip()
            quota, period = int(quota), int(period)
            if quota > 0 and period > 0:
                limits.append(max(1, quota // period))
        except (OSError, ValueError):
            continue
    return min(limits) if limits else None


def available_memory() -> int | None:
    """Remaining bytes under the tightest visible memory cgroup and host limit."""
    from pathlib import Path

    limits = []
    for directory, v2 in _cgroup_paths("memory"):
        limit_name, used_name = (("memory.max", "memory.current") if v2 else
                                 ("memory.limit_in_bytes", "memory.usage_in_bytes"))
        try:
            limit = int((directory / limit_name).read_text().strip())
            used = int((directory / used_name).read_text().strip())
            if 0 < limit < 2**60:
                limits.append(max(0, limit - used))
        except (OSError, ValueError):
            continue
    try:
        values = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
        limits.append(int(values['MemAvailable'].split()[0]) * 1024)
    except (OSError, ValueError, KeyError):
        pass
    return min(limits) if limits else None


def require_memory(nbytes: int, operation: str) -> None:
    """Reject a known large allocation before OOM, reserving 20% of remaining memory.

    This guards calculated numerical workspaces, not arbitrary caller-owned objects.
    Unknown host limits leave the estimate to the caller; they do not prove safety.
    """
    remaining = available_memory()
    if remaining is not None and nbytes > remaining * 0.8:
        raise MemoryError(f"{operation} needs approximately {nbytes / 1024**3:.2f} GiB "
                          f"of additional workspace; {remaining / 1024**3:.2f} GiB remains. "
                          "Use fewer workers, a saved fitted basis, or a larger memory allocation.")


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
    try:
        counts.append(int(os.environ["SLURM_CPUS_PER_TASK"]))
    except (KeyError, ValueError):
        pass
    counts.append(_cgroup_quota())
    counts.append(os.cpu_count())
    usable = [c for c in counts if c]
    return max(1, min(usable)) if usable else max(1, default)


# These runtimes read their limits during import, before a worker initializer runs.
_KERNEL_ENV = ("POLARS_MAX_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
               "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "RAYON_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")


_LOCAL = threading.local()
_SPAWN_ENV_LOCK = threading.RLock()


def kernel_threads(threads: int = 0) -> int:
    """Resolve 0/all, 1/serial, or a positive count within this process's allocation.

    Outer thread-pool workers have a one-core allocation. Explicit excessive counts
    raise rather than silently clamp. Numerical import pools remain serial; native batches
    receive this resolved integer, never a host-wide automatic setting.
    """
    if isinstance(threads, bool):
        raise TypeError("threads must be an integer, not bool")
    threads = operator.index(threads)
    budget = 1 if getattr(_LOCAL, "serial", False) else available_cores()
    if threads < 0:
        raise ValueError("threads must be non-negative (0 means allocated cores)")
    if threads > budget:
        raise ValueError(f"threads must be between 0 and {budget} (0 means allocated cores)")
    return threads or budget


def configure_runtime() -> None:
    """Keep persistent numerical pools serial before imports; limit initialized BLAS too.

    Polars cannot resize an initialized pool. Thread-pool callers detect that case and
    require a fresh interpreter rather than silently oversubscribing it.
    """
    import sys

    os.environ.update(dict.fromkeys(_KERNEL_ENV, "1"))
    if "numpy" in sys.modules or "scipy" in sys.modules:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=1)


@contextmanager
def single_threaded_children():
    """Set import-time limits before spawning workers and restore the parent environment."""
    # ponytail: serialize generation pools; their import limits share the process environment.
    with _SPAWN_ENV_LOCK:
        previous = {key: os.environ.get(key) for key in _KERNEL_ENV}
        os.environ.update(dict.fromkeys(_KERNEL_ENV, "1"))
        try:
            yield
        finally:
            for key, value in previous.items():
                if value is None:
                    os.environ.pop(key, None)
                else:
                    os.environ[key] = value


@contextmanager
def serial_kernels():
    """Bound native work within an outer thread worker, restoring thread-local state."""
    previous = getattr(_LOCAL, "serial", False)
    _LOCAL.serial = True
    try:
        yield
    finally:
        _LOCAL.serial = previous


@contextmanager
def numerical_kernels(threads: int = 1):
    """Temporarily budget BLAS after outer workers have joined; restore caller limits.

    ``0`` selects the allocation, ``1`` is serial. This process-wide context must
    not overlap an outer worker pool; native worker calls remain serial.
    """
    from threadpoolctl import threadpool_limits

    budget = kernel_threads(threads)
    with threadpool_limits(limits=budget, user_api="blas"):
        yield budget


def require_serial_pools() -> None:
    """Reject an immutable multi-thread Polars pool before launching outer thread workers."""
    import polars as pl
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1)
    if pl.thread_pool_size() > 1:
        raise RuntimeError("Polars was initialized with multiple threads before vdjtools; "
                           "start a fresh process and import vdjtools before Polars, "
                           "or use workers=1")
