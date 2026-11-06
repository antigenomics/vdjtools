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
from functools import lru_cache

__all__ = ["available_cores", "numerical_threads"]


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


def physical_cores() -> int:
    """Physical cores visible to this process, capped by its logical allocation.

    Unknown topology falls back to the allocation; SMT is never assumed to be 2x.
    """
    from pathlib import Path
    import sys
    import subprocess

    budget = available_cores()
    try:
        if hasattr(os, "sched_getaffinity"):
            cores = set()
            for cpu in os.sched_getaffinity(0):
                topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
                cores.add(((topology / "physical_package_id").read_text().strip(),
                           (topology / "core_id").read_text().strip()))
            return min(budget, len(cores)) if cores else budget
        if sys.platform == "darwin":
            count = int(subprocess.check_output(
                ["sysctl", "-n", "hw.physicalcpu"], timeout=2, text=True))
            return min(budget, count) if count > 0 else budget
    except (OSError, ValueError, subprocess.SubprocessError):
        pass
    return budget


# These runtimes read their limits during import, before a worker initializer runs.
_KERNEL_ENV = ("POLARS_MAX_THREADS", "OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
               "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "RAYON_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMBA_NUM_THREADS")


_LOCAL = threading.local()
# ponytail: native limits are process-global; independent concurrent phases need processes.
_SPAWN_ENV_LOCK = threading.RLock()
_PHASES = threading.Condition(_SPAWN_ENV_LOCK)
_SERIAL_PHASES = 0
_ACCELERATE_CONFIGURED = False
_ACCELERATE_CAP = None
_OPENBLAS_CAPS = {}


def _check_openblas_capacity(threads=1):
    """Initialize BLAS providers before workers, then guard their startup capacity."""
    __import__('scipy.linalg')  # A lazy second BLAS pool must not appear inside workers.
    from threadpoolctl import threadpool_info

    for pool in threadpool_info():
        if pool.get('internal_api') == 'openblas':
            capacity = _OPENBLAS_CAPS.setdefault(pool['filepath'], pool['num_threads'])
            if threads > capacity:
                raise RuntimeError(
                    f'OpenBLAS started with capacity {capacity}, requested {threads}. '
                    'Start a fresh process, import vdjtools and call '
                    'configure_runtime(blas_threads=N) before NumPy/SciPy, '
                    'or use the CLI --threads N.')


@lru_cache(maxsize=1)
def _accelerate():
    """Apple's thread-local BLAS policy API (macOS 15+), absent on other systems."""
    import sys
    if sys.platform != "darwin":
        return None
    import ctypes
    try:
        library = ctypes.CDLL('/System/Library/Frameworks/Accelerate.framework/Accelerate')
        library.BLASSetThreading.argtypes = [ctypes.c_uint]
        library.BLASSetThreading.restype = ctypes.c_int
        library.BLASGetThreading.argtypes = []
        library.BLASGetThreading.restype = ctypes.c_uint
        return library
    except (OSError, AttributeError):
        return None


def _uses_accelerate():
    import sys
    for name in ('numpy', 'scipy'):
        config = getattr(sys.modules.get(name), '__config__', None)
        blas = getattr(config, 'CONFIG', {}).get('Build Dependencies', {}).get('blas', {})
        if 'accelerate' in str(blas.get('name', '')).lower():
            return True
        legacy = getattr(config, 'get_info', None)
        if legacy and any(word in str(legacy('blas_opt_info')).lower()
                          for word in ('accelerate', 'veclib')):
            return True
    return False


@contextmanager
def _accelerate_threads(threads):
    """Restore Apple's thread-local policy; never claim an unenforceable CPU count."""
    library = _accelerate()
    uses = _uses_accelerate()
    if uses and threads > 1 and (library is None or threads != _ACCELERATE_CAP):
        raise RuntimeError(
            'Accelerate cannot change its initialized BLAS CPU ceiling. Start a fresh '
            'process, import vdjtools first and call configure_runtime(blas_threads=N) '
            'before NumPy/SciPy, or use the CLI --threads N; threads=1 is always supported.')
    if library is None or (threads > 1 and not uses):
        yield
        return
    previous = library.BLASGetThreading()
    if library.BLASSetThreading(1 if threads == 1 else 0) != 0:
        raise RuntimeError('Accelerate rejected the requested BLAS threading policy')
    try:
        yield
    finally:
        if library.BLASSetThreading(previous) != 0:
            raise RuntimeError('Accelerate could not restore the BLAS threading policy')


def kernel_threads(threads: int = 0) -> int:
    """Resolve 0/automatic (physical cores), 1/serial, or a positive count within this process's allocation.

    Outer thread-pool workers have a one-core allocation. Explicit excessive counts
    raise rather than silently clamp. Outer worker import pools remain serial; native batches
    receive this resolved integer, never a host-wide automatic setting.
    """
    if isinstance(threads, bool):
        raise TypeError("threads must be an integer, not bool")
    threads = operator.index(threads)
    child_parent = os.environ.get("_VDJTOOLS_SERIAL_CHILD")
    child = child_parent is not None and child_parent != str(os.getpid())
    budget = 1 if getattr(_LOCAL, "serial", False) or child else available_cores()
    if threads < 0:
        raise ValueError("threads must be non-negative (0 means automatic physical cores)")
    if threads > budget:
        raise ValueError(f"threads must be between 0 and {budget} (0 means automatic physical cores)")
    return threads or min(budget, physical_cores())


def numerical_threads(threads: int = 0) -> int:
    """Resolve a BLAS budget; automatic nested calls borrow their caller's phase.

    An explicit positive count overrides the current phase within the allocation.
    Outer sample workers always have a one-core allocation.
    """
    budget = kernel_threads(threads)
    active = getattr(_LOCAL, 'numerical', 0)
    if threads == 0 and active:
        return min(1 if getattr(_LOCAL, 'serial', False) else available_cores(), active)
    return budget


def configure_runtime(*, blas_threads: int | None = None) -> None:
    """Reserve BLAS startup capacity and keep other persistent import pools serial.

    Polars cannot resize an initialized pool. Thread-pool callers detect that case and
    require a fresh interpreter rather than silently oversubscribing it.
    On Apple Accelerate, ``blas_threads`` selects its immutable parallel ceiling
    before NumPy/SciPy imports; None retains the initial allocation/environment cap.
    Later numerical phases support that ceiling or serial execution.
    """
    import sys
    global _ACCELERATE_CONFIGURED, _ACCELERATE_CAP

    if getattr(_LOCAL, 'numerical', False):
        if blas_threads is not None:
            raise RuntimeError('Configure the BLAS ceiling before starting numerical work')
        return
    if getattr(_LOCAL, "serial", False) or os.environ.get("_VDJTOOLS_SERIAL_CHILD") is not None:
        if blas_threads is not None:
            kernel_threads(blas_threads)
        return  # The parent already owns serial limits; never wait for its pool lock.
    with _SPAWN_ENV_LOCK:
        apple = _accelerate()
        preloaded = 'numpy' in sys.modules or 'scipy' in sys.modules
        if not _ACCELERATE_CONFIGURED:
            if not preloaded and apple:
                ceiling = int(os.environ.get('VECLIB_MAXIMUM_THREADS', kernel_threads()))
                if ceiling < 1:
                    raise ValueError('VECLIB_MAXIMUM_THREADS must be positive')
                _ACCELERATE_CAP = min(available_cores(), ceiling)
            _ACCELERATE_CONFIGURED = True
        if blas_threads is not None:
            budget = kernel_threads(blas_threads)
            if _SERIAL_PHASES and budget != 1 and budget != _ACCELERATE_CAP:
                raise RuntimeError('Configure the BLAS ceiling before starting sample workers')
            if preloaded and budget != 1 and _uses_accelerate() and budget != _ACCELERATE_CAP:
                raise RuntimeError('Configure the Accelerate BLAS ceiling before importing '
                                   'NumPy/SciPy in a fresh process')
            if not preloaded:
                _ACCELERATE_CAP = budget if apple else None
        # Keep serial defaults without imposing a hard cap on later caller-owned contexts.
        os.environ.update(dict.fromkeys((key for key in _KERNEL_ENV
                                        if key != "OMP_THREAD_LIMIT"), "1"))
        if apple:
            os.environ['VECLIB_MAXIMUM_THREADS'] = str(_ACCELERATE_CAP or 1)
            if apple.BLASSetThreading(1) != 0:
                raise RuntimeError('Accelerate rejected serial BLAS initialization')
        # OpenBLAS 0.3.28 can crash when growing a pool initialized at one thread.
        # Reserve the authorized capacity before lazy imports; phases set active limits.
        os.environ["OPENBLAS_NUM_THREADS"] = str(kernel_threads(blas_threads or 0))
        os.environ['OMP_WAIT_POLICY'] = 'PASSIVE'
        if "numpy" in sys.modules or "scipy" in sys.modules:
            from threadpoolctl import threadpool_limits

            _check_openblas_capacity()
            threadpool_limits(limits=1)


def configure_cli_runtime(argv=None):
    """Select the CLI BLAS ceiling before command modules import NumPy/SciPy."""
    import argparse
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument('-t', '--threads', type=int)
    args, _ = parser.parse_known_args(argv)
    configure_runtime(blas_threads=args.threads)


@contextmanager
def single_threaded_children():
    """Set import-time limits before spawning workers and restore the parent environment."""
    # ponytail: serialize generation pools; their import limits share the process environment.
    with _SPAWN_ENV_LOCK:
        previous = {key: os.environ.get(key) for key in (*_KERNEL_ENV, "_VDJTOOLS_SERIAL_CHILD")}
        os.environ.update(dict.fromkeys(_KERNEL_ENV, "1"))
        os.environ["_VDJTOOLS_SERIAL_CHILD"] = str(os.getpid())
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
    if getattr(_LOCAL, "serial", False):
        yield
        return
    from threadpoolctl import threadpool_limits

    previous = getattr(_LOCAL, "serial", False)
    _LOCAL.serial = True
    try:
        # The parent owns process-wide limits; each worker also owns thread-local pools.
        with threadpool_limits(limits=1), _accelerate_threads(1):
            yield
    finally:
        _LOCAL.serial = previous


@contextmanager
def serial_runtime():
    """Allow independent serial-kernel maps together, excluding numerical fits."""
    global _SERIAL_PHASES
    if getattr(_LOCAL, "serial", False):
        with _accelerate_threads(1):
            yield
        return
    if getattr(_LOCAL, 'numerical', False):
        with numerical_kernels(1):
            yield
        return
    from threadpoolctl import threadpool_limits

    with _PHASES:
        _check_openblas_capacity()
        threadpool_limits(limits=1)
        _SERIAL_PHASES += 1
    try:
        yield
    finally:
        with _PHASES:
            _SERIAL_PHASES -= 1
            _PHASES.notify_all()


@contextmanager
def numerical_kernels(threads: int = 1):
    """Temporarily budget BLAS after outer workers have joined; restore caller limits.

    ``0`` borrows an existing numerical phase or selects allocated physical cores; ``1`` is
    serial. Process-wide changes serialize
    with sample and generation pools; calls inside a sample worker stay serial.
    Accelerate supports its startup ceiling or serial; other counts and unknown
    preinitialized ceilings raise, rather than silently running serially.
    """
    from threadpoolctl import threadpool_limits

    budget = numerical_threads(threads)
    if getattr(_LOCAL, "serial", False):
        # The parent owns the serial pools; changing their process-global limits here
        # races sibling workers and waiting for its lock would deadlock the pool.
        yield budget
        return
    with _PHASES:
        _PHASES.wait_for(lambda: _SERIAL_PHASES == 0)
        _check_openblas_capacity(budget)
        with threadpool_limits(limits=1), threadpool_limits(limits=budget, user_api="blas"), \
                _accelerate_threads(budget):
            previous = getattr(_LOCAL, 'numerical', False)
            _LOCAL.numerical = budget
            try:
                yield budget
            finally:
                _LOCAL.numerical = previous


def require_serial_pools(*, outer: bool = True) -> None:
    """Validate immutable pool allocation without changing an owning runtime phase."""
    import polars as pl

    pool = pl.thread_pool_size()
    budget = available_cores()
    if pool > budget:
        raise RuntimeError(f"Polars has {pool} threads but the CPU allocation is {budget}; "
                           "start a fresh process and import vdjtools before Polars")
    if not outer:
        return
    if pool > 1:
        raise RuntimeError("Polars was initialized with multiple threads before vdjtools; "
                           "start a fresh process and import vdjtools before Polars, "
                           "or use workers=1")
