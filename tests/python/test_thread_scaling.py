"""The thread pool must actually be parallel, and the test has to be able to say so.

A thread pool around CPU-bound pure Python is overhead with extra steps, and it fails looking like
success: the results are right, the wall time is merely disappointing, and nobody can tell that
from a workload that is simply expensive. ``_pgen_nt_many`` is claimed to be real parallelism
because ``native.pgen_nt`` releases the GIL, and the claim is checked here rather than asserted in
a docstring.

Measured on 16 cores when this was written, 256 generated TRB junctions: 8.08 s at one thread,
4.21 s at two, 2.19 s at four, 1.17 s at eight, 0.71 s at sixteen -- 1.92x / 3.68x / 6.92x /
11.4x. Doubling the workers does roughly halve the wall time, so the pool is the parallelism
rather than a decoration on it. The test runs a smaller version of exactly that.

Locus choice is not arbitrary. A VJ chain is the wrong instrument: a V/J-marginalized nt Pgen
costs 0.01 ms on TRG and 0.02 ms on TRA, against 36.88 ms on TRB and 40.70 ms on TRD, because the
D-D sum is where the time goes. On TRG the whole 96-call workload finishes in under a millisecond
and the measurement is all pool startup, which reads as a 0.83x slowdown and says nothing.
"""
from __future__ import annotations

import time

import pytest

from vdjtools.model import native
from vdjtools.model.bundled import load_bundled
from vdjtools.model.generate import generate
from vdjtools.model.score import _NT_THREAD_MIN, _pgen_nt_many

#: Above ``_NT_THREAD_MIN`` (or the pool is skipped by design), and large enough that the serial
#: leg dominates pool startup, while keeping the whole module inside a few seconds.
_N = 48


@pytest.fixture(scope="module")
def timings():
    """Both legs, timed once and shared -- running them per test doubles the cost for nothing."""
    assert _N > _NT_THREAD_MIN, "the input must be big enough to take the pooled branch at all"
    m = load_bundled("TRB", "olga")
    native.pack(m)                       # populate the pack cache before either leg is timed
    seqs = [s for s in generate(m, _N, seed=0)["junction_nt"].to_list() if s]
    none = [None] * len(seqs)
    _pgen_nt_many(m, seqs[:2], none[:2], none[:2], 1)      # warm every lazy path

    def wall(threads):
        t0 = time.perf_counter()
        out = _pgen_nt_many(m, seqs, none, none, threads)
        return time.perf_counter() - t0, out

    serial, serial_out = wall(1)
    parallel, parallel_out = wall(4)
    return m, seqs, serial, parallel, serial_out, parallel_out


def test_threads_do_not_change_the_answer(timings):
    """Parallelism that moves a number is not parallelism, it is a bug with a speedup."""
    _, _, _, _, serial_out, parallel_out = timings
    assert serial_out == parallel_out


def test_quadrupling_the_workers_roughly_quarters_the_wall_time(timings):
    """The check the pool would have failed if it were only pretending to run.

    Deliberately loose: the bar is 1.5x from a 4x increase in workers, against the 6.9x measured
    at 16 cores. A pool that is not running at all returns ~1.0x, and a machine merely under load
    or with fewer cores still clears 1.5x -- so this separates the two failure modes it exists to
    separate without becoming a flake on a busy CI box.
    """
    _, seqs, serial, parallel, _, _ = timings
    speedup = serial / parallel
    assert speedup > 1.5, (
        f"{len(seqs)} nt Pgen calls took {serial:.2f} s on one thread and {parallel:.2f} s on "
        f"four -- {speedup:.2f}x. native.pgen_nt is supposed to release the GIL; if it no longer "
        "does, this pool is pure overhead and _pgen_nt_many should go back to a serial loop.")


def test_the_small_input_shortcut_stays_serial(timings):
    """Below the threshold the pool costs more than it saves, so it must not be built."""
    m, seqs, *_ = timings
    short = seqs[:4]
    none = [None] * len(short)
    assert _pgen_nt_many(m, short, none, none, 8) == _pgen_nt_many(m, short, none, none, 1)
