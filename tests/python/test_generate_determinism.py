"""`generate(seed=)` must be reproducible ACROSS PROCESSES, not just within one.

⛔ It was not. Same seed, same wheel, a new interpreter — a different draw. Within a single process
it looked perfect, which is why an in-process test would never have caught it.

Root cause: `collapse_alleles` (which `load_bundled(..., collapse=True)` runs by default) built its
tables with unordered polars `group_by().agg()`. `group_by` is a multithreaded hash aggregation, so
the collapsed table's ROW ORDER varied per process; `_cum` then assigned the same cumulative
interval to a different allele, and the same `rng.random()` selected a different one.

⚠ Not hash randomisation — `PYTHONHASHSEED=0` did not help, which is what ruled that out and
pointed at the aggregation instead. Same class as the nondeterminism recorded against arda's
`correct` stage.
"""
from __future__ import annotations

import subprocess
import sys

import pytest

from vdjtools.model import load_bundled
from vdjtools.model.generate import generate

_SNIPPET = """
from vdjtools.model import load_bundled
from vdjtools.model.generate import generate
m = load_bundled({locus!r}, 'olga')
print(generate(m, 3, seed=11).row(0, named=True)['junction_nt'])
"""


@pytest.mark.parametrize("locus", ["TRA", "TRB"])
def test_same_seed_gives_the_same_draw_in_a_FRESH_PROCESS(locus):
    """The test that the bug required: a *subprocess*, because in-process it always agreed."""
    outs = set()
    for _ in range(3):
        r = subprocess.run([sys.executable, "-c", _SNIPPET.format(locus=locus)],
                           capture_output=True, text=True, timeout=300)
        assert r.returncode == 0, r.stderr[-2000:]
        outs.add(r.stdout.strip())
    assert len(outs) == 1, (
        f"{locus}: seed=11 produced {len(outs)} different draws across processes: {outs}"
    )


@pytest.mark.parametrize("locus", ["TRA", "TRB"])
def test_the_collapsed_table_order_is_stable_in_a_FRESH_PROCESS(locus):
    """The layer underneath, so a failure says WHICH stage drifted rather than just 'output moved'."""
    snippet = (
        "from vdjtools.model import load_bundled\n"
        f"m = load_bundled({locus!r}, 'olga', collapse=True)\n"
        "print(m.tables['v_choice']['v_allele'].to_list())\n"
    )
    outs = set()
    for _ in range(3):
        r = subprocess.run([sys.executable, "-c", snippet], capture_output=True, text=True,
                           timeout=300)
        assert r.returncode == 0, r.stderr[-2000:]
        outs.add(r.stdout.strip())
    assert len(outs) == 1, f"{locus}: collapse_alleles emitted {len(outs)} different row orders"


def test_within_one_process_it_was_always_fine():
    """Pins why the bug survived: repeated calls in one interpreter agreed all along."""
    m = load_bundled("TRB", "olga")
    a = generate(m, 5, seed=3)["junction_nt"].to_list()
    b = generate(m, 5, seed=3)["junction_nt"].to_list()
    assert a == b
