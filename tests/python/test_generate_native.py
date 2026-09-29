"""``generate(engine="native")``: the same Bayes net in C++, and what may and may not differ.

The reference sampler is a per-sequence Python loop of ~20 numpy calls. Measured: TRB 25,428
seq/s, IGH with ``productive_only=True`` **4,029 seq/s** -- and pool generation, not featurisation,
is the dominant serial cost of a synthetic corpus build. The native engine is 149-577x that.

**What must be identical**: the columns, their dtypes, the relationship between ``junction_nt``,
``junction_aa`` and ``productive``, and the output for a given ``(seed, row)`` regardless of thread
count. **What is deliberately not**: the random stream. Every shipped artifact drawn from this
function was built with the reference engine, so its ``seed`` stream is frozen and the native
engine gets its own -- the two agree in distribution, not sequence by sequence.
"""
from __future__ import annotations

import numpy as np
import pytest

from vdjtools.model.bundled import load_bundled
from vdjtools.model.generate import generate
from vdjtools.model.reference import translate


@pytest.fixture(scope="module")
def trb():
    return load_bundled("TRB", "olga")


@pytest.fixture(scope="module")
def trg():
    return load_bundled("TRG", "olga")          # a VJ chain: no D, the other code path


def test_the_two_engines_agree_on_columns_and_dtypes(trb):
    """``engine=`` is a speed choice, so it must not change the frame's shape or its schema.

    ``d2_call`` is the one that made this worth pinning: the reference built it from a list of
    ``None`` and polars inferred dtype ``Null`` -- a column whose type depends on its data, which
    breaks every ``concat`` it later meets. Both engines declare ``Utf8`` now.
    """
    ref = generate(trb, 64, seed=1)
    nat = generate(trb, 64, seed=1, engine="native")
    assert ref.columns == nat.columns
    assert ref.schema == nat.schema
    assert nat.height == 64


@pytest.mark.parametrize("locus", ["TRB", "TRG"])
def test_the_translation_and_the_productive_flag_are_self_consistent(locus):
    """``junction_aa`` is translated in the sampler, so it has to equal what Python would get."""
    nat = generate(load_bundled(locus, "olga"), 256, seed=5, engine="native")
    for row in nat.iter_rows(named=True):
        nt, aa = row["junction_nt"], row["junction_aa"]
        assert aa == translate(nt), f"{nt} translated to {aa}"
        assert row["productive"] == (len(nt) > 0 and len(nt) % 3 == 0 and "*" not in aa)


def test_the_draw_does_not_depend_on_the_thread_count(trb):
    """Each slot is seeded from ``(seed, row)``, so the answer is a function of neither.

    A corpus whose contents depend on the builder's core count cannot be compared with one built
    anywhere else -- this is the property that makes that impossible, not a nicety.
    """
    one = generate(trb, 512, seed=9, engine="native", threads=1)
    many = generate(trb, 512, seed=9, engine="native", threads=8)
    auto = generate(trb, 512, seed=9, engine="native", threads=0)
    assert one.equals(many)
    assert one.equals(auto)


def test_the_same_seed_gives_the_same_draw_and_a_different_seed_does_not(trb):
    assert generate(trb, 128, seed=3, engine="native").equals(
        generate(trb, 128, seed=3, engine="native"))
    assert not generate(trb, 128, seed=3, engine="native").equals(
        generate(trb, 128, seed=4, engine="native"))


def test_productive_only_returns_exactly_n_productive_draws(trb):
    nat = generate(trb, 200, seed=2, engine="native", productive_only=True)
    assert nat.height == 200
    assert nat["productive"].all()
    aa = nat["junction_aa"].to_list()
    assert all("*" not in a for a in aa)
    assert all(len(nt) % 3 == 0 for nt in nat["junction_nt"])


def test_an_unknown_engine_raises_rather_than_picking_one(trb):
    """A silently-ignored engine name would be the worst of both: a knob that does nothing."""
    with pytest.raises(ValueError, match="engine"):
        generate(trb, 4, engine="fast")


def test_the_engines_agree_in_distribution_even_though_the_streams_differ(trb):
    """The substantive agreement check: same marginals, different sequences.

    V usage, nucleotide length and the productive fraction are the three things a bug in the
    ancestral order or in a deletion index would move.

    The bars come from a measurement, not from hope. Against a 200,000-draw native sample, TRB V
    usage agrees with the reference at correlation **0.99020 / 0.99406 / 0.99766** and total
    variation **0.0337 / 0.0247 / 0.0156** for a reference draw of 4,000 / 10,000 / 20,000 -- it
    tightens as 1/sqrt(n_ref), which is what two samplers of the *same* distribution do and what a
    mis-indexed deletion would not. Length means agree to 0.036 nt and the productive fraction to
    0.002 throughout. The bars below sit well outside the n=10,000 noise so this is not a flake,
    and well inside anything a real indexing bug would produce.
    """
    n_ref, n_nat = 10_000, 100_000
    ref = generate(trb, n_ref, seed=11)
    nat = generate(trb, n_nat, seed=11, engine="native")

    genes = sorted(set(ref["v_call"].to_list()) | set(nat["v_call"].to_list()))
    pr = np.array([ref["v_call"].eq(g).sum() / n_ref for g in genes])
    pn = np.array([nat["v_call"].eq(g).sum() / n_nat for g in genes])
    assert np.corrcoef(pr, pn)[0, 1] > 0.985, (
        f"V usage correlation {np.corrcoef(pr, pn)[0, 1]:.5f}")
    assert 0.5 * np.abs(pr - pn).sum() < 0.06, f"V usage TV {0.5 * np.abs(pr - pn).sum():.5f}"

    lr = ref["junction_nt"].str.len_chars().to_numpy()
    ln = nat["junction_nt"].str.len_chars().to_numpy()
    assert abs(lr.mean() - ln.mean()) < 0.5, f"nt length {lr.mean():.3f} vs {ln.mean():.3f}"
    assert abs(lr.std() - ln.std()) < 0.5, f"nt length sd {lr.std():.3f} vs {ln.std():.3f}"

    fr, fn = ref["productive"].mean(), nat["productive"].mean()
    assert abs(fr - fn) < 0.02, f"productive fraction {fr:.4f} vs {fn:.4f}"


def test_a_vj_chain_has_no_d_and_still_draws(trg):
    """TRG has no D segment, which is a different branch of the sampler end to end."""
    nat = generate(trg, 128, seed=4, engine="native")
    assert nat["d_call"].is_null().all()
    assert nat["d2_call"].is_null().all()
    assert nat["junction_nt"].str.len_chars().min() > 0


def test_zero_rows_is_an_empty_frame_not_an_error(trb):
    got = generate(trb, 0, seed=1, engine="native")
    assert got.height == 0
    assert got.columns == generate(trb, 1, seed=1, engine="native").columns
