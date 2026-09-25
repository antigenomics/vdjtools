"""A repeated amino-acid clonotype key, with no ``junction_nt`` to resolve it, must raise.

Filed from a cohort where it did not. Two exports of the same samples -- one collapsed to the
amino-acid key, one not -- went through the same estimators without complaint; richness differed
by 1.5% overall and 5.0% in IGK, and the affected samples then sat 3-5 robust-SD from the rest of
the cohort on exactly the diversity and clonality columns. It was read as biology until the two
exports were diffed.

The frame cannot answer the question it is being asked. Two rows on one amino-acid key are either
two nucleotide clonotypes encoding the same peptide or one clonotype the export duplicated, and
richness, clonality, Shannon and top-clone fraction all differ between those readings. So the
library must refuse, not count rows and pick one.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.io.schema import (
    assert_resolvable,
    collapse_duplicates,
    duplicate_keys,
    has_nt_resolution,
    resolve_duplicates,
)
from vdjtools.signature.blocks import sanitise
from vdjtools.stats.diversity import diversity_cohort, diversity_stats

#: The exact frame filed in TODO.md: four rows, three distinct keys, one key repeated with
#: different counts, no ``junction_nt``. ``observed_diversity`` used to come back 4, silently.
REPRODUCER = pl.DataFrame({
    "junction_aa": ["CASSLGQGAYEQYF", "CASSLGQGAYEQYF", "CASSPRTGELFF", "CASSQDRGNTIYF"],
    "v_call": ["TRBV5-1*01", "TRBV5-1*01", "TRBV7-9*01", "TRBV4-1*01"],
    "j_call": ["TRBJ2-7*01", "TRBJ2-7*01", "TRBJ2-2*01", "TRBJ1-3*01"],
    "c_call": ["TRBC2"] * 4,
    "duplicate_count": [7, 3, 5, 2],
})


def test_the_filed_reproducer_raises_instead_of_returning_four():
    with pytest.raises(ValueError, match="no junction_nt"):
        diversity_stats(REPRODUCER)


def test_summing_is_available_and_gives_the_other_reading():
    """Both legitimate readings must be reachable, and they must differ -- that is the point."""
    got = diversity_stats(REPRODUCER, on_duplicate="sum")
    assert got["observed_diversity"][0] == 3          # was 4 when the rows were counted as written
    assert got["reads"][0] == 17                      # every read is kept either way


def test_a_frame_carrying_junction_nt_is_never_rejected():
    """With ``junction_nt`` the duplicate is real and resolvable, so there is nothing to refuse."""
    resolved = REPRODUCER.with_columns(
        pl.Series("junction_nt", ["TGTGCC", "TGCGCC", "TGTCCC", "TGTCAA"]))
    assert has_nt_resolution(resolved)
    assert duplicate_keys(resolved).height == 0
    assert_resolvable(resolved)
    assert diversity_stats(resolved)["observed_diversity"][0] == 4


def test_an_all_null_junction_nt_column_does_not_count_as_resolution():
    """A schema-conformant frame always HAS the column; an aa-collapsed export fills it with null.

    Presence alone would have let exactly the filed case through the gate.
    """
    nulled = REPRODUCER.with_columns(pl.lit(None, dtype=pl.Utf8).alias("junction_nt"))
    assert not has_nt_resolution(nulled)
    with pytest.raises(ValueError):
        assert_resolvable(nulled)


def test_the_message_names_the_sample_and_both_resolutions():
    """An error a user cannot act on is barely better than the silence it replaced."""
    with pytest.raises(ValueError) as e:
        assert_resolvable(REPRODUCER, name="S1")
    msg = str(e.value)
    assert "S1" in msg and "CASSLGQGAYEQYF" in msg
    assert "junction_nt" in msg and 'on_duplicate="sum"' in msg


def test_collapsing_sums_counts_and_renormalises_frequency():
    with_freq = REPRODUCER.with_columns((pl.col("duplicate_count") / 17).alias("frequency"))
    got = collapse_duplicates(with_freq)
    assert got.height == 3
    assert got.filter(pl.col("junction_aa") == "CASSLGQGAYEQYF")["duplicate_count"][0] == 10
    assert got["frequency"].sum() == pytest.approx(1.0)
    assert got.columns == with_freq.columns          # column order and set survive the collapse


def test_the_signature_path_is_gated_too():
    """``sanitise`` is where both ``vsig`` and ``mir.signature`` ingest a frame."""
    with pytest.raises(ValueError, match="no junction_nt"):
        sanitise(REPRODUCER)
    assert sanitise(REPRODUCER, on_duplicate="sum")[0].height == 3


def test_the_cohort_path_is_gated_per_sample_and_stays_lazy():
    cohort = REPRODUCER.with_columns(pl.lit("S1").alias("sample_id"))
    with pytest.raises(ValueError, match="S1"):
        diversity_cohort(cohort)
    with pytest.raises(ValueError, match="S1"):
        diversity_cohort(cohort.lazy())               # the lazy path must gate identically
    assert diversity_cohort(cohort, on_duplicate="sum")["observed_diversity"][0] == 3


def test_a_clean_cohort_is_untouched_by_the_check():
    """The gate must cost nothing and change nothing on well-formed data."""
    # rows 0 and 1 ARE the duplicate pair, so the clean frame drops one of them.
    clean = REPRODUCER[[0, 2, 3]].with_columns(pl.lit("S1").alias("sample_id"))
    before = diversity_cohort(clean, on_duplicate="sum")
    assert diversity_cohort(clean).equals(before)
    assert collapse_duplicates(clean).equals(clean)


def test_an_unknown_policy_is_rejected_rather_than_ignored():
    for call in (lambda: resolve_duplicates(REPRODUCER, "drop"),
                 lambda: diversity_cohort(
                     REPRODUCER.with_columns(pl.lit("S1").alias("sample_id")), on_duplicate="drop")):
        with pytest.raises(ValueError, match='"error" or "sum"'):
            call()
