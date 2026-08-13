"""SHM summary statistics, and the ``keep=`` plumbing they depend on.

The two are tested together because they fail together: ``v_identity`` is not a canonical column,
so without ``keep=`` every reader drops it and every SHM statistic is ``nan`` for a reason that has
nothing to do with the sample.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.io.read import read_airr, read_parquet, read_vdjtools
from vdjtools.stats import shm_spectrum, shm_summary


def _frame(identity, c_call=None, counts=None):
    n = len(identity)
    return pl.DataFrame({
        "junction_aa": [f"CAR{i}W" for i in range(n)],
        "v_call": ["IGHV3-23*01"] * n,
        "j_call": ["IGHJ4*02"] * n,
        "c_call": c_call if c_call is not None else ["IGHG1*01"] * n,
        "duplicate_count": counts or [10] * n,
        "frequency": [1 / n] * n,
        "v_identity": identity,
    })


class TestMissingIdentityIsNotAnUnmutatedRepertoire:
    """The distinction the whole module rests on."""

    def test_absent_column_gives_nan_not_zero(self):
        s = shm_summary(_frame([0.9, 0.8]).drop("v_identity"))
        assert np.isnan(s["mean_identity"])
        assert np.isnan(s["frac_mutated"]), "a sample with no identity read as fully germline"

    def test_all_null_identity_gives_nan(self):
        s = shm_summary(_frame([None, None]))
        assert np.isnan(s["mean_identity"]) and s["n_identity"] == 0

    def test_isotype_still_works_without_identity(self):
        """Class switching needs no identity, so losing one must not lose the other."""
        s = shm_summary(_frame([0.9] * 4, c_call=["IGHM*01", "IGHG1*01", "IGHA1*01", "IGHG3*01"])
                        .drop("v_identity"))
        assert s["frac_switched"] == pytest.approx(0.75)


class TestGerminalCentreSignal:
    def test_mutated_and_germline_separate(self):
        mutated = shm_summary(_frame([0.85, 0.88, 0.91]))
        germline = shm_summary(_frame([1.0, 0.999, 0.995]))
        assert mutated["frac_mutated"] == pytest.approx(1.0)
        assert germline["frac_mutated"] == pytest.approx(0.0)
        assert mutated["mean_shm"] > germline["mean_shm"]

    def test_switch_gap_is_positive_when_switched_cells_are_more_mutated(self):
        """The germinal-centre mark: switched cells carry more SHM than unswitched ones."""
        s = shm_summary(_frame([1.0, 0.99, 0.88, 0.85],
                               c_call=["IGHM*01", "IGHD*01", "IGHG1*01", "IGHA1*01"]))
        assert s["switch_shm_gap"] > 0.1
        assert s["mean_identity_unswitched"] > s["mean_identity_switched"]

    def test_switch_gap_is_nan_when_one_arm_is_absent(self):
        """An IgG-only sample cannot support the comparison, and must not fake it."""
        s = shm_summary(_frame([0.9, 0.85], c_call=["IGHG1*01", "IGHG3*01"]))
        assert np.isnan(s["switch_shm_gap"])
        assert not np.isnan(s["mean_identity_switched"])

    def test_switched_fraction_excludes_uncalled_from_the_denominator(self):
        """`segment_usage` drops null c_call; dividing by the sample total understates switching."""
        s = shm_summary(_frame([0.9] * 4, c_call=["IGHG1*01", "IGHM*01", None, None]))
        assert s["frac_switched"] == pytest.approx(0.5)

    def test_igd_counts_as_unswitched(self):
        s = shm_summary(_frame([0.9] * 2, c_call=["IGHD*01", "IGHG1*01"]))
        assert s["frac_switched"] == pytest.approx(0.5)

    def test_entropy_is_higher_when_several_mutation_levels_coexist(self):
        """Intraclonal diversification: an active GC holds one lineage at many levels at once."""
        spread = shm_summary(_frame([1.0, 0.95, 0.90, 0.85, 0.80, 0.75]))
        tight = shm_summary(_frame([0.90] * 6))
        assert spread["shm_entropy"] > tight["shm_entropy"]


class TestWeighting:
    def test_freq_weighting_follows_clone_size(self):
        """One huge germline clone among mutated singletons: cells vs lineages disagree."""
        df = _frame([1.0, 0.85, 0.85], counts=[1000, 1, 1])
        by_cell = shm_summary(df, weight="reads")
        by_lineage = shm_summary(df, weight="unique")
        assert by_cell["frac_mutated"] < 0.01
        assert by_lineage["frac_mutated"] == pytest.approx(2 / 3)


class TestSpectrum:
    def test_spectrum_is_a_normalised_distribution(self):
        sp = shm_spectrum(_frame([1.0, 0.9, 0.8, 0.7]), bins=10)
        assert sp.height == 10
        assert sp["fraction"].sum() == pytest.approx(1.0)

    def test_spectrum_without_identity_is_nan_not_empty(self):
        sp = shm_spectrum(_frame([0.9]).drop("v_identity"))
        assert sp.height == 20 and sp["fraction"].is_nan().all()


class TestKeepPlumbing:
    """`keep=` has to survive normalize(), which otherwise selects to the canonical set."""

    @pytest.fixture
    def rows(self):
        return {"junction_aa": ["CASSLGQAYEQYF", "CASSLGQAYEQYF", "CARDRGGYW"],
                "junction": ["TGTGCC", "TGTGCC", "TGTGCA"],
                "v_call": ["IGHV3-23*01"] * 2 + ["IGHV1-2*01"],
                "j_call": ["IGHJ4*02"] * 3,
                "duplicate_count": [5, 3, 2],
                "v_identity": [0.90, 0.94, 0.71]}

    def test_airr_keeps_and_averages_over_collapsed_reads(self, tmp_path, rows):
        p = tmp_path / "x.tsv"
        pl.DataFrame(rows).write_csv(p, separator="\t")
        df = read_airr(p, keep=("v_identity",))
        assert df["v_identity"].dtype == pl.Float64, "TSV path left it Utf8"
        got = df.filter(pl.col("junction_aa") == "CASSLGQAYEQYF")["v_identity"][0]
        # mean of the two collapsed reads, not the first: a first-non-null would hand the whole
        # clonotype one read's sequencing error.
        assert got == pytest.approx(0.92)

    def test_parquet_keeps(self, tmp_path, rows):
        p = tmp_path / "x.parquet"
        pl.DataFrame(rows).write_parquet(p)
        assert "v_identity" in read_parquet(p, keep=("v_identity",)).columns

    def test_native_keeps_an_annotation_column(self, tmp_path):
        p = tmp_path / "n.txt"
        pl.DataFrame({"count": [5, 2], "freq": [0.7, 0.3], "cdr3nt": ["TGTGCC", "TGTGCA"],
                      "cdr3aa": ["CASSL", "CARDR"], "v": ["IGHV3-23", "IGHV1-2"],
                      "d": [".", "."], "j": ["IGHJ4", "IGHJ4"],
                      "shm": [0.11, 0.02]}).write_csv(p, separator="\t")
        df = read_vdjtools(p, keep=("shm",))
        assert df["shm"].dtype == pl.Float64

    @pytest.mark.parametrize("reader,suffix,writer", [
        (read_airr, "tsv", "csv"), (read_parquet, "parquet", "parquet")])
    def test_default_is_unchanged_and_absent_names_are_skipped(
            self, tmp_path, rows, reader, suffix, writer):
        p = tmp_path / f"x.{suffix}"
        df0 = pl.DataFrame(rows)
        df0.write_csv(p, separator="\t") if writer == "csv" else df0.write_parquet(p)
        assert "v_identity" not in reader(p).columns
        assert reader(p, keep=("no_such_column",)).height > 0

    def test_read_then_summarise(self, tmp_path, rows):
        """The path the module docstring advertises, end to end."""
        p = tmp_path / "x.tsv"
        pl.DataFrame(rows).write_csv(p, separator="\t")
        s = shm_summary(read_airr(p, keep=("v_identity",)))
        assert not np.isnan(s["mean_identity"]) and s["n_identity"] == 2
