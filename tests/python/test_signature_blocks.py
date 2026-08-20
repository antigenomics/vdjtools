"""The statistics blocks, and the depth-honesty property they exist to have.

The load-bearing test here is :class:`TestDepthStability`. A signature is only portable if the
same repertoire sequenced to different depths produces the same numbers, so the columns that
claim depth robustness are checked against a synthetic repertoire subsampled across two orders
of magnitude, and the ones that cannot be depth-robust are checked to *stay* in the ``depth``
block where a downstream model can adjust for them.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.io.schema import C_CALL, COUNT, FREQ, J_CALL, JUNCTION_AA, V_CALL
from vdjtools.preprocess import downsample
from vdjtools.signature import blocks as B

AA = list("ACDEFGHIKLMNPQRSTVWY")


def zipf_repertoire(n: int = 4000, seed: int = 0, locus: str = "TRB") -> pl.DataFrame:
    """A synthetic repertoire with a realistically heavy-tailed clone-size distribution."""
    rng = np.random.default_rng(seed)
    counts = np.ceil(rng.zipf(1.6, n).clip(1, 20_000)).astype(int)
    v = "TRBV20-1" if locus == "TRB" else f"{locus}V1-1"
    j = "TRBJ2-2" if locus == "TRB" else f"{locus}J1"
    return pl.DataFrame({
        V_CALL: [v] * n, J_CALL: [j] * n, C_CALL: [None] * n,
        JUNCTION_AA: ["C" + "".join(rng.choice(AA, rng.integers(8, 18))) + "F" for _ in range(n)],
        COUNT: counts.tolist(), FREQ: (counts / counts.sum()).tolist(),
    })


def stats_of(df: pl.DataFrame) -> dict:
    a = df[COUNT].to_numpy()
    return {"n_reads": float(a.sum()), "richness": float(df.height),
            "f1": float((a == 1).sum()), "f2": float((a == 2).sum()),
            "f3plus": float((a >= 3).sum()),
            "top_clone_fraction": float(a.max() / a.sum())}


class TestSanitise:
    def test_drops_non_standard_junctions(self):
        df = zipf_repertoire(50)
        bad = df.with_columns(
            pl.when(pl.int_range(pl.len()) < 5).then(pl.lit("CASS*RSSYEQYF"))
              .otherwise(pl.col(JUNCTION_AA)).alias(JUNCTION_AA))
        clean, dropped = B.sanitise(bad)
        assert clean.height == 45
        assert dropped > 0

    def test_drops_the_legacy_marker(self):
        """A stop codon or '_' is a real rearrangement that encodes no receptor — dropped."""
        for junk in ("CASS*RSSYEQYF", "CASS_RSSYEQYF"):
            df = zipf_repertoire(10).with_columns(pl.lit(junk).alias(JUNCTION_AA))
            clean, dropped = B.sanitise(df)
            assert clean.height == 0, f"{junk!r} survived the filter"
            assert dropped == pytest.approx(1.0)

    def test_raises_on_ambiguity_codes_and_junk(self):
        """These are a damaged table, not a kind of receptor — an exception, not a filter.

        Changed 2026-08-19. They used to be dropped silently, which hid a broken input and
        inflated the reported non-functional fraction with junk. Measured on 6,047,716 rows of
        real clinical AIRR: zero such characters, so the strict default costs nothing.
        """
        for junk in ("CASSXRSSYEQYF", "cassrssyeqytf", ""):
            df = zipf_repertoire(10).with_columns(pl.lit(junk).alias(JUNCTION_AA))
            with pytest.raises(ValueError, match="unparseable"):
                B.sanitise(df)
            clean, dropped = B.sanitise(df, strict=False)
            assert clean.height == 0 and dropped == pytest.approx(1.0)

    def test_reports_dropped_weight_not_dropped_rows(self):
        """Losing one dominant clone matters more than losing fifty singletons."""
        df = pl.DataFrame({
            V_CALL: ["TRBV20-1"] * 3, J_CALL: ["TRBJ2-2"] * 3, C_CALL: [None] * 3,
            JUNCTION_AA: ["CASS*F", "CASSAF", "CASSGF"],
            COUNT: [900, 50, 50], FREQ: [0.9, 0.05, 0.05],
        })
        _, dropped = B.sanitise(df)
        assert dropped == pytest.approx(0.9)

    def test_empty_frame_is_not_an_error(self):
        clean, dropped = B.sanitise(zipf_repertoire(5).head(0))
        assert clean.height == 0 and dropped == 0.0


class TestWorkFrame:
    def test_weights_close_to_one(self):
        w = B.work_frame(zipf_repertoire(200))
        assert w[FREQ].sum() == pytest.approx(1.0)

    def test_concave_weighting_tames_a_dominant_clone(self):
        """Read weighting would let one huge clone be the entire profile.

        A 10,000-read clone beside two singletons owns 0.9998 of the reads but 0.869 of the
        clone weight, so the two rare clones keep about 650x more say than reads would give
        them — without being promoted to the equal footing that presence weighting implies.
        """
        df = pl.DataFrame({
            V_CALL: ["TRBV20-1"] * 3, J_CALL: ["TRBJ2-2"] * 3, C_CALL: [None] * 3,
            JUNCTION_AA: ["CASSAF", "CASSGF", "CASSTF"],
            COUNT: [10_000, 1, 1], FREQ: [0.0, 0.0, 0.0],
        })
        w = B.work_frame(df)[FREQ].to_numpy()
        raw = 10_000 / 10_002
        assert w[0] < raw, "concave weighting did not reduce the dominant clone's share"
        assert w[1] > 50 * (1 / 10_002), "rare clones gained too little"
        assert w[0] > 1 / 3, "weighting collapsed to presence/absence"

    def test_distinct_weighting_is_uniform(self):
        w = B.work_frame(zipf_repertoire(50), weight="distinct")
        assert w[FREQ].to_numpy().std() == pytest.approx(0.0)


class TestDepthStability:
    """A signature is portable only if depth does not move it."""

    @staticmethod
    @pytest.fixture(scope="class")
    def ladder():
        deep = zipf_repertoire(20_000, seed=7)
        return {n: B.sanitise(downsample(deep, n, by="reads"))[0]
                for n in (3_000, 10_000, 30_000)}

    def test_coverage_standardised_diversity_is_stable(self, ladder):
        vals = [B.div_block(s, cstar=0.20)["1D_c"] for s in ladder.values()]
        assert np.isfinite(vals).all(), "diversity masked out on a sample that should support it"
        assert max(vals) - min(vals) < 0.10, f"1D_c drifted with depth: {vals}"

    def test_clonality_is_stable(self, ladder):
        """Built from standardised Hill numbers; observed evenness drifts by ~2x instead."""
        vals = [B.div_block(s, cstar=0.20)["clonality"] for s in ladder.values()]
        assert max(vals) - min(vals) < 0.15, f"clonality drifted with depth: {vals}"

    def test_junction_length_is_stable(self, ladder):
        vals = [B.len_block(B.work_frame(s))["mean"] for s in ladder.values()]
        assert max(vals) - min(vals) < 0.30, f"length drifted with depth: {vals}"

    def test_depth_itself_is_reported_not_hidden(self, ladder):
        """The columns that cannot be depth-free are the ones naming depth."""
        vals = [B.depth_block(s, stats_of(s))["reads"] for s in ladder.values()]
        assert vals == sorted(vals) and max(vals) - min(vals) > 0.5


class TestDiversityRefusesToExtrapolate:
    def test_unreachable_coverage_returns_holes(self):
        """A hole a model can see beats a confident wrong number."""
        s = zipf_repertoire(200, seed=3)
        assert all(np.isnan(v) for v in B.div_block(s, cstar=0.999).values())

    def test_estimable_flags_the_same_condition(self):
        s = zipf_repertoire(200, seed=3)
        assert not B.estimable(s, 0.999)
        assert B.estimable(s, 0.05)

    def test_degenerate_samples_are_holes(self):
        one = zipf_repertoire(1)
        assert not B.estimable(one, 0.2)
        assert all(np.isnan(v) for v in B.div_block(one, cstar=0.2).values())

    def test_full_tier_adds_columns(self):
        s = zipf_repertoire(3000, seed=5)
        assert set(B.div_block(s, 0.05, tier_full=True)) - set(B.div_block(s, 0.05)) == {
            "0D_chao", "d50"}


class TestQC:
    def test_unknown_v_gene_is_reported(self):
        """The silent germline fallback is the likeliest way a collaborator's vector goes wrong."""
        df = zipf_repertoire(100)
        good = B.qc_block(df, B.work_frame(df), "TRB", 0.0)
        bad = B.qc_block(df, B.work_frame(df.with_columns(pl.lit("TRBV999").alias(V_CALL))),
                         "TRB", 0.0)
        assert bad["v_fallback_frac"] > good["v_fallback_frac"]

    def test_adaptive_nomenclature_is_flagged(self):
        df = zipf_repertoire(100).with_columns(pl.lit("TCRBV09-01").alias(V_CALL))
        assert B.qc_block(df, B.work_frame(df), "TRB", 0.0)["v_fallback_frac"] > -1.0

    def test_unknown_locus_masks_rather_than_raises(self):
        df = zipf_repertoire(10)
        out = B.qc_block(df, B.work_frame(df), "NOTALOCUS", 0.0)
        assert np.isnan(out["v_fallback_frac"])


class TestIsotype:
    def test_composition_uses_the_sample_total_as_denominator(self):
        """Two fifths of IGH reads carry no c_call; that share closes the composition."""
        n = 60
        df = zipf_repertoire(n, locus="IGH").with_columns(
            pl.Series(C_CALL, ["IGHM"] * 20 + ["IGHG1"] * 10 + [None] * (n - 30)))
        out = B.iso_block(B.work_frame(df))
        assert set(out) == {"IgM", "IgD", "IgG", "IgA"}
        assert out["IgM"] > out["IgG"] > out["IgD"]

    def test_full_tier_adds_ige(self):
        df = zipf_repertoire(20, locus="IGH")
        assert "IgE" in B.iso_block(B.work_frame(df), tier_full=True)

    def test_empty_locus_is_a_hole(self):
        assert all(np.isnan(v) for v in B.iso_block(zipf_repertoire(5).head(0)).values())


class TestPairAndShm:
    def test_equal_loci_give_a_zero_ratio(self):
        assert B.pair_block({"TRA": 100.0, "TRB": 100.0})["log_TRA_TRB"] == pytest.approx(0.0)

    def test_absent_locus_is_finite(self):
        assert np.isfinite(B.pair_block({})["log_TRA_TRB"])

    def test_shm_masks_out_without_the_column(self):
        """v_identity is not in the canonical schema, so the standard readers drop it."""
        df = B.work_frame(zipf_repertoire(20, locus="IGH"))
        assert np.isnan(B.shm_block(df)["mean_v_identity"])

    def test_shm_reads_the_column_when_it_survives(self):
        df = B.work_frame(zipf_repertoire(20, locus="IGH")).with_columns(
            pl.lit(0.9).alias("v_identity"))
        assert np.isfinite(B.shm_block(df)["mean_v_identity"])


class TestClonality:
    def test_top_clone_share_is_denominator_aware(self):
        shallow = B.clon_block({"n_reads": 3.0, "richness": 3.0, "f1": 3.0, "f2": 0.0,
                                "f3plus": 0.0, "top_clone_fraction": 1 / 3})
        deep = B.clon_block({"n_reads": 3000.0, "richness": 3000.0, "f1": 3000.0, "f2": 0.0,
                             "f3plus": 0.0, "top_clone_fraction": 1 / 3})
        assert shallow["top"] != pytest.approx(deep["top"])

    def test_ships_two_of_three_composition_parts(self):
        out = B.clon_block(stats_of(zipf_repertoire(500)))
        assert set(out) == {"f1", "f2", "top"}


class TestIsotypeNomenclature:
    """The isotype classes are gene names matched by equality, so the allele has to come off first.

    A frame calling ``IGHG1*01`` would otherwise match nothing at all and come back 100%
    uncalled — which is a perfectly plausible composition (some cohorts really are mostly
    uncalled), so nobody downstream would ever see that the block had silently emptied.
    """

    @staticmethod
    def _frame(call: str, n: int = 60):
        return pl.DataFrame({
            "v_call": ["IGHV1-2"] * n, "j_call": ["IGHJ4"] * n, "c_call": [call] * n,
            "junction_aa": ["CASSLGQAYEQYF"] * n, "duplicate_count": [5] * n,
        })

    def test_an_allele_suffixed_call_reads_the_same_as_a_gene_call(self):
        bare = B.iso_block(B.work_frame(self._frame("IGHG1")))
        allele = B.iso_block(B.work_frame(self._frame("IGHG1*01")))
        assert bare == allele
        assert bare["IgG"] > bare["IgM"], "the IgG-only repertoire did not read as IgG"

    def test_a_null_c_call_column_still_computes(self):
        """An all-null ``c_call`` arrives as a Null-dtype column, which no string op accepts."""
        df = self._frame("IGHM").with_columns(pl.lit(None).alias("c_call"))
        out = B.iso_block(B.work_frame(df))
        assert set(out) == {"IgM", "IgD", "IgG", "IgA"}


class TestPgenJunctionDraw:
    """The frozen ``pgen_q05`` and the ``frac_atypical`` measured against it share one draw."""

    @staticmethod
    def _frame(n: int) -> pl.DataFrame:
        # duplicate_count descending, junction length tracking rank: what a real AIRR file looks
        # like, and enough for a head slice to be distinguishable from a random one.
        return pl.DataFrame({
            "junction_aa": ["C" + "A" * (1 + i % 15) + "F" for i in range(n)],
            "duplicate_count": list(range(n, 0, -1)),
            "v_call": ["TRBV20-1"] * n, "j_call": ["TRBJ2-2"] * n,
        })

    def test_the_draw_is_not_the_head_of_a_sorted_frame(self):
        df = self._frame(5000)
        assert B.pgen_junctions(df, "TRB", 2000) != df["junction_aa"].to_list()[:2000]

    def test_the_draw_is_deterministic_across_calls(self):
        df = self._frame(5000)
        assert B.pgen_junctions(df, "TRB", 2000) == B.pgen_junctions(df, "TRB", 2000)

    def test_a_short_frame_is_taken_whole(self):
        df = self._frame(100)
        assert B.pgen_junctions(df, "TRB", 2000) == df["junction_aa"].to_list()


# ----------------------------------------------------- corrupt vs non-functional (2026-08-19)
#
# sanitise separates two things that used to be one. A stop codon and the legacy '_' marker are a
# real rearrangement that does not encode a receptor -- dropped, and their weight is the number
# vsig:qc:*:nonstd_aa_frac reports. Anything else is a damaged file, and raises.


def _frame(junctions):
    return pl.DataFrame({"junction_aa": junctions,
                         "duplicate_count": [10] * len(junctions)},
                        schema={"junction_aa": pl.Utf8, "duplicate_count": pl.Int64})


@pytest.mark.parametrize("junction", ["CASSXRSSYEQYF", "CASSBRSSYEQYF", "CASSZRSSYEQYF",
                                      "casslrssyeqyf", "CASS1RSSYEQYF", ""])
def test_sanitise_raises_on_corrupt_junctions(junction):
    with pytest.raises(ValueError, match="unparseable"):
        B.sanitise(_frame(["CASSIRSSYEQYF", junction]))


@pytest.mark.parametrize("junction", ["CASS*YEQYF", "CASS_YEQYF"])
def test_sanitise_drops_nonfunctional_without_raising(junction):
    keep, dropped = B.sanitise(_frame(["CASSIRSSYEQYF", junction]))
    assert keep.height == 1
    assert dropped == pytest.approx(0.5)


def test_sanitise_strict_false_restores_the_old_dropping_behaviour():
    keep, dropped = B.sanitise(_frame(["CASSIRSSYEQYF", "CASSXRSSYEQYF"]), strict=False)
    assert keep.height == 1
    assert dropped == pytest.approx(0.5)


def test_sanitise_does_not_filter_on_length():
    # neither predicate has ever had a length bound; a 2-mer and a 62-mer both survive
    keep, dropped = B.sanitise(_frame(["CF", "C" + "A" * 60 + "F"]))
    assert keep.height == 2 and dropped == 0.0


def test_a_null_junction_is_dropped_not_corrupt():
    # null is absence, not a damaged character -- it stays a drop so a partly-annotated table
    # still yields a vector
    keep, dropped = B.sanitise(_frame(["CASSIRSSYEQYF", None]))
    assert keep.height == 1 and dropped == pytest.approx(0.5)
