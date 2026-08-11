"""The assembled vsig vector: the contract a collaborator actually joins on.

The blocks are tested elsewhere. What is tested here is the promise the signature makes to
someone who never sees our code — that column *i* means the same thing in their matrix as in
ours, that a locus they did not sequence produces a hole rather than a zero, and that asking for
a narrower tier gives a prefix of the wider one rather than a differently-computed vector.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.signature import layout as L
from vdjtools.signature.assemble import vsig, vsig_cohort

AA = list("ACDEFGHIKLMNPQRSTVWY")


def frame(n=800, v="TRBV20-1", j="TRBJ2-2", c=None, seed=0) -> pl.DataFrame:
    rng = np.random.default_rng(seed)
    return pl.DataFrame(
        {
            "v_call": [v] * n, "j_call": [j] * n, "c_call": [c] * n,
            "junction_aa": ["C" + "".join(rng.choice(AA, 12)) + "F" for _ in range(n)],
            "duplicate_count": np.ceil(rng.zipf(1.5, n).clip(1, 900)).astype(int).tolist(),
        },
        # c_call must be typed even when it is all-null, or concatenating a chain that has
        # constant-gene calls with one that does not fails on the schema rather than on anything real
        schema_overrides={"c_call": pl.Utf8},
    )


@pytest.fixture(scope="module")
def sample():
    return {"TRB": frame(2000, seed=1),
            "IGH": frame(800, "IGHV1-2", "IGHJ4", "IGHM", seed=2)}


class TestContract:
    @pytest.mark.parametrize("tier", L.TIERS)
    def test_columns_match_the_layout_exactly(self, sample, tier):
        assert list(vsig(sample, tier=tier)) == L.columns(tier, "vsig")

    def test_a_narrower_tier_is_a_prefix_not_a_recomputation(self, sample):
        core, std = vsig(sample, tier="core"), vsig(sample, tier="standard")
        for k, v in core.items():
            assert v == std[k] or (np.isnan(v) and np.isnan(std[k])), f"{k} changed with tier"

    def test_every_value_is_a_float(self, sample):
        assert all(isinstance(v, float) for v in vsig(sample).values())

    def test_unknown_tier_raises(self, sample):
        with pytest.raises(ValueError, match="tier must be"):
            vsig(sample, tier="enormous")


class TestHolesNotZeros:
    """A model that reads 'not sequenced' as 'zero' will invent biology that is not there."""

    def test_absent_locus_is_nan_with_a_mask(self, sample):
        v = vsig(sample)
        assert v["vsig:mask:TRA:present"] == 0.0
        for col in ("vsig:div:TRA:1D_c", "vsig:len:TRA:mean", "vsig:depth:TRA:reads"):
            assert np.isnan(v[col]), f"{col} should be a hole"

    def test_present_locus_is_finite(self, sample):
        v = vsig(sample)
        assert v["vsig:mask:TRB:present"] == 1.0
        assert np.isfinite(v["vsig:div:TRB:1D_c"])
        assert np.isfinite(v["vsig:len:TRB:mean"])

    def test_a_trb_only_user_still_gets_a_full_width_vector(self):
        """The portability claim: a collaborator with one chain is not a special case."""
        v = vsig({"TRB": frame(1500, seed=3)})
        assert list(v) == L.columns("standard", "vsig")
        assert v["vsig:mask:TRB:present"] == 1.0
        assert v["vsig:mask:IGH:present"] == 0.0
        assert np.isfinite(v["vsig:div:TRB:1D_c"])

    def test_too_undersampled_for_diversity_is_a_hole_not_a_guess(self):
        """An all-singleton repertoire has Good-Turing coverage 0: nothing was seen twice.

        Every clonotype appearing exactly once is the signature of a repertoire sampled so far
        below saturation that no coverage-standardised estimate is a measurement — reaching any
        target level means extrapolating past everything observed. Refuse rather than guess.
        """
        singletons = frame(400, seed=12).with_columns(pl.lit(1).alias("duplicate_count"))
        v = vsig({"TRB": singletons}, cstar=0.2)
        assert v["vsig:mask:TRB:present"] == 1.0, "the locus was sequenced"
        assert v["vsig:mask:TRB:estimable"] == 0.0, "but not deeply enough to support diversity"
        assert np.isnan(v["vsig:div:TRB:1D_c"])
        assert np.isfinite(v["vsig:depth:TRB:reads"]), "depth itself is still measurable"

    def test_a_tiny_but_saturated_sample_is_still_estimable(self):
        """The guard is about saturation, not about size — and a floor was deliberately not added.

        Three clonotypes over twelve reads has coverage 0.92: shallow, but not undersampled, and
        the target depth stays within twice what was observed. A blanket minimum-clonotype rule
        would discard it, which is a blood assumption rather than a tissue one.
        """
        v = vsig({"TRB": frame(3, seed=4)}, cstar=0.9)
        assert v["vsig:mask:TRB:estimable"] == 1.0
        assert np.isfinite(v["vsig:div:TRB:1D_c"])

    def test_estimable_mask_tracks_the_diversity_columns(self, sample):
        for cstar in (0.05, 0.2, 0.999):
            v = vsig(sample, cstar=cstar)
            assert v["vsig:mask:TRB:estimable"] == float(np.isfinite(v["vsig:div:TRB:1D_c"]))

    def test_an_all_junk_locus_does_not_crash(self):
        junk = frame(50, seed=5).with_columns(pl.lit("CASS*RSSYEQYF").alias("junction_aa"))
        v = vsig({"TRB": junk})
        assert v["vsig:mask:TRB:present"] == 1.0, "the locus was sequenced"
        assert v["vsig:mask:TRB:estimable"] == 0.0, "but nothing usable survived"
        assert np.isnan(v["vsig:div:TRB:1D_c"])


class TestInputForms:
    def test_a_single_frame_with_a_locus_column_works(self):
        both = pl.concat([frame(500, seed=6),
                          frame(500, "IGHV1-2", "IGHJ4", "IGHM", seed=7)])
        from vdjtools.io.schema import add_locus
        v = vsig(add_locus(both))
        assert v["vsig:mask:TRB:present"] == 1.0 and v["vsig:mask:IGH:present"] == 1.0

    def test_dict_and_frame_forms_agree(self):
        d = {"TRB": frame(600, seed=8)}
        from vdjtools.io.schema import add_locus
        a, b = vsig(d), vsig(add_locus(d["TRB"]))
        for k in a:
            assert a[k] == b[k] or (np.isnan(a[k]) and np.isnan(b[k])), k


class TestIGHOnlyColumns:
    def test_isotype_needs_a_constant_call(self):
        with_c = vsig({"IGH": frame(600, "IGHV1-2", "IGHJ4", "IGHM", seed=9)})
        assert with_c["vsig:mask:IGH:c_call"] == 1.0
        assert np.isfinite(with_c["vsig:iso:IGH:IgM"])

    def test_shm_masks_out_because_the_readers_drop_v_identity(self):
        v = vsig({"IGH": frame(600, "IGHV1-2", "IGHJ4", "IGHM", seed=10)})
        assert v["vsig:mask:IGH:shm"] == 0.0
        assert np.isnan(v["vsig:shm:IGH:mean_v_identity"])


class TestCohort:
    def test_cohort_frame_is_in_layout_order(self, sample):
        F = vsig_cohort({"a": sample, "b": {"TRB": frame(900, seed=11)}})
        assert F.columns == ["sample_id", *L.columns("standard", "vsig")]
        assert F.height == 2

    def test_empty_cohort_is_not_an_error(self):
        assert vsig_cohort({}).height == 0

    def test_pair_and_n_loci_are_sample_level(self, sample):
        v = vsig(sample)
        assert v["vsig:qc:-:n_loci_present"] == 2.0
        assert np.isfinite(v["vsig:pair:-:log_IGH_TRB"])
        assert v["vsig:pair:-:log_IGH_TRB"] != 0.0


class TestFullTier:
    """The full tier must compute what it declares, or say plainly that it cannot."""

    def test_composition_and_pgen_are_computed_not_holes(self):
        v = vsig({"TRB": frame(2000, seed=20)}, tier="full")
        for col in ("vsig:aa:TRB:C", "vsig:aa:TRB:G",
                    "vsig:pchem:TRB:all_hydropathy", "vsig:pchem:TRB:center_charge",
                    "vsig:pgen:TRB:mean_log10", "vsig:pgen:TRB:sd_log10"):
            assert np.isfinite(v[col]), f"{col} is declared at full tier but came back a hole"

    def test_residue_composition_is_bounded(self):
        """Arcsine, so every residue column lives in [0, pi/2] whatever the depth."""
        v = vsig({"TRB": frame(2000, seed=21)}, tier="full")
        vals = [v[f"vsig:aa:TRB:{c}"] for c in "ACDEFGHIKLMNPQRSTVWY"]
        assert all(0.0 <= x <= np.pi / 2 for x in vals)

    def test_frac_atypical_needs_a_reference_quantile(self):
        """Being 'atypical' is meaningless without something to be atypical against."""
        assert np.isnan(vsig({"TRB": frame(500, seed=22)}, tier="full")["vsig:pgen:TRB:frac_atypical"])
        with_ref = vsig({"TRB": frame(500, seed=22)}, tier="full", pgen_q05={"TRB": -12.0})
        assert np.isfinite(with_ref["vsig:pgen:TRB:frac_atypical"])

    def test_pgen_is_reproducible_across_calls(self):
        """The subsample is CRC32-seeded, not hash()-seeded, so it survives a new process."""
        a = vsig({"TRB": frame(5000, seed=23)}, tier="full")["vsig:pgen:TRB:mean_log10"]
        b = vsig({"TRB": frame(5000, seed=23)}, tier="full")["vsig:pgen:TRB:mean_log10"]
        assert a == b

    def test_strict_refuses_a_tier_it_cannot_fully_compute(self):
        with pytest.raises(ValueError, match="not computed yet"):
            vsig({"TRB": frame(200, seed=24)}, tier="full", strict=True)

    def test_strict_accepts_a_tier_it_can(self):
        assert vsig({"TRB": frame(200, seed=25)}, tier="core", strict=True)
