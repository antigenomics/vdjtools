"""Properties the signature transforms exist to guarantee.

These are property tests rather than golden-value tests: what matters is not that a logit
returns a particular number, but that it stays finite at the boundary, stays monotone, and
refuses to pretend a proportion of three reads is as well determined as one of five hundred.
"""
from __future__ import annotations

import numpy as np
import pytest

from vdjtools.signature import transform as T


class TestDenominatorAwareness:
    """The property the shallow-depth regime depends on."""

    def test_logit_of_zero_depends_on_the_denominator(self):
        shallow, deep = T.logit(0.0, 3), T.logit(0.0, 500)
        assert shallow > deep, "0/3 and 0/500 must not map to the same place"
        assert deep < -6

    def test_logit_of_one_depends_on_the_denominator(self):
        assert T.logit(1.0, 3) < T.logit(1.0, 500)

    def test_arcsine_of_zero_depends_on_the_denominator(self):
        assert T.arcsine(0.0, 3) > T.arcsine(0.0, 500)

    def test_shrinkage_is_toward_the_middle(self):
        """A tiny denominator pulls an extreme proportion in, never pushes it out."""
        assert abs(T.logit(0.0, 3)) < abs(T.logit(0.0, 500))

    def test_apply_refuses_a_proportion_without_its_denominator(self):
        for code in ("logit", "arcsine"):
            with pytest.raises(ValueError, match="needs a denominator"):
                T.apply(code, 0.5)


class TestBoundaries:
    def test_finite_at_both_ends(self):
        for m in (0, 1, 3, 10_000):
            assert np.isfinite(T.logit(0.0, m)) and np.isfinite(T.logit(1.0, m))
            assert np.isfinite(T.arcsine(0.0, m)) and np.isfinite(T.arcsine(1.0, m))

    def test_arcsine_is_defined_on_an_empty_sample(self):
        assert np.isfinite(T.arcsine(0.0, 0))

    def test_log_transforms_are_finite_at_zero(self):
        assert T.log10(0.0) == 0.0
        assert T.log1p(0.0) == 0.0
        assert np.isfinite(T.log10(np.array([0.0, 1.0, 1e6]))).all()

    def test_logit_is_symmetric_about_a_half(self):
        assert T.logit(0.5, 10) == pytest.approx(0.0)


class TestMonotonicity:
    @pytest.mark.parametrize("code", ["log10", "log1p"])
    def test_unary_monotone(self, code):
        y = T.apply(code, np.linspace(0.1, 1000, 200))
        assert np.all(np.diff(y) >= 0)

    @pytest.mark.parametrize("code", ["logit", "arcsine"])
    def test_binary_monotone(self, code):
        y = T.apply(code, np.linspace(0, 1, 200), 50)
        assert np.all(np.diff(y) > 0)

    def test_arcsine_is_bounded(self):
        y = T.arcsine(np.linspace(0, 1, 200), 50)
        assert y.min() >= 0 and y.max() <= np.pi / 2 + 1e-9


class TestCLR:
    def test_coordinates_sum_to_zero(self):
        c = T.clr({"a": 4.0, "b": 2.0, "c": 1.0})
        assert sum(c.values()) == pytest.approx(0.0, abs=1e-12)

    def test_ratios_are_preserved(self):
        c = T.clr({"a": 4.0, "b": 2.0, "c": 1.0})
        assert c["a"] - c["b"] == pytest.approx(np.log(2))

    def test_invariant_to_the_total(self):
        a = T.clr({"x": 1.0, "y": 3.0}, m=100)
        b = T.clr({"x": 10.0, "y": 30.0}, m=100)
        assert a["x"] == pytest.approx(b["x"])

    def test_zero_replacement_is_multiplicative(self):
        """A structural zero must not disturb the ratios among the parts that were seen."""
        dense = T.clr({"a": 3.0, "b": 1.0, "c": 5.0})
        sparse = T.clr({"a": 3.0, "b": 1.0, "c": 0.0}, m=4)
        assert dense["a"] - dense["b"] == pytest.approx(np.log(3))
        assert sparse["a"] - sparse["b"] == pytest.approx(np.log(3))

    def test_a_zero_part_is_the_smallest_coordinate(self):
        c = T.clr({"a": 3.0, "b": 1.0, "c": 0.0}, m=4)
        assert c["c"] < c["b"] < c["a"]

    def test_zero_replacement_depends_on_the_denominator(self):
        shallow = T.clr({"a": 3.0, "b": 1.0, "c": 0.0}, m=4)
        deep = T.clr({"a": 300.0, "b": 100.0, "c": 0.0}, m=400)
        assert shallow["c"] > deep["c"], "a zero seen on 400 reads is more certainly a zero"

    def test_an_all_zero_composition_is_flat(self):
        assert all(v == 0.0 for v in T.clr({"a": 0.0, "b": 0.0, "c": 0.0}, m=10).values())

    def test_array_input_returns_an_array(self):
        out = T.clr(np.array([4.0, 2.0, 1.0]))
        assert isinstance(out, np.ndarray) and out.sum() == pytest.approx(0.0, abs=1e-12)

    def test_rejects_degenerate_input(self):
        with pytest.raises(ValueError, match="at least 2 parts"):
            T.clr({"a": 1.0})
        with pytest.raises(ValueError, match="non-negative"):
            T.clr({"a": 1.0, "b": -1.0})

    def test_apply_refuses_clr(self):
        with pytest.raises(ValueError, match="whole composition"):
            T.apply("clr", 0.5)

    def test_subcomposition_is_not_the_full_coordinate(self):
        """Why the emitter must CLR the whole composition and then select coordinates."""
        full = T.clr({"a": 4.0, "b": 2.0, "c": 1.0, "d": 1.0})
        sub = T.clr({"a": 4.0, "b": 2.0})
        assert full["a"] != pytest.approx(sub["a"])


class TestReferenceRescaling:
    def test_clips_symmetrically(self):
        assert T.reference_z(1e6, 0.0, 1.0) == T.DEFAULT_CLIP
        assert T.reference_z(-1e6, 0.0, 1.0) == -T.DEFAULT_CLIP

    def test_a_zero_scale_does_not_divide_by_zero(self):
        assert T.reference_z(3.0, 1.0, 0.0) == 2.0

    def test_centres_and_scales(self):
        assert T.reference_z(5.0, 3.0, 2.0) == 1.0

    def test_robust_statistics_ignore_holes(self):
        x = np.array([[1.0], [2.0], [3.0], [np.nan], [np.inf]])
        loc, scale = T.robust_loc_scale(x)
        assert loc[0] == 2.0
        assert scale[0] > 0

    def test_robust_scale_resists_an_outlier(self):
        clean = np.array([[1.0], [2.0], [3.0], [4.0], [5.0]])
        spiked = np.array([[1.0], [2.0], [3.0], [4.0], [1e9]])
        assert T.robust_loc_scale(spiked)[1][0] == pytest.approx(T.robust_loc_scale(clean)[1][0])

    def test_magnitude_scaling_leaves_the_origin_alone(self):
        """A sample with no deviation must land at the origin, not at minus-the-median."""
        assert np.allclose(T.magnitude_scale(np.zeros(5), 2.0), 0.0)

    def test_magnitude_scaling_keeps_relative_size(self):
        small = T.magnitude_scale(np.array([0.1, 0.0, 0.0]), 2.0)
        big = T.magnitude_scale(np.array([10.0, 0.0, 0.0]), 2.0)
        assert np.linalg.norm(big) > 50 * np.linalg.norm(small)


def test_unknown_transform_raises():
    with pytest.raises(ValueError, match="unknown transform"):
        T.apply("sqrtish", 1.0)
