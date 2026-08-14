"""The signature column contract: names, nesting, index subsets, and no silent drift.

These are contract tests, not numeric ones — nothing here computes a feature. They exist
because a collaborator joins on these column names, so a change to any of them is a breaking
change that must be a deliberate version bump rather than an accident.
"""
from __future__ import annotations

import pytest

from vdjtools.signature import layout


def test_tiers_are_nested():
    for sig in ("vsig", "rsig", None):
        prev: set[str] = set()
        for tier in layout.TIERS:
            cols = set(layout.columns(tier, sig))
            assert prev <= cols, f"{sig} {tier} drops a column from the previous tier"
            prev = cols


def test_tiers_are_index_subsets_of_full():
    """A narrower tier must be sliceable out of the full matrix, not recomputed."""
    for sig in ("vsig", "rsig", None):
        full = layout.columns("full", sig)
        for tier in layout.TIERS:
            assert [full[i] for i in layout.index(tier, sig)] == layout.columns(tier, sig)


def test_no_duplicate_columns():
    full = layout.columns("full")
    assert len(full) == len(set(full))


def test_signatures_do_not_collide():
    """depth and div exist under both signatures; the sig prefix is what keeps them apart."""
    v = set(layout.columns("full", "vsig"))
    r = set(layout.columns("full", "rsig"))
    assert not v & r
    assert {"depth", "div"} <= {b.name for b in layout.registry("vsig")}
    assert {"depth", "div"} <= {b.name for b in layout.registry("rsig")}


def test_every_column_parses():
    for c in layout.columns("full"):
        sig, block, locus, feature = layout.parse(c)
        assert sig in ("vsig", "rsig")
        assert locus in (*layout.LOCI, layout.NO_LOCUS)
        assert block and feature


def test_parse_rejects_malformed():
    with pytest.raises(ValueError, match="malformed"):
        layout.parse("vsig:div:TRB")


def test_unknown_tier_and_block_raise():
    with pytest.raises(ValueError, match="tier must be"):
        layout.columns("tiny")
    with pytest.raises(ValueError, match="unknown block"):
        layout.columns("core", blocks_=("nosuchblock",))


def test_describe_matches_columns():
    for tier in layout.TIERS:
        d = layout.describe(tier)
        assert d["column"].to_list() == layout.columns(tier)


def test_mask_is_exempt_and_contrast_is_magnitude_scaled():
    """The two rescaling exemptions are contract, not implementation detail.

    A mask is already 0/1 and must not be z-scored; the contrast block carries its meaning in
    its magnitude, so per-column centring would erase the very deficiency it encodes.
    """
    by = {(b.sig, b.name): b for b in layout.registry()}
    assert by[("vsig", "mask")].exempt
    assert by[("rsig", "contrast")].magnitude
    assert not by[("rsig", "contrast")].exempt


def test_attributable_is_declared_only_for_geometry():
    """Only a block with a clonotype pre-image may claim attributability."""
    attributable = {b.name for b in layout.registry() if b.attributable}
    assert attributable == {"contrast", "phiv", "phij", "phic"}
    assert all(b.sig == "rsig" for b in layout.registry() if b.attributable)


def test_register_rejects_a_duplicate_feature():
    dup = layout.Block("vsig", "div", layout.feats("core", "log10", "1D_c"))
    with pytest.raises(ValueError, match="already declares"):
        layout.register(dup)


def test_every_feature_declares_a_known_transform():
    d = layout.describe("full")
    assert set(d["transform"].unique()) <= set(layout.TRANSFORMS)


def test_transform_is_per_feature_not_per_block():
    """A heterogeneous block is the normal case, not an exception."""
    clon = next(b for b in layout.registry("vsig") if b.name == "clon")
    assert clon.transform("f1") == "clr"
    assert clon.transform("top") == "logit"


def test_block_rejects_an_unknown_transform():
    with pytest.raises(ValueError, match="unknown transform"):
        layout.feats("core", "sqrtish", "x")
