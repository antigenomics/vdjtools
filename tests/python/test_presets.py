"""Feature presets: they resolve, they are documented, and they mean what the docs say.

A preset is a promise to a collaborator — "use `transfer` and your model will survive another lab's
samples". These tests pin the parts of that promise that code can check.
"""
from __future__ import annotations

import pytest

from vdjtools.signature import layout as L
from vdjtools.signature import presets as P

ALL = sorted(P.PRESETS)


@pytest.mark.parametrize("name", ALL)
class TestEveryPreset:
    def test_resolves_to_real_columns(self, name):
        cols = P.columns(name)
        assert cols, f"{name} selects nothing"
        contract = set(L.columns("full"))
        assert set(cols) <= contract, f"{name} invented columns"

    def test_column_order_matches_the_layout(self, name):
        """Positional means positional: two users must get the same order, not the same set."""
        cols = P.columns(name)
        order = {c: i for i, c in enumerate(L.columns("full"))}
        assert cols == sorted(cols, key=order.__getitem__)

    def test_is_documented(self, name):
        p = P.get(name)
        for field in ("summary", "features", "how", "use_cases"):
            text = getattr(p, field)
            assert text and len(text) > 20, f"{name}.{field} is not documentation"

    def test_rank_is_one_of_three(self, name):
        assert P.get(name).rank in ("recommended", "specific", "avoid")

    def test_resolution_is_deterministic(self, name):
        assert P.columns(name) == P.columns(name)


class TestTheRelationshipsTheDocsClaim:
    def test_narrower_recommended_presets_are_smaller(self):
        assert P.get("compact").n_columns < P.get("classify").n_columns < P.get("full").n_columns

    def test_transfer_is_a_strict_subset_of_classify(self):
        """`transfer` is documented as `classify` minus the batch-heavy blocks."""
        assert set(P.columns("transfer")) < set(P.columns("classify"))

    def test_the_two_halves_do_not_overlap(self):
        assert not set(P.columns("geometry")) & set(P.columns("statistics"))

    def test_statistics_needs_no_geometry(self):
        """The whole point: usable without mirpy installed."""
        assert all(c.startswith("vsig:") for c in P.columns("statistics"))

    def test_geometry_is_all_rsig(self):
        assert all(c.startswith("rsig:") for c in P.columns("geometry"))

    def test_nuisance_is_disjoint_from_every_feature_preset(self):
        """The floor must not leak into a set it is supposed to be the control for."""
        floor = set(P.columns("nuisance"))
        for name in ALL:
            # `full` is defined as every contract column, so it contains the floor by definition
            # and says so in its notes. Every other preset must not, or a model trained on it
            # could be reading sequencing depth while looking like a biological result.
            if name in ("nuisance", "full"):
                continue
            assert not floor & set(P.columns(name)), f"{name} contains nuisance columns"

    def test_bcell_carries_no_tcr_specific_column(self):
        tcr = [c for c in P.columns("bcell")
               if any(f":{t}:" in c for t in ("TRA", "TRB", "TRG", "TRD"))]
        assert not tcr, tcr[:3]

    def test_bcell_keeps_the_cross_locus_ratios(self):
        """log(IGH/TRB) is a B-cell readout — how B-dominated the sample is — so it stays."""
        assert any(":pair:" in c for c in P.columns("bcell"))

    def test_presets_are_distinct(self):
        seen = {name: tuple(P.columns(name)) for name in ALL}
        assert len(set(seen.values())) == len(seen)


class TestErrorsAndTable:
    def test_unknown_preset_names_the_valid_ones(self):
        with pytest.raises(KeyError) as e:
            P.get("nope")
        assert "transfer" in str(e.value)

    def test_table_covers_every_preset(self):
        t = P.table()
        assert t.height == len(ALL)
        assert set(t["preset"]) == set(ALL)

    def test_table_is_ranked_recommended_first(self):
        assert P.table()["rank"][0] == "recommended"

    def test_at_least_one_of_each_rank(self):
        ranks = {P.get(n).rank for n in ALL}
        assert ranks == {"recommended", "specific", "avoid"}
