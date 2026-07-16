"""match="fuzzy" is a SEARCH (candidate keeps identity, gains incidence), not clustering."""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.biomarker import association
from vdjtools.io import schema as S


def _cohort():
    """20 subjects. CASSPAAAF is carried EXACTLY by 4 cases; its 1mm neighbour CASSPAAAY by
    6 more cases and 1 control. Exact -> weak; fuzzy -> the union (10 cases / 1 control)."""
    rows = []
    for i in range(4):                       # exact carriers of the candidate (cases)
        rows.append(("CASSPAAAF", "TRBV1", f"case{i}"))
    for i in range(4, 10):                   # carriers of the 1mm neighbour only (cases)
        rows.append(("CASSPAAAY", "TRBV1", f"case{i}"))
    rows.append(("CASSPAAAY", "TRBV1", "ctl0"))          # one control has the neighbour
    for i in range(10):                      # filler so every subject exists
        rows.append(("CAAAAAAAW", "TRBV9", f"case{i}"))
    for i in range(10):
        rows.append(("CAAAAAAAW", "TRBV9", f"ctl{i}"))
    aa, v, sid = zip(*rows)
    return pl.DataFrame({S.JUNCTION_AA: list(aa), S.V_CALL: list(v), "sample_id": list(sid),
                         S.COUNT: [1] * len(rows)})


def _design():
    return pl.DataFrame({"sample_id": [f"case{i}" for i in range(10)]
                                      + [f"ctl{i}" for i in range(10)],
                         "_pos": [True] * 10 + [False] * 10})


@pytest.mark.parametrize("match,expected", [("exact", 4), ("fuzzy", 10)])
def test_fuzzy_search_gains_incidence_and_keeps_identity(match, expected):
    res = association(_cohort(), _design(), test="fisher",
                                  key=(S.JUNCTION_AA, S.V_CALL), match=match,
                                  min_incidence=1, alternative="greater")
    hit = res.filter(pl.col(S.JUNCTION_AA) == "CASSPAAAF")
    # identity survives: the candidate is still its own row, not folded into a group
    assert hit.height == 1, f"{match}: candidate lost its identity"
    assert int(hit["n_pos_present"][0]) == expected


def test_fuzzy_universe_is_not_restricted_to_candidates():
    """The neighbour is NOT a candidate; fuzzy must still find it (regression: the candidate
    semi-join used to shrink the search universe, silently under-counting incidence)."""
    cand = pl.DataFrame({S.JUNCTION_AA: ["CASSPAAAF"], S.V_CALL: ["TRBV1"]})
    res = association(_cohort(), _design(), test="fisher",
                                  key=(S.JUNCTION_AA, S.V_CALL), match="fuzzy",
                                  min_incidence=1, candidates=cand, alternative="greater")
    hit = res.filter(pl.col(S.JUNCTION_AA) == "CASSPAAAF")
    assert int(hit["n_pos_present"][0]) == 10, "fuzzy universe was clipped to the candidate set"


def test_fuzzy_v_must_match_exactly():
    """V is part of the key -> a 1mm neighbour on a different V is NOT the same feature."""
    c = _cohort().with_columns(
        pl.when(pl.col("sample_id").str.starts_with("case") & (pl.col(S.JUNCTION_AA) == "CASSPAAAY"))
          .then(pl.lit("TRBV7")).otherwise(pl.col(S.V_CALL)).alias(S.V_CALL))
    res = association(c, _design(), test="fisher", key=(S.JUNCTION_AA, S.V_CALL),
                                  match="fuzzy", min_incidence=1, alternative="greater")
    hit = res.filter((pl.col(S.JUNCTION_AA) == "CASSPAAAF") & (pl.col(S.V_CALL) == "TRBV1"))
    assert int(hit["n_pos_present"][0]) == 4, "V-mismatched neighbours were wrongly counted"
