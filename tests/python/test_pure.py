"""Pure-Python unit tests — no OLGA/arda/models needed, so they always run.

Covers the load-bearing helpers (genetic code, palindrome cut segments) and the graph/schema
validators, whose error paths the oracle-gated integration tests never exercise.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjtools.model.events import Event, EventKind, validate_graph
from vdjtools.model.reference import cut_segment, cut_segment_d, reverse_complement, translate
from vdjtools.model.schema import Manifest, validate_tables


def test_reverse_complement():
    assert reverse_complement("GTGT") == "ACAC"
    assert reverse_complement("CTAA") == "TTAG"
    assert reverse_complement("") == ""


def test_translate():
    assert translate("TGT") == "C" and translate("TGC") == "C"          # Cys
    assert translate("ATG") == "M"                                       # Met
    assert translate("TAA") == translate("TAG") == translate("TGA") == "*"  # stops
    assert translate("TGTGCC") == "CA"                                   # two codons
    assert translate("TGTG") == "C"                                      # trailing partial dropped


def test_cut_segment():
    # V: append RC of the last max_pal nt at the 3' end.
    assert cut_segment("TGTGCAGAGGTGT", "V", 4) == "TGTGCAGAGGTGT" + reverse_complement("GTGT")
    # J: prepend RC of the first max_pal nt at the 5' end.
    assert cut_segment("CTAACTAT", "J", 4) == reverse_complement("CTAA") + "CTAACTAT"
    # short segment: only min(len, max_pal) palindromic nt.
    assert cut_segment("TG", "V", 4) == "TG" + reverse_complement("TG")
    with pytest.raises(ValueError, match="V.*or.*J"):
        cut_segment("ACGT", "D", 4)


def test_cut_segment_d():
    d = cut_segment_d("GGGTAC", 4, 4)
    assert d == reverse_complement("GGGT") + "GGGTAC" + reverse_complement("GTAC")


def test_event_graph_valid_and_errors():
    validate_graph({"a": Event("a", EventKind.GENE_CHOICE),
                    "b": Event("b", EventKind.GENE_CHOICE, ("a",))})  # no raise
    with pytest.raises(ValueError, match="cycle"):
        validate_graph({"a": Event("a", EventKind.GENE_CHOICE, ("b",)),
                        "b": Event("b", EventKind.GENE_CHOICE, ("a",))})
    with pytest.raises(ValueError, match="unknown parent"):
        validate_graph({"b": Event("b", EventKind.GENE_CHOICE, ("a",))})


def test_manifest_chain_type_and_roundtrip():
    events = {"v_choice": Event("v_choice", EventKind.GENE_CHOICE),
              "j_choice": Event("j_choice", EventKind.GENE_CHOICE, ("v_choice",))}
    m = Manifest(locus="TRA", organism="human", chain_type="VJ", events=events)
    r = Manifest.from_json(m.to_json())
    assert r.chain_type == "VJ" and r.events == events
    with pytest.raises(ValueError, match="VDJ.*VJ"):
        Manifest(locus="TRA", organism="human", chain_type="XX", events=events)


def _mini_model():
    events = {
        "v_choice": Event("v_choice", EventKind.GENE_CHOICE),
        "v_3_del": Event("v_3_del", EventKind.DELETION, ("v_choice",)),
    }
    manifest = Manifest(locus="TRA", organism="human", chain_type="VJ", events=events)
    tables = {
        "v_choice": pl.DataFrame({"v_allele": ["a", "b"], "p": [0.5, 0.5]}),
        # gene 'a' has a real deletion distribution; gene 'b' is all-zero (undefined conditional).
        "v_3_del": pl.DataFrame({
            "v_allele": ["a", "a", "b", "b"],
            "ndel": [0, 1, 0, 1],
            "p": [0.7, 0.3, 0.0, 0.0],
        }),
    }
    return manifest, tables


def test_validate_tables_ok_and_zero_group():
    manifest, tables = _mini_model()
    validate_tables(manifest, tables)  # no raise: 'a' sums to 1, 'b' sums to 0 (allowed)


def test_validate_tables_errors():
    manifest, tables = _mini_model()
    with pytest.raises(ValueError, match="missing the marginal"):
        validate_tables(manifest, {"v_choice": tables["v_choice"]})
    with pytest.raises(ValueError, match="columns"):
        bad = tables["v_choice"].rename({"v_allele": "wrong"})
        validate_tables(manifest, {**tables, "v_choice": bad})
    with pytest.raises(ValueError, match="sum to neither"):
        bad = tables["v_choice"].with_columns(p=pl.Series("p", [0.5, 0.7]))  # sums to 1.2
        validate_tables(manifest, {**tables, "v_choice": bad})
