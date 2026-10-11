"""The vectorized feature groups, and the two faults found while vectorizing them.

Both faults were silent. The first shipped wrong numbers: ``pchem_group`` asked
``physchem_profile`` to group by locus, and a locus is derived **per row** from the gene name --
so a TRD repertoire, which legitimately carries TRAV V genes because TRA and TRD share their V
segments, split into two groups and the tidy-frame loop let the last one win. Every TRD physchem
column was the weighted mean over one subset of the repertoire rather than over the repertoire.

The second is a performance fault that reads as correctness here: the build must be bit-identical
at every ``n_jobs``, because a corpus whose values depend on the builder's core count cannot be
compared with one built anywhere else.
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.features.physchem import DEFAULT_PROPERTIES, load_property_table
from vdjtools.signature import corpus as C
from vdjtools.signature import features as FE
from vdjtools.signature.layout import AMINO_ACIDS

#: A TRD frame as a real one looks: TRDV calls and TRAV calls side by side. The junctions differ in
#: composition on purpose, so a group that is dropped changes the answer.
TRD_MIXED = pl.DataFrame({
    "junction_aa": ["CAKWWWWWWF", "CAKWWWWWWF", "CAKGGGGGGF", "CAKGGGGGGF"],
    "v_call": ["TRDV1*01", "TRDV2*01", "TRAV13-1*01", "TRAV38-2*01"],
    "j_call": ["TRDJ1*01"] * 4,
    "duplicate_count": [10, 10, 10, 10],
})


def _expected_pchem_all(df: pl.DataFrame, prop: str) -> float:
    """Weighted mean over EVERY clonotype of its own residue-mean of ``prop``."""
    tbl = {r["amino_acid"]: r[prop] for r in load_property_table().iter_rows(named=True)}
    w = df["frequency"].to_numpy()
    per = np.array([np.mean([tbl[c] for c in s]) for s in df["junction_aa"]])
    return float((per * w).sum() / w.sum())


def test_a_trd_frame_carrying_trav_calls_is_not_split_by_derived_locus():
    """The whole repertoire, not whichever derived-locus group happened to be written last."""
    wf = FE.work_frame(FE.sanitise(TRD_MIXED)[0])
    got = FE.pchem_group(wf)
    for prop in ("hydropathy", "volume", "kf1"):
        assert got[f"all_{prop}"] == pytest.approx(_expected_pchem_all(wf, prop), rel=1e-12), prop


def test_dropping_either_half_of_that_frame_gives_a_different_answer():
    """Guards the test above: if both halves agreed, it could not detect the bug it exists for."""
    wf = FE.work_frame(FE.sanitise(TRD_MIXED)[0])
    whole = FE.pchem_group(wf)["all_hydropathy"]
    for rows in ([0, 1], [2, 3]):
        half = FE.pchem_group(FE.work_frame(TRD_MIXED[rows]))["all_hydropathy"]
        assert abs(half - whole) > 0.1, "the two V-gene groups must differ for this to be a test"


def test_the_center_region_skips_a_junction_that_has_no_centre():
    """``_region_expr`` returns null below five residues; the vectorized path must skip the same."""
    short = pl.DataFrame({"junction_aa": ["CAF", "CAKWWWWWWF"], "v_call": ["TRBV2*01"] * 2,
                          "j_call": ["TRBJ1-1*01"] * 2, "duplicate_count": [10, 10]})
    wf = FE.work_frame(FE.sanitise(short)[0])
    got = FE.pchem_group(wf)
    only_long = FE.pchem_group(FE.work_frame(short[[1]]))
    for prop in DEFAULT_PROPERTIES:
        assert got[f"center_{prop}"] == pytest.approx(only_long[f"center_{prop}"], rel=1e-12)
    # and the 'all' region keeps both, so the two regions genuinely disagree here
    assert got["all_volume"] != pytest.approx(got["center_volume"], rel=1e-6)


def test_the_composition_groups_count_what_a_hand_count_counts():
    """``aa`` and ``kmer`` are weighted counts; the bincount path must reproduce the obvious loop."""
    df = pl.DataFrame({"junction_aa": ["CAAF", "CAAAF"], "v_call": ["TRBV2*01"] * 2,
                       "j_call": ["TRBJ1-1*01"] * 2, "duplicate_count": [3, 1]})
    wf = FE.work_frame(FE.sanitise(df)[0])
    w = wf["frequency"].to_numpy()
    seqs = wf["junction_aa"].to_list()

    acc = dict.fromkeys(AMINO_ACIDS, 0.0)
    pairs: dict[str, float] = {}
    n_tok = 0
    for s, wi in zip(seqs, w):
        for ch in s:
            acc[ch] += wi
        for i in range(len(s) - 1):
            pairs[s[i:i + 2]] = pairs.get(s[i:i + 2], 0.0) + wi
            n_tok += 1
    m = sum(len(s) for s in seqs)
    tot = sum(acc.values())

    aa = FE.aa_group(wf)
    from vdjtools.signature.transform import arcsine
    for a, v in acc.items():
        assert aa[a] == pytest.approx(float(arcsine(v / tot, m)), rel=1e-12)
    km = FE.kmer_group(wf)
    ptot = sum(pairs.values())
    for p, v in pairs.items():
        assert km[p] == pytest.approx(float(arcsine(v / ptot, max(float(n_tok), 1.0))), rel=1e-12)
    # a k-mer must not straddle two clonotypes: "FC" would be the join of "CAAF" and "CAAAF"
    assert km["FC"] == pytest.approx(min(km.values()), rel=0), "a k-mer crossed a junction boundary"


def test_an_unrecognised_gene_lands_in_the_residual_not_in_a_gene():
    """``usage`` and ``spec`` close their composition with a residual cell; a stranger goes there."""
    df = pl.DataFrame({"junction_aa": ["CAKWWF", "CAKGGF"],
                       "v_call": ["TRBV2*01", "TCRBV99-9*01"],   # the second is Adaptive-style junk
                       "j_call": ["TRBJ1-1*01"] * 2, "duplicate_count": [10, 10]})
    wf = FE.work_frame(FE.sanitise(df)[0])
    vocab = FE.gene_vocab("TRB")
    usage = FE.usage_group(wf, "v_call", vocab["vus"])
    assert np.isfinite(list(usage.values())).all()
    # the recognised gene carries weight; every other vocabulary gene is a structural zero and so
    # must share one value, strictly below it
    assert usage["TRBV2"] == max(usage.values())
    others = {v for k, v in usage.items() if k != "TRBV2"}
    assert len(others) == 1 and max(others) < usage["TRBV2"]
    qc = FE.qc_channel(df, wf, "TRB", 0.0)
    assert qc["v_fallback_frac"] == pytest.approx(0.5, abs=1e-9)   # one of two, by equal weight


def test_the_corpus_build_is_bit_identical_across_worker_counts(tmp_path):
    """A corpus whose values depend on the builder's core count is not comparable to anything.

    Two workers, not sixteen: the point is that the pooled path and the in-process path agree, and
    a 2-worker pool exercises every join the 16-worker one does.
    """
    kw = dict(loci=("TRG",), n_samples=8, size=60, seed=5, n_components=3)
    serial, _ = C.synthesize("memory", n_jobs=1, **kw)
    pooled, _ = C.synthesize("memory", n_jobs=2, **kw)
    a, b = serial.save(tmp_path / "a"), pooled.save(tmp_path / "b")
    assert a.read_bytes() == b.read_bytes(), "the corpus depends on how many workers built it"


def test_the_pool_chunk_count_is_a_constant_and_not_the_core_count():
    """``POOL_CHUNKS`` decides the pools' contents, so it must not be read off the machine."""
    import inspect

    src = inspect.getsource(C.build_pools)
    assert "POOL_CHUNKS" in src
    assert "available_cores" not in src, "the chunk count must not depend on this machine"
    assert isinstance(C.POOL_CHUNKS, int) and C.POOL_CHUNKS >= 1
