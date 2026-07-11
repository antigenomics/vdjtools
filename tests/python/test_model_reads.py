"""Model-reads bootstrap: unique-clonotype dedup + a skip-guarded end-to-end.

The dedup (:func:`vdjtools.model.data.unique_clonotypes`) is the model's clonotype identity — same
V allele, J allele, junction; isotype/CIGAR/error variants collapse — and is unit-tested here on a
synthetic arda clonotype frame (no network / arda needed). The end-to-end fetch+annotate test runs
only when ``huggingface_hub`` + the ``arda`` CLI + the dataset are all available, else skips cleanly.
"""
import shutil

import polars as pl
import pytest

from vdjtools.model import data


def _clones(rows):
    return pl.DataFrame(
        rows,
        schema=["v_call", "j_call", "junction", "junction_aa", "c_call", "locus",
                "d_call", "d2_call", "duplicate_count"],
        orient="row",
    )


def test_unique_clonotypes_collapses_cigar_and_isotype():
    # Same (V, J, junction) but different isotype (IGHM/IGHG) and read support -> ONE clonotype.
    df = _clones([
        ("IGHV3-23*01", "IGHJ4*02", "TGTGCG", "CA", "IGHM", "IGH", "IGHD3-10*01", None, 5),
        ("IGHV3-23*01", "IGHJ4*02", "TGTGCG", "CA", "IGHG1*01", "IGH", "IGHD3-10*01", None, 3),
        ("IGHV3-23*01", "IGHJ4*02", "TGTAAA", "CK", "IGHM", "IGH", None, None, 2),  # diff junction
    ])
    out = data.unique_clonotypes(df)
    assert out.height == 2
    row = out.filter(pl.col("junction") == "TGTGCG").row(0, named=True)
    assert row["count"] == 8  # 5 + 3, isotypes merged
    assert out.filter(pl.col("junction") == "TGTAAA")["count"].item() == 2


def test_unique_clonotypes_alleles_kept_junction_null_dropped():
    df = _clones([
        ("TRBV20-1*01", "TRBJ2-1*01", "TGTGCC", "CA", None, "TRB", None, None, 4),
        ("TRBV20-1*02", "TRBJ2-1*01", "TGTGCC", "CA", None, "TRB", None, None, 1),  # diff V allele
        ("TRBV20-1*01", "TRBJ2-1*01", None, None, None, "TRB", None, None, 9),       # null junction
    ])
    out = data.unique_clonotypes(df)
    assert out.height == 2  # two V alleles kept distinct; null-junction row dropped
    assert set(out["v_call"]) == {"TRBV20-1*01", "TRBV20-1*02"}


def test_unique_clonotypes_naive_igm_filter():
    df = _clones([
        ("IGHV1-2*01", "IGHJ6*02", "TGTGCG", "CA", "IGHM", "IGH", None, None, 3),
        ("IGHV1-2*01", "IGHJ6*02", "TGTAAA", "CK", "IGHG1*01", "IGH", None, None, 7),
    ])
    out = data.unique_clonotypes(df, naive_igm_only=True)
    assert out.height == 1 and out["junction"].item() == "TGTGCG"  # only the IgM clonotype survives


def test_unique_clonotypes_str_count_and_empty_dcalls():
    # Real arda output: duplicate_count can arrive typed as str, and absent D / D2 is "" not null.
    df = pl.DataFrame(
        [
            ("TRBV20-1*01", "TRBJ2-1*01", "TGTGCC", "CA", "", "TRB", "", "", "5"),
            ("TRBV20-1*01", "TRBJ2-1*01", "TGTGCC", "CA", "", "TRB", "TRBD1*01", "", "3"),
        ],
        schema=["v_call", "j_call", "junction", "junction_aa", "c_call", "locus",
                "d_call", "d2_call", "duplicate_count"],
        orient="row",
    )
    assert df["duplicate_count"].dtype == pl.Utf8  # str, as arda writes it
    out = data.unique_clonotypes(df)
    assert out.height == 1
    assert out["count"].item() == 8            # "5" + "3" summed numerically
    assert out["d2_call"].item() is None       # absent second D "" -> null, not a spurious D-D


def test_unique_clonotypes_preserves_d2_for_dd():
    df = _clones([
        ("IGHV3-23*01", "IGHJ4*02", "TGTAAACCC", "KP", "IGHM", "IGH", "IGHD3-10*01", "IGHD6-19*01", 6),
    ])
    out = data.unique_clonotypes(df)
    assert out["d2_call"].item() == "IGHD6-19*01"  # the second D (D-D) is carried through


@pytest.mark.slow
def test_end_to_end_fetch_annotate():
    pytest.importorskip("huggingface_hub")
    if shutil.which("arda") is None:
        pytest.skip("arda CLI not on PATH (needs arda[rnaseq])")
    try:
        fq = data.fetch_fastq("human", "TRB", "nonfunctional")
    except Exception as e:
        pytest.skip(f"dataset unavailable: {e}")
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        clones = data.annotate_reads(fq, out_dir=d, prefix="t", organism="human", cap=2000)
        uniq = data.unique_clonotypes(clones)
        assert uniq.height > 0
        assert {"v_call", "j_call", "junction", "count"} <= set(uniq.columns)
        assert uniq["junction"].null_count() == 0
