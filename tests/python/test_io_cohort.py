"""Tests for cohort-scale I/O: iter_samples streaming + hive-partitioned Parquet.

Pins the at-scale ingest path — one sample in RAM at a time, a partitioned Parquet
dataset scanned as a single LazyFrame with sample_id recovered from the path and
metadata joined lazily. These are the invariants a 100k-sample cohort relies on.
"""
import polars as pl

from vdjtools import io as vio
from vdjtools.io import schema as S


def _write_cohort(base):
    """Three tiny AIRR samples on disk; return the metadata frame."""
    base.mkdir(parents=True, exist_ok=True)
    rows = {
        "s1": (["TRBV5-1", "TRBV7-9"], ["CASSL", "CASST"], [10, 5]),
        "s2": (["TRBV5-1", "TRBV20-1"], ["CASSL", "CASSF"], [8, 8]),
        "s3": (["TRBV7-9", "TRBV20-1"], ["CASST", "CASSF"], [3, 20]),
    }
    for sid, (v, cdr3, cnt) in rows.items():
        pl.DataFrame({"v_call": v, "j_call": ["TRBJ2-1"] * len(v),
                      "junction_aa": cdr3, "duplicate_count": cnt}
                     ).write_csv(base / f"{sid}.tsv", separator="\t")
    return pl.DataFrame({"sample_name": ["s1", "s2", "s3"], "age": ["20", "40", "80"]})


def test_iter_samples_streams_one_at_a_time(tmp_path):
    """iter_samples yields (sample_id, frame) lazily, in metadata order, tagged."""
    meta = _write_cohort(tmp_path / "raw")
    out = list(vio.iter_samples(meta, tmp_path / "raw", file_template="{sample}.tsv"))
    assert [sid for sid, _ in out] == ["s1", "s2", "s3"]
    _, f = out[0]
    assert "sample_id" in f.columns and f["sample_id"][0] == "s1"
    assert "age" in f.columns and f["age"][0] == "20"       # metadata attached


def test_ingest_cohort_layout(tmp_path):
    """ingest_cohort writes sample_id=<id>/part.parquet + a single metadata.parquet;
    the clonotype partitions do NOT carry the broadcast metadata (age)."""
    meta = _write_cohort(tmp_path / "raw")
    out = vio.ingest_cohort(meta, tmp_path / "raw", tmp_path / "cohort",
                            file_template="{sample}.tsv")
    parts = sorted(p.name for p in out.glob("sample_id=*"))
    assert parts == ["sample_id=s1", "sample_id=s2", "sample_id=s3"]
    assert (out / "metadata.parquet").exists()
    one = pl.read_parquet(out / "sample_id=s1" / "part.parquet")
    assert "age" not in one.columns and "sample_id" not in one.columns  # not broadcast


def test_scan_cohort_streamed_groupby_and_lazy_metadata(tmp_path):
    meta = _write_cohort(tmp_path / "raw")
    out = vio.ingest_cohort(meta, tmp_path / "raw", tmp_path / "cohort",
                            file_template="{sample}.tsv")
    lf = vio.scan_cohort(out)
    assert isinstance(lf, pl.LazyFrame)
    # sample_id recovered from the path exactly once (no hive/column collision).
    assert lf.collect_schema().names().count("sample_id") == 1
    # all 6 clonotype rows present; metadata.parquet was NOT pulled into the scan.
    assert lf.select(pl.len()).collect().item() == 6
    # metadata joined lazily on sample_id.
    ages = (lf.select(["sample_id", "age"]).unique().collect()
            .sort("sample_id")["age"].to_list())
    assert ages == ["20", "40", "80"]
    # the canonical streamed feature-matrix pass.
    usage = (lf.group_by(["sample_id", S.V_CALL])
             .agg(pl.col(S.COUNT).sum())
             .collect(engine="streaming"))
    assert usage.height == 6                                # 2 genes x 3 samples


def test_scan_cohort_without_metadata_join(tmp_path):
    meta = _write_cohort(tmp_path / "raw")
    out = vio.ingest_cohort(meta, tmp_path / "raw", tmp_path / "cohort",
                            file_template="{sample}.tsv")
    lf = vio.scan_cohort(out, join_metadata=False)
    assert "age" not in lf.collect_schema().names()
    assert "sample_id" in lf.collect_schema().names()
