"""Cohort-scale I/O: a hive-partitioned Parquet dataset scanned as one LazyFrame.

The single-file readers and :func:`vdjtools.io.read_samples` materialise a whole
cohort in RAM — fine for tens-to-hundreds of samples, fatal for 100k. This module
provides the at-scale path:

1. :func:`ingest_cohort` streams every sample through the readers **one at a time**
   (never holding two at once) and writes each to its own Parquet partition
   ``<out_dir>/sample_id=<id>/part.parquet``, with the sample metadata sheet stored
   once as ``<out_dir>/metadata.parquet`` (not broadcast into every clonotype row).
2. :func:`scan_cohort` opens the whole cohort as a single :class:`polars.LazyFrame`
   — ``sample_id`` recovered from the partition path, metadata joined lazily — so a
   cohort far larger than RAM is analysed with one streaming ``group_by``:

   >>> lf = scan_cohort("cohort/")                       # doctest: +SKIP
   >>> usage = (lf.group_by(["sample_id", "v_call"])     # doctest: +SKIP
   ...            .agg(pl.col("duplicate_count").sum())
   ...            .collect(engine="streaming"))

Every analysis in the package is already ``group_by(...).agg(...)``, so the cohort
feature matrix is one streamed pass; only the final (bounded) sample×feature pivot
is materialised.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl

from .batch import iter_samples
from .schema import LOCUS, SCHEMA

#: Reserved partition-key column recovered from the ``sample_id=<id>`` path.
SAMPLE_ID = "sample_id"


def ingest_cohort(metadata: pl.DataFrame, base_dir: str | os.PathLike,
                  out_dir: str | os.PathLike, sample_col: str = "sample_name",
                  file_template: str = "{sample}.tsv.gz", fmt: str = "auto") -> Path:
    """Convert a sample cohort into a hive-partitioned Parquet dataset (streaming).

    Reads each sample in turn via :func:`vdjtools.io.iter_samples` (O(1-sample) RAM)
    and writes the canonical clonotype frame to ``out_dir/sample_id=<id>/part.parquet``.
    The metadata sheet is written once to ``out_dir/metadata.parquet`` with its sample
    column renamed to ``sample_id`` — it is *not* broadcast into the clonotype rows.
    Run this once; thereafter analyse the cohort with :func:`scan_cohort`.

    Args:
        metadata: Metadata frame (e.g. from :func:`vdjtools.io.read_metadata`).
        base_dir: Directory holding the per-sample clonotype files.
        out_dir: Destination directory for the partitioned dataset (created if absent).
        sample_col: Metadata column carrying the sample name.
        file_template: Filename template; ``{sample}`` is substituted per sample.
        fmt: Reader format (``"auto"`` / ``"vdjtools"`` / ``"airr"`` / ``"parquet"``).

    Returns:
        The ``out_dir`` :class:`~pathlib.Path` (feed straight into :func:`scan_cohort`).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    # ponytail: sample ids must be filesystem/hive-safe (no '/', '='); real sample
    #           names are; sanitize here if that ever stops holding.
    for sid, frame in iter_samples(metadata, base_dir, sample_col, file_template,
                                   fmt, add_metadata=False):
        part = out / f"{SAMPLE_ID}={sid}"
        part.mkdir(parents=True, exist_ok=True)
        # Drop sample_id/file_name — the partition path already encodes sample_id.
        frame.drop(SAMPLE_ID, "file_name").write_parquet(part / "part.parquet")
    meta = (metadata.rename({sample_col: SAMPLE_ID})
            if sample_col != SAMPLE_ID else metadata)
    meta.write_parquet(out / "metadata.parquet")
    return out


def scan_cohort(out_dir: str | os.PathLike, *, join_metadata: bool = True) -> pl.LazyFrame:
    """Scan a hive-partitioned cohort (from :func:`ingest_cohort`) as one LazyFrame.

    Nothing is read until ``.collect()``; ``group_by`` / ``filter`` / column selection
    push down into the Parquet scan, so a cohort far larger than memory is reduced
    with ``.collect(engine="streaming")``. ``sample_id`` is recovered from the
    partition path; the metadata columns are joined lazily on ``sample_id``.

    Args:
        out_dir: The partitioned dataset directory written by :func:`ingest_cohort`.
        join_metadata: If ``True`` (default) and ``metadata.parquet`` is present,
            left-join the per-sample metadata onto the clonotype rows lazily.

    Returns:
        A :class:`polars.LazyFrame` over the whole cohort (canonical clonotype
        columns + ``sample_id`` + metadata).
    """
    out = Path(out_dir)
    # Empty cohort (no partitions written): return an empty LazyFrame with the
    # canonical clonotype schema rather than letting the scan raise on a no-match glob.
    if not any(out.glob(f"{SAMPLE_ID}=*/*.parquet")):
        return pl.DataFrame(schema={**SCHEMA, LOCUS: pl.Utf8, SAMPLE_ID: pl.Utf8}).lazy()
    # Explicit sample_id=*/*.parquet glob so the top-level metadata.parquet is never
    # pulled into the clonotype scan; hive parsing recovers sample_id. Force it to
    # Utf8 (hive_schema) so numeric sample ids ('1','2') don't infer as Int64 and
    # break the join with the all-Utf8 metadata — and partition pruning is preserved.
    clones = pl.scan_parquet(str(out / f"{SAMPLE_ID}=*" / "*.parquet"),
                             hive_partitioning=True, hive_schema={SAMPLE_ID: pl.String})
    meta_path = out / "metadata.parquet"
    if join_metadata and meta_path.exists():
        meta = pl.scan_parquet(meta_path).with_columns(pl.col(SAMPLE_ID).cast(pl.Utf8))
        clones = clones.join(meta, on=SAMPLE_ID, how="left")
    return clones
