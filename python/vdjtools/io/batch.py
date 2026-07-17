"""Format auto-detection and metadata-driven batch reading.

Sits on top of the single-file readers in :mod:`vdjtools.io.read`: sniffs a file's
format, dispatches to the right reader, and loads a whole cohort described by a
metadata table into one long clonotype frame (or a per-sample dict).
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import polars as pl

from . import convert
from .read import read_airr, read_parquet, read_vdjtools
from .schema import COLUMNS, LOCUS

# Legacy third-party formats (vdjtools.io.convert), keyed by the fmt string used below.
_CONVERTERS = {
    "mixcr": convert.read_mixcr, "migec": convert.read_migec,
    "immunoseq": convert.read_immunoseq, "imgt": convert.read_imgt,
    "vidjil": convert.read_vidjil, "rtcr": convert.read_rtcr,
    "trust4": convert.read_trust4, "arda": convert.read_arda,
    "mitcr": convert.read_mitcr,
}


def _header_columns(path: str | os.PathLike) -> list[str]:
    """Return the (lower-cased) header column names of a TSV without reading rows."""
    cols = pl.read_csv(Path(path), separator="\t", n_rows=0, infer_schema_length=0).columns
    return [c.lower() for c in cols]


def _peek(path: str | os.PathLike, n: int = 8) -> str:
    """Return the first ``n`` characters (gzip-transparent), for JSON detection."""
    with open(path, "rb") as fb:
        gz = fb.read(2) == b"\x1f\x8b"
    opener = gzip.open if gz else open
    with opener(path, "rt", errors="ignore") as f:
        return f.read(n)


def sniff_format(path: str | os.PathLike) -> str:
    """Detect a clonotype file's format from its header.

    Args:
        path: Path to a clonotype table.

    Returns:
        The detected format string, one of: ``"parquet"`` (``.parquet`` / ``.pq``
        extension), ``"vidjil"`` (``.vidjil`` / ``.json`` or a leading ``{``),
        ``"imgt"``, ``"migec"``, ``"mitcr"``, ``"rtcr"``, ``"mixcr"``, ``"immunoseq"``, ``"trust4"``
        (each by its signature header columns), ``"arda"`` (AIRR + arda's ``d2_call``),
        ``"vdjtools"`` (native table / MigMap — ``cdr3aa`` / ``count`` + ``cdr3nt``), or
        ``"airr"`` (AIRR Rearrangement — ``v_call`` / ``junction_aa`` / ``junction_nt`` /
        ``cdr3_aa``).

    Raises:
        ValueError: If no known format signature is recognised.
    """
    if Path(path).suffix.lower() in (".parquet", ".pq"):
        return "parquet"
    if Path(path).suffix.lower() in (".vidjil", ".json") or _peek(path).lstrip()[:1] == "{":
        return "vidjil"  # Vidjil .vidjil JSON
    cols = set(_header_columns(path))
    # Specific third-party formats first (their signature columns don't collide with each other).
    if "v-gene and allele" in cols and "junction" in cols:
        return "imgt"
    if "cdr3 nucleotide sequence" in cols and "v segments" in cols:
        return "migec"
    # MiTCR / tcR (R package): the dot-separated dialect. Checked before `vdjtools`, whose
    # `{count, cdr3nt}` picks don't collide, but keep it with the other third-party formats.
    if {"read.count", "cdr3.nucleotide.sequence"} <= cols:
        return "mitcr"
    if "number of reads" in cols and "junction nucleotide sequence" in cols:
        return "rtcr"
    if cols & {"all v hits", "allvhitswithscore"} and \
            cols & {"n. seq. cdr3", "nseqcdr3", "nseqimputedcdr3"}:
        return "mixcr"
    if "count (templates/reads)" in cols or ({"rearrangement", "amino_acid"} <= cols
                                             and "v_gene" in cols):
        return "immunoseq"
    if "cid_full_length" in cols and "cdr3nt" in cols:
        return "trust4"  # TRUST4 *_report.tsv (its own cid_full_length column)
    # arda AIRR output — standard AIRR names plus arda's tandem-D ``d2_call`` column;
    # match before plain AIRR so it routes to read_arda (nulls arda's ``""`` empty calls).
    if "d2_call" in cols and cols & {"v_call", "junction_aa", "junction", "cdr3_aa"}:
        return "arda"
    # Native vdjtools / migmap (cdr3nt/cdr3aa headers) and AIRR Rearrangement.
    if "cdr3aa" in cols or {"count", "cdr3nt"} <= cols:
        return "vdjtools"
    if cols & {"v_call", "junction_aa", "junction_nt", "cdr3_aa"}:
        return "airr"
    raise ValueError(f"unrecognised clonotype format; header columns: {sorted(cols)}")


def read(path: str | os.PathLike, fmt: str = "auto",
         n_rows: int | None = None) -> pl.DataFrame:
    """Read a clonotype table, auto-detecting the format by default.

    Args:
        path: Path to a native vdjtools, AIRR Rearrangement, or Parquet table, or a
            third-party tool export (MiXcr, MiGec, MiTCR/tcR, immunoSEQ, IMGT/HighV-QUEST,
            Vidjil, RTCR, TRUST4, arda) — see :mod:`vdjtools.io.convert` (``.gz`` ok for the
            text formats).
        fmt: ``"auto"`` (sniff the header / extension), ``"vdjtools"``, ``"airr"``,
            ``"parquet"``, or a legacy format string (``"mixcr"``, ``"migec"``, ``"mitcr"``,
            ``"immunoseq"``, ``"imgt"``, ``"vidjil"``, ``"rtcr"``, ``"trust4"``,
            ``"arda"``).
        n_rows: If given, read at most this many data rows (preview huge files).

    Returns:
        Canonical clonotype frame.

    Raises:
        ValueError: If ``fmt`` is unknown or auto-detection fails.
    """
    if fmt == "auto":
        fmt = sniff_format(path)
    if fmt == "vdjtools":
        return read_vdjtools(path, n_rows=n_rows)
    if fmt == "airr":
        return read_airr(path, n_rows=n_rows)
    if fmt == "parquet":
        return read_parquet(path, n_rows=n_rows)
    if fmt == "vidjil":
        return convert.read_vidjil(path)  # whole-JSON reader (no row cap)
    if fmt in _CONVERTERS:
        return _CONVERTERS[fmt](path, n_rows=n_rows)
    raise ValueError(
        f"fmt must be 'auto', 'vdjtools', 'airr', 'parquet', or a legacy format "
        f"({', '.join(sorted(_CONVERTERS))}); got {fmt!r}"
    )


def read_metadata(path: str | os.PathLike) -> pl.DataFrame:
    """Read a sample metadata TSV.

    The literal string ``"nan"`` (and empty strings) are treated as null. All columns
    are read as strings so metadata joins are stable. A leading ``#`` on the first
    column name (some metadata sheets comment out the header line, e.g.
    ``#file_name\tsample_id\t...``) is stripped.

    Args:
        path: Path to a metadata TSV (one row per sample).

    Returns:
        A ``pl.DataFrame`` of metadata, all-Utf8, with ``"nan"`` → null.
    """
    df = pl.read_csv(Path(path), separator="\t", infer_schema_length=0,
                     null_values=["nan", "NaN", ""])
    if df.columns and df.columns[0].startswith("#"):
        df = df.rename({df.columns[0]: df.columns[0].lstrip("#")})
    return df


def iter_samples(metadata: pl.DataFrame, base_dir: str | os.PathLike,
                 sample_col: str = "sample_name", file_template: str = "{sample}.tsv.gz",
                 fmt: str = "auto", add_metadata: bool = True):
    """Yield ``(sample_id, frame)`` one sample at a time — O(1-sample) RAM.

    The streaming counterpart to :func:`read_samples`. Because it reads and yields
    each sample in turn (never accumulating), a caller can reduce-and-discard —
    per-sample diversity/usage, or sink each sample to a Parquet partition — and so
    process a 100k-sample cohort whose concatenation would not fit in memory. Each
    yielded frame carries the canonical columns + ``locus`` + the reserved
    ``sample_id`` / ``file_name``, and (if ``add_metadata``) that row's metadata.

    Args:
        metadata: Metadata frame (e.g. from :func:`read_metadata`).
        base_dir: Directory holding the per-sample clonotype files.
        sample_col: Metadata column carrying the sample name.
        file_template: Filename template; ``{sample}`` is substituted with the
            sample name.
        fmt: Format for the readers (``"auto"`` / ``"vdjtools"`` / ``"airr"`` /
            ``"parquet"``).
        add_metadata: If ``True``, attach the metadata columns to every clonotype row.

    Yields:
        ``(sample_id, frame)`` tuples in metadata order.
    """
    base = Path(base_dir)
    keep = [*COLUMNS, LOCUS]
    # Exclude the canonical clonotype columns (and the reserved sample tags): a metadata column
    # named e.g. ``locus`` / ``frequency`` / ``duplicate_count`` must NOT overwrite clonotype data.
    reserved = {"sample_id", "file_name", *COLUMNS, LOCUS}
    meta_cols = [c for c in metadata.columns if c not in reserved and c != sample_col]
    for row in metadata.iter_rows(named=True):
        sample = row[sample_col]
        fname = file_template.format(sample=sample)
        clones = read(base / fname, fmt=fmt)
        clones = clones.select([c for c in keep if c in clones.columns])
        clones = clones.with_columns(pl.lit(sample).alias("sample_id"),
                                     pl.lit(fname).alias("file_name"))
        if add_metadata:
            clones = clones.with_columns(pl.lit(row[mc]).alias(mc) for mc in meta_cols)
        yield str(sample), clones


def read_samples(metadata: pl.DataFrame, base_dir: str | os.PathLike,
                 sample_col: str = "sample_name", file_template: str = "{sample}.tsv.gz",
                 fmt: str = "auto", add_metadata: bool = True, as_dict: bool = False):
    """Read a batch of samples described by a metadata table into one frame.

    Eagerly materialises the whole cohort — convenient for tens-to-hundreds of
    samples, but it holds every sample in RAM at once. For large cohorts (thousands+
    of samples) use :func:`iter_samples` to stream, or :func:`vdjtools.io.cohort`
    to persist a hive-partitioned Parquet dataset and scan it lazily.

    For each metadata row the file ``base_dir/file_template.format(sample=<sample>)``
    is read into the canonical schema, tagged with the reserved columns ``sample_id``
    (the ``sample_col`` value) and ``file_name``, and — if ``add_metadata`` — joined
    with that row's metadata columns.

    Args:
        metadata: Metadata frame (e.g. from :func:`read_metadata`).
        base_dir: Directory holding the per-sample clonotype files.
        sample_col: Metadata column carrying the sample name (default
            ``"sample_name"``).
        file_template: Filename template; ``{sample}`` is substituted with the sample
            name (default ``"{sample}.tsv.gz"``).
        fmt: Format for the readers (``"auto"`` / ``"vdjtools"`` / ``"airr"`` /
            ``"parquet"``).
        add_metadata: If ``True``, attach the metadata columns to every clonotype row.
        as_dict: If ``True``, return ``{sample_id: frame}`` instead of one long frame.

    Returns:
        One concatenated long ``pl.DataFrame`` (canonical columns + ``locus`` +
        ``sample_id`` + ``file_name`` + metadata), or a ``dict[str, pl.DataFrame]``
        if ``as_dict``.
    """
    frames = dict(iter_samples(metadata, base_dir, sample_col, file_template,
                               fmt, add_metadata))
    if as_dict:
        return frames
    return pl.concat(list(frames.values()), how="vertical_relaxed")
