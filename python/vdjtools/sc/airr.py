"""The single-cell interchange layer: a flat AIRR Rearrangement table.

scirpy, dandelion and scRepertoire all read the **same** thing — one AIRR Rearrangement
row per contig, carrying ``sequence_id`` and ``cell_id``. None of them consumes AIRR
*Cell* objects; cell-level state lives in ``adata.obs`` / ``Dandelion.metadata`` / Seurat
``meta.data``. So the interop surface here is one emitter (:func:`to_airr`) and one
inverse (:func:`from_airr`); every bridge in :mod:`vdjtools.sc` is a thin adapter on top,
and the AIRR *Cell* Data File (:func:`vdjtools.sc.write_airr_cell`) stays a separate,
spec-faithful export rather than an interop path.

Two spellings are reconciled here, in exactly one place:

* vdjtools calls the nucleotide junction ``junction_nt``; AIRR — and therefore every
  downstream tool — calls it ``junction``.
* ``consensus_count`` is what scRepertoire's AIRR parser reads as the read count, while
  scirpy and dandelion prefer ``umi_count``. Both are emitted, so one file feeds all three.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..io.schema import (
    C_CALL,
    COUNT,
    D_CALL,
    J_CALL,
    JUNCTION_AA,
    JUNCTION_NT,
    LOCUS,
    V_CALL,
)
from .read import CELL_ID, CLONE_ID, PRODUCTIVE, SEQUENCE_ID, UMI_COUNT

#: AIRR spelling of the nucleotide junction (vdjtools stores it as ``junction_nt``).
JUNCTION = "junction"
#: Read count backing the contig consensus; scRepertoire's AIRR parser reads this one.
CONSENSUS_COUNT = "consensus_count"

#: Columns of the emitted AIRR Rearrangement table, in order.
AIRR_COLUMNS: list[str] = [
    SEQUENCE_ID, CELL_ID, LOCUS,
    V_CALL, D_CALL, J_CALL, C_CALL,
    JUNCTION, JUNCTION_AA, PRODUCTIVE,
    COUNT, UMI_COUNT, CONSENSUS_COUNT, CLONE_ID,
]


def _opt(df: pl.DataFrame, name: str, dtype: pl.DataType) -> pl.Expr:
    """``pl.col(name)`` when present, else a typed all-null literal under that name."""
    return (pl.col(name) if name in df.columns else pl.lit(None, dtype=dtype)).alias(name)


def to_airr(cells: pl.DataFrame) -> pl.DataFrame:
    """Convert the canonical sc long frame to a flat AIRR Rearrangement table.

    This is the contract every downstream bridge is built on. ``sequence_id`` is
    synthesised as ``<cell_id>_contig_<n>`` when the frame does not carry one — dandelion
    derives ``cell_id`` back out of exactly that pattern when a file omits it, so the
    round-trip closes either way.

    ``consensus_count`` is set from ``duplicate_count`` (the read count): for a 10x contig
    the reads backing the consensus *are* the duplicate observations, and emitting both
    spellings is what lets a single file satisfy scRepertoire (which reads
    ``consensus_count``) and scirpy/dandelion (which prefer ``umi_count``) at once.

    Args:
        cells: Single-cell long frame — :data:`vdjtools.sc.read.SC_COLUMNS`, as returned by
            :func:`~vdjtools.sc.read_10x` or :func:`~vdjtools.sc.read_airr_cell`. Missing
            optional columns are emitted as nulls.

    Returns:
        A ``pl.DataFrame`` with columns :data:`AIRR_COLUMNS`, one row per contig.

    Raises:
        ValueError: If ``cells`` has no ``cell_id`` column.
    """
    if CELL_ID not in cells.columns:
        raise ValueError(f"to_airr: frame has no {CELL_ID!r} column; got {cells.columns}")

    have_seq_id = SEQUENCE_ID in cells.columns
    # Per-cell contig ordinal, used only where a sequence_id is missing.
    ordinal = (pl.int_range(pl.len(), dtype=pl.UInt32).over(CELL_ID) + 1).cast(pl.Utf8)
    synthetic = pl.col(CELL_ID).cast(pl.Utf8) + "_contig_" + ordinal
    seq_id = pl.coalesce(pl.col(SEQUENCE_ID), synthetic) if have_seq_id else synthetic

    count = _opt(cells, COUNT, pl.Int64)
    return cells.select(
        seq_id.alias(SEQUENCE_ID),
        pl.col(CELL_ID).cast(pl.Utf8),
        _opt(cells, LOCUS, pl.Utf8),
        _opt(cells, V_CALL, pl.Utf8), _opt(cells, D_CALL, pl.Utf8),
        _opt(cells, J_CALL, pl.Utf8), _opt(cells, C_CALL, pl.Utf8),
        _opt(cells, JUNCTION_NT, pl.Utf8).alias(JUNCTION),   # AIRR spelling
        _opt(cells, JUNCTION_AA, pl.Utf8),
        _opt(cells, PRODUCTIVE, pl.Boolean),
        count,
        _opt(cells, UMI_COUNT, pl.Int64),
        count.alias(CONSENSUS_COUNT),
        _opt(cells, CLONE_ID, pl.Utf8),
    )


def from_airr(airr: pl.DataFrame) -> pl.DataFrame:
    """Convert a flat AIRR Rearrangement table back to the canonical sc long frame.

    The inverse of :func:`to_airr`, and the shared tail of every ``from_*`` bridge
    (:func:`~vdjtools.sc.from_scirpy`, :func:`~vdjtools.sc.from_dandelion`) — they each
    reduce their container to an AIRR frame and hand it here.

    Args:
        airr: AIRR Rearrangement frame with at least ``cell_id``. Both the AIRR
            ``junction`` and the vdjtools ``junction_nt`` spellings are accepted, as are
            ``consensus_count`` / ``duplicate_count`` for the read count.

    Returns:
        A ``pl.DataFrame`` in the canonical layout (:data:`vdjtools.sc.read.SC_COLUMNS`).

    Raises:
        ValueError: If ``airr`` has no ``cell_id`` column.
    """
    from .read import SC_COLUMNS

    if CELL_ID not in airr.columns:
        raise ValueError(f"from_airr: frame has no {CELL_ID!r} column; got {airr.columns}")
    cols = airr.columns

    junction_nt = (pl.col(JUNCTION_NT) if JUNCTION_NT in cols
                   else pl.col(JUNCTION) if JUNCTION in cols
                   else pl.lit(None, dtype=pl.Utf8))
    count = (pl.col(COUNT) if COUNT in cols
             else pl.col(CONSENSUS_COUNT) if CONSENSUS_COUNT in cols
             else pl.lit(None, dtype=pl.Int64))
    return airr.select(
        pl.col(CELL_ID).cast(pl.Utf8),
        _opt(airr, SEQUENCE_ID, pl.Utf8),
        _opt(airr, LOCUS, pl.Utf8),
        _opt(airr, V_CALL, pl.Utf8), _opt(airr, D_CALL, pl.Utf8),
        _opt(airr, J_CALL, pl.Utf8), _opt(airr, C_CALL, pl.Utf8),
        _opt(airr, JUNCTION_AA, pl.Utf8),
        junction_nt.cast(pl.Utf8).alias(JUNCTION_NT),
        count.cast(pl.Int64).alias(COUNT),
        _opt(airr, UMI_COUNT, pl.Int64),
        _opt(airr, CLONE_ID, pl.Utf8),
        _opt(airr, PRODUCTIVE, pl.Boolean),
    ).select(SC_COLUMNS)


def write_airr(cells: pl.DataFrame, path: str | Path) -> Path:
    """Write the sc frame as an AIRR Rearrangement TSV.

    The result is readable by :func:`~vdjtools.sc.read_airr_cell`, by scirpy
    (``ir.io.read_airr``), by dandelion (``ddl.read_airr``) and by scRepertoire
    (``loadContigs(format="AIRR")``).

    Args:
        cells: Single-cell long frame (:data:`vdjtools.sc.read.SC_COLUMNS`).
        path: Destination ``.tsv``.

    Returns:
        The path written.
    """
    out = Path(path)
    to_airr(cells).write_csv(out, separator="\t")
    return out


#: Columns scRepertoire's ``.parseAIRR`` requires, in the order it expects them.
_SCREPERTOIRE_AIRR: list[str] = [
    CELL_ID, LOCUS, CONSENSUS_COUNT, V_CALL, D_CALL, J_CALL, C_CALL, JUNCTION, JUNCTION_AA,
]

#: 10x contig-annotation columns scRepertoire's ``.parse10x`` consumes.
_SCREPERTOIRE_10X: list[tuple[str, str]] = [
    (CELL_ID, "barcode"), (SEQUENCE_ID, "contig_id"), (LOCUS, "chain"),
    (V_CALL, "v_gene"), (D_CALL, "d_gene"), (J_CALL, "j_gene"), (C_CALL, "c_gene"),
    (JUNCTION_AA, "cdr3"), (JUNCTION, "cdr3_nt"),
    (COUNT, "reads"), (UMI_COUNT, "umis"), (PRODUCTIVE, "productive"),
]


def write_screpertoire(cells: pl.DataFrame, path: str | Path, *,
                       format: str = "airr") -> Path:
    """Write a file R's scRepertoire can load, for ``combineTCR`` / ``combineBCR``.

    Export only — no R code ships with vdjtools. Load the result with
    ``scRepertoire::loadContigs(path, format = "AIRR")`` (or ``"10X"``).

    NOTE: two footguns on the R side, neither of which this function can prevent.
    ``.parseAIRR`` reads the read count from **``consensus_count``**, not ``umi_count`` —
    which is why :func:`to_airr` emits both. And ``combineTCR(samples=, ID=)`` rewrites
    barcodes to ``sample_ID_barcode``, the usual cause of a silent barcode-join failure
    against a Seurat ``meta.data``; pass the same ``samples``/``ID`` there and to
    ``combineExpression``, or omit both.

    Args:
        cells: Single-cell long frame (:data:`vdjtools.sc.read.SC_COLUMNS`).
        path: Destination — ``.tsv`` for ``format="airr"``, ``.csv`` for ``format="10x"``.
        format: ``"airr"`` writes an ``airr_rearrangement.tsv``-shaped file (the columns
            ``.parseAIRR`` requires); ``"10x"`` writes a
            ``filtered_contig_annotations.csv``-shaped file.

    Returns:
        The path written.

    Raises:
        ValueError: If ``format`` is not ``"airr"`` or ``"10x"``.
    """
    out = Path(path)
    airr = to_airr(cells)
    if format == "airr":
        airr.select(_SCREPERTOIRE_AIRR).write_csv(out, separator="\t")
    elif format == "10x":
        # .parse10x drops chain == "Multi", non-productive rows, and cdr3 == "None".
        airr.select(
            [pl.col(src).alias(dst) for src, dst in _SCREPERTOIRE_10X]
        ).write_csv(out)
    else:
        raise ValueError(f"format must be 'airr' or '10x'; got {format!r}")
    return out
