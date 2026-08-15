"""Bridge the single-cell frame to and from dandelion.

A ``Dandelion`` object is two tables: ``.data``, the contig-level AIRR frame indexed by
``sequence_id``, and ``.metadata``, one row per cell. The first is exactly what
:func:`~vdjtools.sc.to_airr` emits, so the bridge is a format hand-off rather than a
translation — dandelion's own scirpy converter round-trips through the same flat AIRR
table.

:func:`read_h5ddl` is the useful part: dandelion persists to ``.h5ddl``, which is plain
HDF5 (h5py structured arrays, with a sibling Zarr store for distances), so a dandelion
result is readable **without installing dandelion**.

NOTE: dandelion also ships a polars backend (``ddl.set_backend("polars")``,
``DandelionPolars``) whose ``.data``/``.metadata`` take polars frames directly. If you are
on that backend, hand :func:`~vdjtools.sc.to_airr` output straight over; the pandas
conversion here is only for the default backend.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from .airr import from_airr, to_airr
from .read import CELL_ID, SEQUENCE_ID


def to_dandelion(cells: pl.DataFrame):
    """Wrap the sc long frame as a :class:`dandelion.Dandelion`.

    ``sequence_id`` and ``umi_count`` are both required by dandelion's loader;
    :func:`~vdjtools.sc.to_airr` guarantees the first (synthesising ``<cell>_contig_<n>``
    when absent) and carries the second.

    Args:
        cells: Single-cell long frame (:data:`vdjtools.sc.read.SC_COLUMNS`).

    Returns:
        A ``Dandelion`` whose ``.data`` is the contig-level AIRR table; ``.metadata`` is
        built by dandelion's own ``update_metadata()``.

    Raises:
        ImportError: If dandelion is not installed (``pip install sc-dandelion``).
    """
    try:
        import dandelion as ddl
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError(
            "to_dandelion needs dandelion: pip install sc-dandelion "
            "(read_h5ddl needs only h5py)"
        ) from e

    from .anndata import _to_pandas

    return ddl.Dandelion(data=_to_pandas(to_airr(cells)))


def from_dandelion(vdj) -> pl.DataFrame:
    """Convert a :class:`dandelion.Dandelion` back to the canonical sc long frame.

    Args:
        vdj: A ``Dandelion`` object (or anything exposing a contig-level ``.data`` frame).

    Returns:
        A ``pl.DataFrame`` in the canonical layout (:data:`vdjtools.sc.read.SC_COLUMNS`).

    Raises:
        AttributeError: If ``vdj`` has no ``.data``.
    """
    from .anndata import _from_pandas

    data = vdj.data if hasattr(vdj, "data") else vdj
    df = data if isinstance(data, pl.DataFrame) else _from_pandas(data)
    return from_airr(_ensure_cell_id(df))


def _ensure_cell_id(df: pl.DataFrame) -> pl.DataFrame:
    """Derive ``cell_id`` from ``<cell>_contig_<n>`` when a frame omits it.

    This mirrors dandelion's own ``load_data``, which splits ``sequence_id`` on
    ``"_contig"`` — so a table that went out through :func:`to_dandelion` comes back keyed
    the same way even if the cell column was dropped in between.
    """
    if CELL_ID in df.columns or SEQUENCE_ID not in df.columns:
        return df
    return df.with_columns(
        pl.col(SEQUENCE_ID).cast(pl.Utf8).str.split("_contig").list.first().alias(CELL_ID)
    )


def read_h5ddl(path: str | Path) -> pl.DataFrame:
    """Read the contig table out of a dandelion ``.h5ddl`` file, without dandelion.

    ``.h5ddl`` is plain HDF5: ``data`` is a h5py structured array of the contig-level AIRR
    table (``metadata`` holds the cell-level one, and distances may live in a sibling Zarr
    store). Only ``h5py`` is needed to read it.

    Args:
        path: Path to a ``.h5ddl`` file.

    Returns:
        A ``pl.DataFrame`` in the canonical layout (:data:`vdjtools.sc.read.SC_COLUMNS`).

    Raises:
        ImportError: If ``h5py`` is not installed.
        KeyError: If the file has no ``data`` dataset.
    """
    try:
        import h5py
    except ImportError as e:  # pragma: no cover - optional dep
        raise ImportError("read_h5ddl needs h5py: pip install h5py") from e

    with h5py.File(Path(path), "r") as fh:
        if "data" not in fh:
            raise KeyError(f"{path!r} has no 'data' dataset; is it a .h5ddl file?")
        arr = fh["data"][:]
        names = arr.dtype.names or ()
        cols = {}
        for name in names:
            col = arr[name]
            # h5py returns fixed-width bytes for the string columns; decode to str.
            cols[name] = ([v.decode() if isinstance(v, bytes) else v for v in col]
                          if col.dtype.kind in "SO" else col)
    return from_airr(_ensure_cell_id(pl.DataFrame(cols)))
