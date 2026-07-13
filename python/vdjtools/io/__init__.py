"""vdjtools.io — AIRR Rearrangement + native vdjtools I/O (polars).

Readers return the canonical clonotype frame defined in :mod:`vdjtools.io.schema`.
"""
from . import schema
from .batch import iter_samples, read, read_metadata, read_samples, sniff_format
from .cohort import ingest_cohort, scan_cohort
from .convert import (
    read_arda,
    read_imgt,
    read_immunoseq,
    read_migec,
    read_mixcr,
    read_rtcr,
    read_trust4,
    read_vidjil,
)
from .read import read_airr, read_parquet, read_vdjtools
from .schema import (
    C_CALL,
    JUNCTION_AA,
    JUNCTION_NT,
    COLUMNS,
    COUNT,
    D_CALL,
    FREQ,
    J_CALL,
    LOCUS,
    SCHEMA,
    V_CALL,
    add_locus,
    locus_of,
    normalize,
    recompute_frequency,
)

__all__ = [
    "schema",
    "read",
    "read_airr",
    "read_vdjtools",
    "read_parquet",
    "read_mixcr",
    "read_migec",
    "read_immunoseq",
    "read_imgt",
    "read_vidjil",
    "read_rtcr",
    "read_trust4",
    "read_arda",
    "read_metadata",
    "read_samples",
    "iter_samples",
    "ingest_cohort",
    "scan_cohort",
    "sniff_format",
    "normalize",
    "add_locus",
    "locus_of",
    "recompute_frequency",
    "SCHEMA",
    "COLUMNS",
    "V_CALL",
    "D_CALL",
    "J_CALL",
    "C_CALL",
    "JUNCTION_AA",
    "JUNCTION_NT",
    "COUNT",
    "FREQ",
    "LOCUS",
]
