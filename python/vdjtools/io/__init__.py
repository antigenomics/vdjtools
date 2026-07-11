"""vdjtools.io — AIRR Rearrangement + native vdjtools I/O (polars).

Readers return the canonical clonotype frame defined in :mod:`vdjtools.io.schema`.
"""
from . import schema
from .batch import read, read_metadata, read_samples, sniff_format
from .read import read_airr, read_vdjtools
from .schema import (
    C_CALL,
    CDR3_AA,
    CDR3_NT,
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
    "read_metadata",
    "read_samples",
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
    "CDR3_AA",
    "CDR3_NT",
    "COUNT",
    "FREQ",
    "LOCUS",
]
