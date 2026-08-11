"""vdjtools.signature — the repertoire-signature contract, and the statistics half of it.

See :mod:`vdjtools.signature.layout` for what a signature is and how columns are named.
"""
from .layout import (
    LOCI,
    NO_LOCUS,
    TIERS,
    Block,
    blocks,
    columns,
    describe,
    index,
    parse,
    register,
)

__all__ = [
    "LOCI",
    "NO_LOCUS",
    "TIERS",
    "Block",
    "blocks",
    "columns",
    "describe",
    "index",
    "parse",
    "register",
]
