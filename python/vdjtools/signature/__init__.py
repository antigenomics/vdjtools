"""vdjtools.signature — the repertoire-signature contract, and the statistics half of it.

See :mod:`vdjtools.signature.layout` for what a signature is and how columns are named.
"""
from .layout import (
    LOCI,
    NO_LOCUS,
    TIERS,
    TRANSFORMS,
    Block,
    columns,
    describe,
    feats,
    index,
    parse,
    register,
    registry,
)
from .assemble import DEFAULT_CSTAR, vsig, vsig_cohort
from .transform import (
    DEFAULT_CLIP,
    arcsine,
    clr,
    log1p,
    log10,
    logit,
    magnitude_scale,
    reference_z,
    robust_loc_scale,
)

__all__ = [
    "DEFAULT_CLIP",
    "DEFAULT_CSTAR",
    "LOCI",
    "NO_LOCUS",
    "TIERS",
    "TRANSFORMS",
    "Block",
    "arcsine",
    "registry",
    "clr",
    "columns",
    "describe",
    "feats",
    "index",
    "log10",
    "log1p",
    "logit",
    "magnitude_scale",
    "parse",
    "reference_z",
    "register",
    "robust_loc_scale",
    "vsig",
    "vsig_cohort",
]
