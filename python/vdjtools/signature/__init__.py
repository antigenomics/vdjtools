"""Repertoire signatures: raw features, a corpus-fitted rotation, and the emitted vector.

The public stages:

* :mod:`~vdjtools.signature.layout` -- the column contract. Which raw features exist, what each
  one's support is, which columns pass through as channels. Computes nothing, imports nothing heavy.
* :mod:`~vdjtools.signature.transform` -- the variance-stabilising transforms, applied where a raw
  feature is computed so its denominator is still in scope.
* :mod:`~vdjtools.signature.features` -- the raw features and channels for one sample. A pure
  function of that sample plus the germline vocabulary.
* :mod:`~vdjtools.signature.corpus` -- the only module that needs a corpus: build one, winsorize it,
  fit the rotation and the per-PC scaling, write and read the artifact.

The internal ``bounded`` module supplies sparse fitting workspaces to ``corpus``.

:func:`vsig` and :func:`vsig_cohort` in :mod:`~vdjtools.signature.signature` put these stages together.
The geometry half (``rsig``) lives in ``mir.signature`` and registers its groups into this same
layout registry; nothing here imports ``mir``.
"""
from .corpus import DEFAULT_COMPONENTS, MODES, WINSOR_PS, Corpus, fit, synthesize
from .features import PAIRS, WEIGHTS, gene_vocab, raw_and_channels, sanitise, work_frame
from .layout import (
    AMINO_ACIDS,
    LOCI,
    NO_LOCUS,
    PC_BLOCK,
    SPECTRATYPE_LENGTHS,
    SUPPORTS,
    TRANSFORMS,
    Channel,
    RawGroup,
    channel_columns,
    channels,
    feats,
    feature_definition,
    validate_feature_definition,
    parse,
    pc_columns,
    raw_columns,
    raw_groups,
    register_channel,
    register_raw,
    signature_columns,
    support_of,
)
from .signature import raw_features_cohort, vsig, vsig_cohort
from .transform import arcsine, clr, log1p, log10, logit

__all__ = [
    "DEFAULT_COMPONENTS", "MODES", "WINSOR_PS", "Corpus", "fit", "raw_and_channels", "synthesize",
    "vsig", "vsig_cohort", "raw_features_cohort",
    "AMINO_ACIDS", "LOCI", "NO_LOCUS", "PAIRS", "PC_BLOCK", "SPECTRATYPE_LENGTHS", "SUPPORTS",
    "TRANSFORMS", "WEIGHTS", "Channel", "RawGroup", "arcsine", "channel_columns", "channels",
    "clr", "feats", "feature_definition", "validate_feature_definition", "gene_vocab", "log10", "log1p", "logit", "parse", "pc_columns",
    "raw_columns", "raw_groups", "register_channel", "register_raw", "sanitise",
    "signature_columns", "support_of", "work_frame",
]
