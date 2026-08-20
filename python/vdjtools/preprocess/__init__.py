"""vdjtools.preprocess — Downsampling, error-correction, filtering, batch-effect correction.

Free functions over the canonical clonotype frame (see :mod:`vdjtools.io.schema`).
"""
from .batch import apply_vj_correction, correct_vj_usage
from .correct import correct
from .decontaminate import decontaminate
from .downsample import downsample, select_top
from .filter import (
    filter_by_sample,
    filter_frequency,
    filter_functional,
    filter_length,
    filter_functional_genes,
    filter_productive,
    filter_segment,
    productive_mask,
)
from .join import join_samples
from .pool import pool_samples, resolve_key

__all__ = [
    "downsample",
    "select_top",
    "filter_productive",
    "filter_length",
    "filter_functional_genes",
    "productive_mask",
    "filter_functional",          # deprecated alias for filter_productive
    "filter_frequency",
    "filter_segment",
    "filter_by_sample",
    "correct",
    "decontaminate",
    "pool_samples",
    "join_samples",
    "resolve_key",
    "correct_vj_usage",
    "apply_vj_correction",
]
