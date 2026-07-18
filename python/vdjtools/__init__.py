"""vdjtools — TCR/BCR immune-repertoire analysis (v2; Python + C++).

A clean-room rewrite of the legacy Groovy/Java vdjtools, standardised on the AIRR
schema and polars DataFrames with minimal object-orientation, built on the
antigenomics ecosystem (``seqtree``, ``vdjmatch``, ``arda``).

Native hot loops — the V(D)J Pgen dynamic program, the generation sampler, and the
EM E-step — live in the compiled :mod:`vdjtools._core` extension, imported lazily by
:mod:`vdjtools.model`. Everything else is pure polars/numpy, and generic sequence
primitives (Hamming/edit distance, fuzzy search) are used straight from the
``seqtree`` / ``vdjmatch`` / ``arda`` dependencies rather than duplicated here.
Subpackages (``io``, ``model``, ``stats``, ``features``, ``overlap``, ``preprocess``,
``biomarker``, ``sc``) are imported explicitly by the caller, so ``import vdjtools``
never pays the cost of the compiled extension or heavy optional dependencies until a
feature that needs them is used.
"""

__version__ = "3.0.0"
__all__ = ["__version__"]
