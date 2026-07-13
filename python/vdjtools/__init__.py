"""vdjtools — TCR/BCR immune-repertoire analysis (v2; Python + C++).

A clean-room rewrite of the legacy Groovy/Java vdjtools, standardised on the AIRR
schema and polars DataFrames with minimal object-orientation, built on the
antigenomics ecosystem (``seqtree``, ``vdjmatch``, ``arda``).

Native hot loops — the V(D)J Pgen dynamic program, the generation sampler, and the
EM E-step — live in the compiled :mod:`vdjtools._core` extension. Everything else is
pure polars/numpy. Subpackages (``io``, ``model``, ``stats``, ``features``, ``overlap``,
``preprocess``, ``biomarker``, ``sc``) are imported explicitly by the caller so that
``import vdjtools`` never pays the cost of heavy optional dependencies (arda/mmseqs2,
vdjmatch/seqtree) until a feature that needs them is used.
"""
from ._core import hamming

__version__ = "2.2.0"
__all__ = ["hamming", "__version__"]
