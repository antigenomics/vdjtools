vdjtools
========

TCR/BCR immune-repertoire analysis — v2, a clean-room rewrite in Python + C++.

Standardised on the AIRR schema and polars DataFrames with minimal object-orientation,
built on the antigenomics ecosystem (`seqtree <https://github.com/antigenomics/seqtree>`_,
`vdjmatch <https://github.com/antigenomics/vdjmatch>`_,
`arda <https://github.com/antigenomics/arda>`_).

.. warning::

   v2.0.0 is a pre-release under active development (latest ``v2.0.0-alpha.2`` — the native V(D)J
   model engine; ``pip install --pre vdjtools``). The legacy Groovy/Java vdjtools (v1.x) lives on the
   ``legacy-1.x`` branch and its releases remain available under the repository's tags.

Quickstart — recombination model engine
----------------------------------------

Precomputed models for all 7 human loci ship in the wheel (no OLGA or download needed):

.. code-block:: python

   from vdjtools.model import load_bundled, native
   from vdjtools.model.generate import generate

   model = load_bundled("TRB", source="olga")        # or source="learned" (fit to real repertoires)
   native.pgen_nt(model, "TGTGCCAGCAGC...")           # nucleotide generation probability (native C++)
   native.pgen_aa(model, "CASSLAPGATNEKLFF")          # amino-acid Pgen (codon-marginalised)
   native.pgen_aa(model, "CASSLAPGATNEKLFF", mismatches=1)  # + the whole Hamming-1 ball
   generate(model, 1000)                              # sample a repertoire -> polars DataFrame

Matches OLGA's Pgen to machine precision across all 7 loci, adds tandem-D (D-D) support, and
learns models from your own reads (:func:`vdjtools.model.infer.infer_native`). Explore any model's
recombination Bayes net interactively with ``marimo edit notebooks/model_explorer.py``.

Capabilities (rolling out by phase — see the project ROADMAP):

- **Model** — native V(D)J recombination model: Pgen (nt / aa / 1-mismatch / V/J-agnostic),
  sequence generation, EM inference, and tandem-D (D-D) support — a native pybind11 core that
  supersedes OLGA and IGoR. Concordant with OLGA on all 7 loci; bundled precomputed models.
- **Stats** — diversity, spectratype, V/J/VJ usage.
- **Features** — CDR physicochemical profiles, k-mer / V+k-mer summaries.
- **Overlap** — sample overlap and TCRnet (via vdjmatch/seqtree).
- **Preprocess** — downsampling, error-correction, batch-effect correction.
- **Biomarker** — incidence-based association (Fisher) and metaclonotype grouping.
- **Single-cell** — AIRR Cell / 10x interoperability.

.. toctree::
   :hidden:

   self
   api
