vdjtools
========

TCR/BCR immune-repertoire analysis — v2, a clean-room rewrite in Python + C++.

Standardised on the AIRR schema and polars DataFrames with minimal object-orientation,
built on the antigenomics ecosystem (`seqtree <https://github.com/antigenomics/seqtree>`_,
`vdjmatch <https://github.com/antigenomics/vdjmatch>`_,
`arda <https://github.com/antigenomics/arda>`_).

.. warning::

   v2.0.0 is under active development. The legacy Groovy/Java vdjtools (v1.x) lives on the
   ``legacy-1.x`` branch and its releases remain available under the repository's tags.

Capabilities (rolling out by phase — see the project ROADMAP):

- **Model** — native V(D)J recombination model: generation-probability (Pgen), sequence
  generation, and EM model inference (supersedes OLGA and IGoR).
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
