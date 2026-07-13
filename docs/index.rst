vdjtools
========

TCR/BCR immune-repertoire analysis — v2, a clean-room rewrite in Python + C++.

Standardised on the AIRR schema and polars DataFrames with minimal object-orientation,
built on the antigenomics ecosystem (`seqtree <https://github.com/antigenomics/seqtree>`_,
`vdjmatch <https://github.com/antigenomics/vdjmatch>`_,
`arda <https://github.com/antigenomics/arda>`_).

.. note::

   **v2.2.0** — the native V(D)J model engine plus the full analytics suite (diversity, overlap/TCRnet,
   preprocessing, biomarkers, single-cell), CDR features, and legacy-format ingestion (MiXcr, MiGec,
   immunoSEQ, IMGT/HighV-QUEST, Vidjil, RTCR). Clonotype columns follow the AIRR **junction**
   convention (``junction_nt`` / ``junction_aa``). The legacy Groovy/Java vdjtools (v1.x) lives on the
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
   native.pgen_aa_batch(model, seqs, mismatches=1, threads=0)  # many CDR3s, thread-parallel (~11x)
   generate(model, 1000)                              # sample a repertoire -> polars DataFrame

Matches OLGA's Pgen to machine precision across all 7 loci, adds tandem-D (D-D) support, and
learns models from your own reads (:func:`vdjtools.model.infer.infer_native`). Explore any model's
recombination Bayes net interactively with ``marimo edit notebooks/model_explorer.py``.

Command line
------------

``pip install vdjtools`` installs the ``vdjtools`` command — the model engine (OLGA/IGoR-style) and
the repertoire analytics (sample files or a metadata table, like the legacy tool):

.. code-block:: bash

   vdjtools models                                # list the bundled models
   vdjtools generate -m TRB -n 1000 -o gen.tsv    # sample sequences   (cf. olga-generate_sequences)
   vdjtools pgen seqs.tsv -m TRB -o pgen.tsv      # Pgen per CDR3       (cf. olga-compute_pgen)
   vdjtools pgen seqs.tsv -m TRB --mismatches 1   # + the Hamming-1 ball; --v-col/--j-col to condition

   vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
   vdjtools overlap       *.tsv -o overlap.tsv
   vdjtools segment-usage *.tsv --segment v -o usage.tsv
   vdjtools spectratype   *.tsv -o spectra.tsv

Native vdjtools and AIRR inputs are auto-detected; every command writes TSV to ``-o`` (or stdout).
Run ``vdjtools <command> --help`` for options.

Performance
-----------

The Pgen / generation / EM / diversity hot paths are a native C++ (pybind11) core; everything else is
polars. On an Apple M3 (single thread, bundled human TRB model): nucleotide Pgen **~0.5 ms/seq**
(**9× OLGA**, single-D VDJ), amino-acid Pgen **~0.6–0.9 ms/seq**
(**8.6× OLGA**, exact to 1e-15), the Hamming-1 ball **~15 ms/seq** (**8.7×**), sequence generation
**~32 000 seq/s**. Both nt and aa Pgen use the same transfer-matrix DP (an in-frame CDR3 is an aa
query with one codon fixed per position); nt is exact vs OLGA on all 7 loci. Batched Pgen / 1-mismatch
over many CDR3s parallelises over sequences (:func:`vdjtools.model.native.pgen_aa_batch`, ~11× on 16
cores, bitwise-identical to serial); the EM E-step parallelises over reads (~6.7× on 8 threads), and
diversity / rarefaction run on a native iNEXT kernel (bootstrap + parallel batch). Memory stays light —
**~63 MB** resident for ``import vdjtools`` plus one model, **~123 MB** with all seven bundled models.

Capabilities (see the :doc:`API reference <api>` and the project ROADMAP):

- **IO** — canonical AIRR **junction** clonotype frame (``junction_nt`` / ``junction_aa``); readers
  for native vdjtools, AIRR TSV, Parquet, and converters for MiXcr, MiGec, immunoSEQ (v1/v2),
  IMGT/HighV-QUEST, Vidjil, and RTCR (:mod:`vdjtools.io.convert`); metadata-driven batch + cohorts.
- **Model** — native V(D)J recombination model: Pgen (nt / aa / 1-mismatch / V/J-agnostic /
  thread-parallel batch), sequence generation, EM inference, and tandem-D (D-D) support — a native
  pybind11 core that supersedes OLGA and IGoR. Concordant with OLGA on all 7 loci; bundled models.
- **Stats** — diversity, spectratype, V/J/VJ usage.
- **Features** — CDR physicochemical profiles, k-mer / V+k-mer summaries.
- **Overlap** — sample overlap and TCRnet (via vdjmatch/seqtree), similarity-aware overlap, clustering.
- **Preprocess** — downsampling, error-correction, batch-effect correction, pooling/joining.
- **Biomarker** — incidence-based association (Fisher) and metaclonotype grouping.
- **Single-cell** — AIRR Cell / 10x interoperability, chain pairing + QC, and paired α/β Pgen.

.. toctree::
   :hidden:

   self
   api
