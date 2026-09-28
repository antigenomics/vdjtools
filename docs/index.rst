vdjtools
========

.. rst-class:: lead

   Immune-repertoire analysis for T- and B-cell receptor sequencing: read any format, measure
   diversity and overlap, correct batch effects, score generation probability against a native
   V(D)J model, and reduce a whole repertoire to one comparable feature vector.

vdjtools is a Python library and a command-line tool. It reads the output of every common
repertoire pipeline into one canonical `AIRR <https://docs.airr-community.org/>`_ table backed by
`polars <https://pola.rs>`_, and every analysis takes and returns that same table -- so results
chain together instead of needing a converter between each step.

.. grid:: 1 2 2 2
   :gutter: 3
   :margin: 4 4 0 0

   .. grid-item-card:: Get started
      :link: getting-started
      :link-type: doc

      Install, load a real sample, and get your first diversity and overlap numbers. Start here if
      you have not used vdjtools before.

   .. grid-item-card:: Guides
      :link: usage
      :link-type: doc

      One page per task: loading data, pre-processing, diversity, overlap, dynamics, biomarkers,
      single cell, and the recombination model workshop.

   .. grid-item-card:: Command reference
      :link: cli
      :link-type: doc

      Every command and every flag, with the defaults. No Python needed.

   .. grid-item-card:: Python API
      :link: api
      :link-type: doc

      Every module, class and function, with signatures and types.

   .. grid-item-card:: Worked examples
      :link: notebooks
      :link-type: doc

      Fifteen runnable notebooks, several reproducing a published result --- repertoire ageing,
      CMV and HLA association, vaccination time courses, single-cell pairing.

   .. grid-item-card:: Glossary
      :link: glossary
      :link-type: doc

      Every term this documentation uses, and which convention vdjtools follows where the field
      disagrees --- junction against CDR3, productive against functional.

Install
-------

.. code-block:: bash

   pip install vdjtools

That is the whole installation. Wheels are prebuilt for CPython 3.10-3.13 on Linux, macOS
(Apple Silicon) and Windows, with the native C++ extension included; a source install compiles it.

Everything the documentation describes works from that one command -- there are no extras to opt
into for core functionality. The three engines vdjtools delegates to are base dependencies, imported
lazily so that ``import vdjtools`` stays light:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Engine
     - Powers
   * - `arda <https://github.com/antigenomics/arda>`_
     - the germline reference (V/D/J segments and CDR3 anchors), the model engine, annotation
   * - `seqtree <https://github.com/antigenomics/seqtree>`_
     - fuzzy search and e-values: error correction, similarity overlap, TCRnet
   * - `vdjmatch <https://github.com/antigenomics/vdjmatch>`_
     - sample overlap, TCRnet, metaclonotypes

Two optional extras exist because each has a working alternative:

.. code-block:: bash

   pip install "vdjtools[overlap]"   # scikit-learn, for cluster_samples(method="mds")
   pip install "vdjtools[sc]"        # single-cell bridges: anndata, awkward, mudata, pyyaml

MMseqs2 is needed only for arda's alignment and annotation path, never for germline lookup, Pgen,
generation or the analytics.

Your first result
-----------------

.. tab-set::

   .. tab-item:: Python

      .. code-block:: python

         from vdjtools import io as vio, stats, overlap

         sample = vio.read("clones.tsv")          # MiXcr, immunoSEQ, AIRR, Parquet -- detected
         stats.diversity_stats(sample)            # richness, Chao1, Shannon, Simpson, d50
         stats.segment_usage(sample, "v")         # V-gene usage

         cohort = vio.read_samples(vio.read_metadata("metadata.txt"), base_dir="samples/")
         overlap.overlap_matrix(cohort)           # pairwise repertoire overlap

   .. tab-item:: Command line

      .. code-block:: bash

         vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
         vdjtools segment-usage *.tsv --segment v -o usage.tsv
         vdjtools overlap       *.tsv -o overlap.tsv

Every command writes to ``-o`` -- TSV, or Parquet when the path ends in ``.parquet`` or ``.pq`` --
or to stdout, so commands pipe. :doc:`getting-started` walks through this with real data.

Choose your path
----------------

.. list-table::
   :header-rows: 1
   :widths: 52 48

   * - You want to
     - Go to
   * - Get a first result from a file you already have
     - :doc:`getting-started`
   * - Clean, filter, downsample or batch-correct a cohort
     - :doc:`preprocessing`
   * - Measure diversity, overlap, gene usage or spectratype
     - :doc:`usage`
   * - Score generation probability, or fit a model to your own reads
     - :doc:`model`
   * - Reduce each sample to one fixed feature vector for a classifier
     - :doc:`signature`
   * - Work with 10x or single-cell data, or hand off to scirpy or dandelion
     - :doc:`singlecell`
   * - Look up a command, a flag or a function
     - :doc:`cli` | :doc:`api`
   * - Understand a term used in these pages
     - :doc:`glossary`

What is in the box
------------------

.. list-table::
   :header-rows: 1
   :widths: 20 56 24

   * - Module
     - What it does
     - Guide
   * - :mod:`vdjtools.io`
     - Readers for native vdjtools, AIRR Rearrangement TSV and Parquet; format-detecting converters
       for MiXcr, MiGec, Adaptive immunoSEQ, IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4 and arda;
       metadata-driven batches and hive-partitioned cohorts
     - :doc:`preprocessing`
   * - :mod:`vdjtools.stats`
     - Diversity (observed richness, Chao1, Efron-Thisted, Shannon, Simpson, d50), coverage-based
       Hill-number rarefaction with bootstrap intervals, spectratype, V/D/J/C usage
     - :doc:`usage`
   * - :mod:`vdjtools.model`
     - Native V(D)J recombination model: Pgen over nucleotides and amino acids, the Hamming-1 ball,
       sequence generation, EM inference, tandem-D support, and a full model workshop
     - :doc:`model`
   * - :mod:`vdjtools.overlap`
     - Exact and similarity-aware repertoire overlap, TCRnet neighbourhood enrichment, sample
       clustering
     - :doc:`usage`
   * - :mod:`vdjtools.preprocess`
     - Format conversion, the three filtering axes, frequency handling, downsampling, error
       correction, V/J-usage batch correction, pooling and joining
     - :doc:`preprocessing`
   * - :mod:`vdjtools.features`
     - Junction physicochemical profiles and k-mer summaries
     - :doc:`usage`
   * - :mod:`vdjtools.biomarker`
     - Incidence-based association against binary, per-HLA-allele or stratified conditions;
       co-occurrence pairing; metaclonotypes
     - :doc:`usage`
   * - :mod:`vdjtools.dynamics`
     - Longitudinal clonotype tracking between timepoints: the paired within-donor expansion test,
       the VDJtrack recapture model, an edgeR NB-exact caller
     - :doc:`usage`
   * - :mod:`vdjtools.signature`
     - One fixed, named, already-standardised feature vector per repertoire, rotated through a
       published corpus
     - :doc:`signature`
   * - :mod:`vdjtools.sc`
     - Single-cell ingestion, chain pairing and QC, paired Pgen, and round-trip interoperability
       with scirpy, dandelion and scRepertoire
     - :doc:`singlecell`

Performance
-----------

The Pgen, generation, EM and diversity hot paths are a native C++ core reached through pybind11;
everything else is polars. Amino-acid Pgen agrees with OLGA to 1e-15 on all seven human loci.

Single thread, Apple M3, bundled human TRB model:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Operation
     - Throughput
     - Relative to OLGA
   * - Nucleotide Pgen, single-D VDJ
     - 0.5 ms per sequence
     - **9x**
   * - Amino-acid Pgen
     - 0.6-0.9 ms per sequence
     - **8.6x**
   * - Amino-acid Pgen over the Hamming-1 ball
     - 15 ms per sequence
     - **8.7x**
   * - Sequence generation, raw draws
     - 32,100 per second
     - not applicable
   * - Sequence generation, ``productive_only=True``
     - 8,960 per second
     - not applicable

The productive filter costs 3.6x on TRB and 5.1x on IGH (19,900 raw draws per second against 3,930
productive), because out-of-frame and stop-codon draws are discarded and redrawn rather than
repaired. Quote whichever of the two your pipeline uses.

Batched Pgen over many junctions parallelises over sequences
(:func:`~vdjtools.model.native.pgen_aa_batch`, 11x on 16 cores, and identical to the serial result to the last bit);
the EM E-step parallelises over reads (6.7x on 8 threads). Memory stays modest: 63 MB resident for
``import vdjtools`` plus one model, 123 MB with all seven bundled models loaded.

Citing
------

For the v2 rewrite, cite this repository. For the original tool, cite
`Shugay et al., PLoS Computational Biology 2015 <https://doi.org/10.1371/journal.pcbi.1004503>`_.
The VDJtrack recapture model in :mod:`vdjtools.dynamics` is Pavlova, Zvyagin and Shugay (2024).

Licensed GPL-3.0-or-later. The legacy Groovy/Java v1.x tool lives on the
`legacy-1.x <https://github.com/antigenomics/vdjtools/tree/legacy-1.x>`_ branch, with its releases
under the repository tags ``v0.0.1`` through ``1.2.1``.

.. toctree::
   :hidden:
   :caption: Getting started
   :maxdepth: 2

   self
   getting-started
   glossary

.. toctree::
   :hidden:
   :caption: Guides
   :maxdepth: 2

   usage
   preprocessing
   signature
   model
   singlecell

.. toctree::
   :hidden:
   :caption: Reference
   :maxdepth: 2

   cli
   api
   channels

.. toctree::
   :hidden:
   :caption: How it works
   :maxdepth: 2

   signature-methods

.. toctree::
   :hidden:
   :caption: Worked examples
   :maxdepth: 2

   notebooks
