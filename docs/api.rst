API reference
=============

Every public subpackage is documented below: the native model engine (``vdjtools.model``),
IO/schema (``vdjtools.io``), statistics (``vdjtools.stats``), CDR3 features
(``vdjtools.features``), overlap/TCRnet (``vdjtools.overlap``), preprocessing
(``vdjtools.preprocess``), biomarker association (``vdjtools.biomarker``), repertoire dynamics (``vdjtools.dynamics``), single-cell interop
(``vdjtools.sc``), and the command-line interface (``vdjtools.cli``).

vdjtools
--------

.. automodule:: vdjtools
   :members:
   :undoc-members:
   :show-inheritance:

Model engine (``vdjtools.model``)
---------------------------------

The native V(D)J recombination engine: a model is a directory of tidy ``polars`` marginal
tables plus a ``manifest.json`` declaring the recombination Bayes net. It supersedes OLGA
(generation probability, sampling) and IGoR (EM inference), adds tandem-D (D-D) support, and
exposes information-theoretic diagnostics.

Model container, schema and events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.model
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.schema
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.events
   :members:
   :undoc-members:
   :show-inheritance:

Import, germline reference and stitching
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.io
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.bundled
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.reference
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.stitch
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.data
   :members:
   :undoc-members:
   :show-inheritance:

Generation probability, sampling and inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.pgen
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.native
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.generate
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.infer
   :members:
   :undoc-members:
   :show-inheritance:

Tandem-D (D-D) extension
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.dd
   :members:
   :undoc-members:
   :show-inheritance:

Model diagnostics
~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.analyze
   :members:
   :undoc-members:
   :show-inheritance:


Input/output and schema (``vdjtools.io``)
-----------------------------------------

The canonical clonotype frame (AIRR Rearrangement column names + polars dtypes) and every reader that emits it: native vdjtools tables, AIRR Rearrangement TSV, and Parquet, plus format auto-detection, metadata-driven batch loading, and hive-partitioned cohort scans.

``vdjtools.io.schema``
~~~~~~~~~~~~~~~~~~~~~~

Canonical clonotype schema and coercion helpers.

.. automodule:: vdjtools.io.schema
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.io.read``
~~~~~~~~~~~~~~~~~~~~

Single-file readers (native vdjtools, AIRR, Parquet).

.. automodule:: vdjtools.io.read
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.io.batch``
~~~~~~~~~~~~~~~~~~~~~

Format sniffing and metadata-driven batch / streaming reads.

.. automodule:: vdjtools.io.batch
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.io.cohort``
~~~~~~~~~~~~~~~~~~~~~~

Hive-partitioned Parquet cohort writer / lazy scanner.

.. automodule:: vdjtools.io.cohort
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.io.convert``
~~~~~~~~~~~~~~~~~~~~~~~~~

Converters for third-party repertoire formats (MiXcr, MiGec, Adaptive immunoSEQ,
IMGT/HighV-QUEST, Vidjil, RTCR) to the canonical clonotype frame.

.. automodule:: vdjtools.io.convert
   :members:
   :undoc-members:
   :show-inheritance:

Repertoire statistics (``vdjtools.stats``)
------------------------------------------

Diversity estimators (observed, Chao1/ChaoE, Shannon/Simpson, d50, Efron-Thisted) with exact, resampled, rarefaction and quantile variants; spectratype; and V/J/VJ segment usage.

``vdjtools.stats.diversity``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diversity indices (exact, resampled, quantile).

.. automodule:: vdjtools.stats.diversity
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.stats.rarefaction``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rarefaction / extrapolation curves.

.. automodule:: vdjtools.stats.rarefaction
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.stats.inext``
~~~~~~~~~~~~~~~~~~~~~~~~

iNEXT-style Hill-number rarefaction and extrapolation.

.. automodule:: vdjtools.stats.inext
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.stats.spectratype``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CDR3 length spectratype (nt / aa, weighted / unweighted).

.. automodule:: vdjtools.stats.spectratype
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.stats.usage``
~~~~~~~~~~~~~~~~~~~~~~~~

V / J / VJ segment-usage vectors and matrices.

.. automodule:: vdjtools.stats.usage
   :members:
   :undoc-members:
   :show-inheritance:

CDR3 features (``vdjtools.features``)
-------------------------------------

Per-clonotype and sample-level CDR3 sequence features: amino-acid physicochemical region profiles and k-mer / V+k-mer summaries.

``vdjtools.features.physchem``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Physicochemical CDR3 region profiles (Kidera factors, charge, hydropathy, ...).

.. automodule:: vdjtools.features.physchem
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.features.kmer``
~~~~~~~~~~~~~~~~~~~~~~~~~~

k-mer and V+k-mer occurrence summaries.

.. automodule:: vdjtools.features.kmer
   :members:
   :undoc-members:
   :show-inheritance:

Repertoire overlap and TCRnet (``vdjtools.overlap``)
----------------------------------------------------

Sample overlap and TCRnet built on vdjmatch / seqtree, including sequence-similarity-aware (TINA / Leinster-Cobbold) overlap, pairwise-distance matrices, clustering / MDS, and tracking.

``vdjtools.overlap.metrics``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exact-match overlap metrics (F, D, Jaccard, Morisita-Horn, ...).

.. automodule:: vdjtools.overlap.metrics
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.overlap.similarity``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sequence-similarity-weighted (TINA / Leinster-Cobbold) overlap.

.. automodule:: vdjtools.overlap.similarity
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.overlap.fuzzy``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Fuzzy (mismatch-tolerant) overlap via vdjmatch / seqtree.

.. automodule:: vdjtools.overlap.fuzzy
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.overlap.tcrnet``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TCRnet neighbourhood-enrichment degree statistics.

.. automodule:: vdjtools.overlap.tcrnet
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.overlap.alice``
~~~~~~~~~~~~~~~~~~~~~~~~~~

The same neighbourhood-enrichment test against a V(D)J *generation model* (Pgen null) rather than a control repertoire — the complement of ``tcrnet``.

.. automodule:: vdjtools.overlap.alice
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.overlap.cluster``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pairwise-distance matrices, clustering and MDS.

.. automodule:: vdjtools.overlap.cluster
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.overlap.track``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Clonotype tracking across a sample series.

.. automodule:: vdjtools.overlap.track
   :members:
   :undoc-members:
   :show-inheritance:

Preprocessing (``vdjtools.preprocess``)
---------------------------------------

Sample preprocessing and operations: downsampling, frequency error-correction, decontamination, segment/frequency/functional filters, VJ-usage batch-effect correction, and pooling / joining.

``vdjtools.preprocess.downsample``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read / unique / frequency downsampling (numpy multinomial).

.. automodule:: vdjtools.preprocess.downsample
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.preprocess.filter``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Segment, frequency and functional (coding) filters.

.. automodule:: vdjtools.preprocess.filter
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.preprocess.correct``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Frequency-based sequencing-error correction.

.. automodule:: vdjtools.preprocess.correct
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.preprocess.decontaminate``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cross-sample contamination removal.

.. automodule:: vdjtools.preprocess.decontaminate
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.preprocess.pool``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pool samples into a joint clonotype table.

.. automodule:: vdjtools.preprocess.pool
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.preprocess.join``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Join samples on shared clonotypes.

.. automodule:: vdjtools.preprocess.join
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.preprocess.batch``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VJ-usage batch-effect correction.

.. automodule:: vdjtools.preprocess.batch
   :members:
   :undoc-members:
   :show-inheritance:

Repertoire dynamics (``vdjtools.dynamics``)
-------------------------------------------

Paired **within-donor** clonotype testing across timepoints (Ayestaran 2024): per-pair effective sample size, deterministic downscale, two-tailed Fisher, and a five-way classification into emergent / expanded / persistent / contracted / vanishing. The sibling of ``vdjtools.biomarker`` — that tests incidence across *subjects*, this tests frequency across *timepoints* — and the complement of ``vdjtools.overlap.tcrnet`` / ``alice``, which measure breadth where this measures magnitude.

``vdjtools.dynamics.paired``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per-pair ``N_eff`` estimation and the paired clonotype test.

.. automodule:: vdjtools.dynamics.paired
   :members:
   :undoc-members:
   :show-inheritance:

Biomarker association (``vdjtools.biomarker``)
----------------------------------------------

Incidence-based clonotype-association testing across a cohort of repertoires (Emerson 2017, Howie 2015, De Witt 2018, Vlasova 2026): feature-vs-condition association (Fisher / χ² / Bayesian / permutation; binary, category, or Cochran–Mantel–Haenszel stratified conditions) and feature-vs-feature co-occurrence (α-β pairing, same-chain co-specificity), with exact or 1-mismatch matching and metaclonotype grouping.

``vdjtools.biomarker.association``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

General incidence association and candidate selection.

.. automodule:: vdjtools.biomarker.association
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.biomarker.cooccurrence``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Feature-vs-feature co-occurrence (in-silico α-β pairing, same-chain co-specificity).

.. automodule:: vdjtools.biomarker.cooccurrence
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.biomarker.condition``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phenotype-design builders (binary, categorical, HLA alleles, zygosity, CMH strata).

.. automodule:: vdjtools.biomarker.condition
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.biomarker.stats``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vectorised 2×2 test kernels (Fisher, χ², Bayesian, CMH, permutation, FDR).

.. automodule:: vdjtools.biomarker.stats
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.biomarker.fisher``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Fisher-exact incidence association (V/J-match, exact / 1mm) — the Emerson-2017 shortcut.

.. automodule:: vdjtools.biomarker.fisher
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.biomarker.metaclonotype``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Metaclonotype grouping (fuzzy CDR3 + V/J).

.. automodule:: vdjtools.biomarker.metaclonotype
   :members:
   :undoc-members:
   :show-inheritance:

Single-cell interop (``vdjtools.sc``)
-------------------------------------

Single-cell AIRR Cell / 10x paired-chain interop: contig ingestion, chain resolution and pairing, doublet / mispairing QC, cluster evaluation, and an AnnData bridge.

``vdjtools.sc.read``
~~~~~~~~~~~~~~~~~~~~

10x / AIRR-Cell ingestion and AIRR Data File export.

.. automodule:: vdjtools.sc.read
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.sc.pair``
~~~~~~~~~~~~~~~~~~~~

Chain resolution, alpha/beta pairing and mispairing flags.

.. automodule:: vdjtools.sc.pair
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.sc.pgen``
~~~~~~~~~~~~~~~~~~~~~

Paired-chain generation probability (``Pgen(α)·Pgen(β)``) via the native model.

.. automodule:: vdjtools.sc.pgen
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.sc.cluster_eval``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clustering-quality evaluation against antigen labels.

.. automodule:: vdjtools.sc.cluster_eval
   :members:
   :undoc-members:
   :show-inheritance:

``vdjtools.sc.anndata``
~~~~~~~~~~~~~~~~~~~~~~~

AnnData bridge for the single-cell frame.

.. automodule:: vdjtools.sc.anndata
   :members:
   :undoc-members:
   :show-inheritance:

Command-line interface (``vdjtools.cli``)
-----------------------------------------

The unified ``vdjtools`` typer CLI (pgen, generate, diversity, spectratype, segment-usage, overlap, models).

``vdjtools.cli``
~~~~~~~~~~~~~~~~

The ``vdjtools`` command-line application.

.. automodule:: vdjtools.cli
   :members:
   :undoc-members:
   :show-inheritance:

