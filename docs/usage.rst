User guide
==========

Worked, copy-paste-runnable examples for every analysis module. All functions take and
return `polars <https://pola.rs>`_ ``DataFrame`` s on the canonical clonotype schema
(AIRR **junction** column names; see :mod:`vdjtools.io.schema`), so results chain together
and drop straight into plotting or ``.write_csv``.

A sample to play with
---------------------

Real workflows load a repertoire with ``io.read("sample.tsv")`` (see `Loading data`_).
So this page is runnable with no downloads, the snippet below draws a synthetic counted
sample straight from a bundled model — every later example reuses ``demo_sample``:

.. code-block:: python

   import polars as pl
   from vdjtools import io as vio
   from vdjtools.io import schema as S
   from vdjtools.model import load_bundled
   from vdjtools.model.generate import generate

   def demo_sample(locus="TRB", n=4000, seed=0):
       """A canonical clonotype frame sampled from the bundled OLGA model."""
       seqs = generate(load_bundled(locus, "olga"), n, seed=seed, productive_only=True)
       counted = (seqs.group_by(["junction_nt", "junction_aa", "v_call", "d_call", "j_call"])
                      .len().rename({"len": S.COUNT}))
       return vio.normalize(counted, recompute_freq=True)

   sample = demo_sample(seed=1)
   sample.columns
   # ['v_call', 'd_call', 'j_call', 'c_call', 'junction_aa', 'junction_nt',
   #  'duplicate_count', 'frequency']

Loading data
------------

:mod:`vdjtools.io` reads native vdjtools, AIRR Rearrangement TSV, and Parquet, and
auto-detects and converts the common third-party formats (MiXcr v1–4, MiGec, Adaptive
immunoSEQ v1/v2, IMGT/HighV-QUEST, Vidjil, RTCR). Every reader returns the same canonical
frame:

.. code-block:: python

   from vdjtools import io as vio

   vio.sniff_format("clones.txt")        # -> 'mixcr' | 'immunoseq' | 'airr' | 'vdjtools' | ...
   df = vio.read("clones.txt")           # auto-detect + convert to the canonical frame
   df = vio.read_immunoseq("adaptive.tsv")   # or call a converter explicitly

Load a whole cohort from a metadata sheet (one row per sample, plus any phenotype columns),
joining the metadata onto every clonotype:

.. code-block:: python

   meta = vio.read_metadata("metadata.txt")          # sample_name, disease_status, hla, ...
   cohort = vio.read_samples(meta, base_dir="samples/")   # one long frame, metadata attached
   by_id = vio.read_samples(meta, base_dir="samples/", as_dict=True)   # or {sample_id: frame}

For large cohorts, ``ingest_cohort`` writes a hive-partitioned Parquet store that
``scan_cohort`` reads back lazily (``pl.LazyFrame``), so you never hold every sample in
memory at once.

Repertoire statistics
---------------------

:mod:`vdjtools.stats` covers diversity, rarefaction/extrapolation, spectratype, and segment
usage. ``diversity_stats`` returns one row of estimators (observed richness, Chao1, ChaoE,
Efron–Thisted, Shannon, normalized Shannon, inverse Simpson, d50):

.. code-block:: python

   from vdjtools import stats

   stats.diversity_stats(sample)
   # columns: reads, observed_diversity, chao1, chaoE, efron_thisted,
   #          shannon_wiener, normalized_shannon_wiener, inverse_simpson, d50

   # iNEXT-style Hill-number rarefaction + extrapolation with bootstrap CIs
   stats.inext(sample, q=(0, 1, 2))         # order_q, m, method, sample_coverage, qD, qD_lo, qD_hi

   # V / J / VJ segment usage, and the CDR3-length spectratype
   stats.segment_usage(sample, "v")         # locus, v_call, weight
   stats.spectratype(sample, kind="aa")     # length spectrum, reads-weighted

Pass ``weight="unique"`` to any usage/spectratype call to weight by clonotype instead of
reads; ``inext_batch`` / ``rarefaction_batch`` run a whole cohort on the native parallel
kernel.

CDR3 features
-------------

:mod:`vdjtools.features` summarises CDR3 sequence content. ``physchem_profile`` gives the
mean amino-acid physicochemical properties (Kidera factors, charge, hydropathy, …) per group
and CDR3 region; ``kmer_profile`` counts k-mers:

.. code-block:: python

   from vdjtools import features

   features.physchem_profile(sample, region="all", group_by=("v_call",))
   # v_call, region, property, mean_value   (long format, one row per property)

   features.kmer_profile(sample, k=3)       # locus, kmer, weight
   features.v_kmer_c_profile(sample, k=3)   # V-anchored k-mer occurrences

Overlap and TCRnet
------------------

:mod:`vdjtools.overlap` compares samples. Exact-match overlap is pure polars; fuzzy /
similarity-aware overlap and TCRnet delegate to the vdjmatch + seqtree engine
(``pip install "vdjtools[overlap]"``):

.. code-block:: python

   from vdjtools import overlap

   a, b, c = demo_sample(seed=1), demo_sample(seed=2), demo_sample(seed=3)

   overlap.overlap_metrics(a, b)            # {'D':.., 'F':.., 'F2':.., 'R':.., 'd1':.., 'd2':.., 'd12':..}
   overlap.similarity_overlap(a, b)         # TINA / Leinster-Cobbold sequence-similarity overlap
   overlap.tcrnet(a)                        # per-clonotype neighbourhood enrichment (E, p_enrichment)

   # all-pairs distance matrix -> 2-D embedding for a cohort
   dist = overlap.pairwise_distances({"A": a, "B": b, "C": c}, metric="F")
   overlap.cluster_samples(dist, method="mds")     # sample, mds1, mds2

   # frequency trajectories of the top clonotypes across an ordered series
   overlap.track_clonotypes({"t0": a, "t1": b, "t2": c}, top=50)

Preprocessing
-------------

:mod:`vdjtools.preprocess` normalises samples before comparison — downsampling to a common
depth, error-correction, filtering, and pooling/joining:

.. code-block:: python

   from vdjtools import preprocess

   preprocess.downsample(sample, 1000)              # resample to 1000 reads (numpy multinomial)
   preprocess.filter_functional(sample)             # drop out-of-frame / stop-codon clonotypes
   preprocess.filter_frequency(sample, min_freq=1e-4)
   preprocess.correct(sample, max_mismatches=2)     # collapse likely sequencing errors

   # combine samples: pooled clonotype table, or an incidence/frequency join
   preprocess.pool_samples([a, b, c])
   preprocess.join_samples([a, b, c], min_samples=2)   # clonotypes seen in >=2 samples

``correct_vj_usage`` applies a VJ-usage batch-effect correction across a cohort.

Biomarker association
---------------------

:mod:`vdjtools.biomarker` finds clonotypes whose **incidence** (presence across subjects)
associates with a phenotype, by Fisher-exact test (the Emerson-2017 design):

.. code-block:: python

   import polars as pl
   from vdjtools.biomarker import fisher_association, metaclonotypes

   # cohort: long frame with a sample_id column; phenotype: one row per subject
   cohort = pl.DataFrame({
       "sample_id":   ["p0","p1","p2","n0","n1","n2"],
       "v_call":      ["TRBV1"]*6, "j_call": ["TRBJ1"]*6,
       "junction_aa": ["CASSXF","CASSXF","CASSXF","CASSXF","CASSBG","CASSBG"],
       "duplicate_count": [10]*6,
   })
   pheno = pl.DataFrame({"sample_id": ["p0","p1","p2","n0","n1","n2"],
                         "cmv": [True, True, True, False, False, False]})

   fisher_association(cohort, pheno, pheno_col="cmv")
   # per feature: incidence, n_pos_present, n_neg_present, direction, log2_or, p_value

   metaclonotypes(cohort)                 # group near-identical CDR3s (1-mismatch + V/J) -> meta_id

Use ``match="1mm"`` on ``fisher_association`` to pool single-mismatch neighbours of each
clonotype (needs the ``[overlap]`` engine).

Single-cell
-----------

:mod:`vdjtools.sc` ingests 10x / AIRR-Cell contigs into a flat ``cell_id``-keyed frame,
resolves and pairs chains with doublet / mispairing QC, and scores paired α/β generation
probability:

.. code-block:: python

   from vdjtools import sc

   cells = sc.read_10x("filtered_contig_annotations.csv")   # -> cell_id-keyed Rearrangement frame
   cells = sc.resolve_chains(cells)                         # pick the productive chain per locus
   paired = sc.pair_chains(cells, locus_pair="TRA_TRB")     # one row per α/β cell

   sc.paired_pgen(paired)                 # adds pgen_alpha, pgen_beta, pgen_paired (= product)

``write_airr_cell`` exports the AIRR Cell / Receptor format; ``to_anndata`` bridges to the
scverse ecosystem.

Recombination model
-------------------

:mod:`vdjtools.model` is the native V(D)J engine — generation probability, sampling, and EM
inference. Precomputed models for all 7 human loci ship in the wheel:

.. code-block:: python

   from vdjtools.model import load_bundled, native
   from vdjtools.model.generate import generate

   model = load_bundled("TRB", "olga")           # or "learned" (fit to real repertoires)
   native.pgen_aa(model, "CASSLAPGATNEKLFF")      # amino-acid Pgen (matches OLGA to 1e-15)
   native.pgen_aa(model, "CASSLAPGATNEKLFF", mismatches=1)     # + the Hamming-1 ball
   native.pgen_aa_batch(model, seqs, threads=0)   # many CDR3s, thread-parallel (~11x)
   generate(model, 1000)                          # sample a repertoire -> DataFrame

Learn a model from your own out-of-frame reads with :func:`vdjtools.model.infer.infer_native`,
and explore any model's recombination Bayes net interactively with
``marimo edit notebooks/model_explorer.py``. See the :doc:`API reference <api>` for the full
surface.

Command line
------------

Every workflow above has a CLI counterpart; inputs are auto-detected and results are written
as TSV to ``-o`` (or stdout):

.. code-block:: bash

   vdjtools models                                # list the bundled models
   vdjtools generate -m TRB -n 1000 -o gen.tsv    # (cf. olga-generate_sequences)
   vdjtools pgen seqs.tsv -m TRB -o pgen.tsv      # (cf. olga-compute_pgen)
   vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
   vdjtools overlap       *.tsv -o overlap.tsv
   vdjtools segment-usage *.tsv --segment v -o usage.tsv
   vdjtools spectratype   *.tsv -o spectra.tsv

Run ``vdjtools <command> --help`` for options.
