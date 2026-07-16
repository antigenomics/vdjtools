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
auto-detects and converts the common third-party formats (MiXcr v1–4 incl. the C-gene / BCR
isotype, MiGec, Adaptive immunoSEQ v1/v2, IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4, and arda's
AIRR annotation output). Every reader returns the same canonical frame:

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

Cross-batch V/J-usage bias is corrected in two steps — batch-correct the usage, then push it
back onto each sample's clonotype table:

.. code-block:: python

   # cohort = one long frame with sample_id + batch columns (see io.read_samples)
   usage = preprocess.correct_vj_usage(cohort, batch_col="batch", transform="sigmoid")
   fixed = preprocess.apply_vj_correction(sampleA, usage, sample_id="A0")   # resampled table

``transform="location"`` (default) is a ComBat location adjustment; ``transform="sigmoid"`` is the
σ-standardised, grand-mean-preserving correction of Vlasova et al. 2026 (*Genome Medicine* 18:20).
``apply_vj_correction`` roulette-wheel resamples the clonotype table to the corrected usage
(``resample=False`` for deterministic expected counts). The batch mean/σ use the plain log-normal
statistics by default (paper-faithful); pass ``winsor_q=0.025`` only for the noisy
usage-as-features regime (many shallow / RNA-seq repertoires, e.g. as UMAP features). Validated on
the paper's own FMBA covid TCRβ cohort (deep repertoires, ~3.3M reads/sample): V-usage variance
explained by batch drops from η²≈0.11 to ≈0.002 while the grand-mean usage and per-sample read
depth are preserved.

Biomarker association
---------------------

:mod:`vdjtools.biomarker` tests each clonotype **feature**'s incidence (presence across
subjects) against a condition — the incidence-contingency framework of Emerson 2017, Howie
2015, De Witt 2018 and Vlasova 2026. ``association`` shares one streamed subject-incidence
table across five tests (Fisher, χ², Bayesian log-odds, Beta-Binomial Bayes factor,
permutation) and three condition types (binary, category, stratified):

.. code-block:: python

   import polars as pl
   from vdjtools import biomarker
   from vdjtools.biomarker import association, condition

   # cohort: long frame with a sample_id column; meta: one row per subject
   cohort = pl.DataFrame({
       "sample_id":   ["p0","p1","p2","n0","n1","n2"],
       "v_call":      ["TRBV1"]*6, "j_call": ["TRBJ1"]*6,
       "junction_aa": ["CASSXF","CASSXF","CASSXF","CASSXF","CASSBG","CASSBG"],
       "duplicate_count": [10]*6,
   })
   meta = pl.DataFrame({"sample_id": ["p0","p1","p2","n0","n1","n2"],
                        "cmv": ["+","+","+","-","-","-"], "hla": ["A*02"]*3 + ["A*01"]*3})

   # binary condition, several tests at once (long output with a `test` column)
   association(cohort, condition.binary(meta, "cmv"),
               test=["fisher", "chi2", "bayes_bf"])

   # category: one-vs-rest per HLA allele (Emerson/DeWitt); add a `level` column
   association(cohort, condition.hla_alleles(meta, ["hla"]), level_col="_level")

   # paired: CMV association conditioned on HLA, combined by Cochran–Mantel–Haenszel
   association(cohort, condition.stratified(meta, "cmv", "hla"), stratum_col="_stratum")

The **match scope** is set by ``key`` (``(junction_aa,)`` / ``+v`` / ``+v+j``) and ``match``
(``"exact"`` or ``"1mm"``, the latter pooling single-mismatch neighbours via ``metaclonotypes``,
needs the ``[overlap]`` engine). Candidate features are all public clonotypes
(``min_incidence`` count or ``min_incidence_frac`` fraction) or an explicit ``candidates`` list
(:func:`~vdjtools.biomarker.select_candidates`). ``fisher_association`` is kept as the
Emerson-2017 Fisher shortcut with the original column schema.

**Co-occurrence** tests feature-vs-feature incidence across the subjects profiled for both
chains — in-silico α-β pairing (Howie 2015, Vlasova 2026) and same-chain co-specificity
(De Witt 2018) — by the lift ``θ = n·n_AB/(n_A·n_B)`` + Fisher/χ² + FDR:

.. code-block:: python

   biomarker.cooccurrence(cohort, chain_a="TRA", chain_b="TRB", evalue=True)   # α-β pairs
   biomarker.cooccurrence(cohort, chain_a="TRB", chain_b=None)                  # same-chain pairs

   biomarker.metaclonotypes(cohort)       # group near-identical CDR3s (1-mismatch + V/J) -> meta_id

.. warning::

   **Cross-subject co-occurrence is not evidence of physical chain pairing.**
   :func:`~vdjtools.biomarker.cooccurrence` tests whether two clonotypes occur in the same
   *subjects* more often than chance. Unlike the randomised wells of a pairSEQ experiment
   (Howie et al., *Sci Transl Med* 2015), subjects carry HLA type, germline variation,
   ancestry, sequencing depth and infection history — any of which makes two clonotypes
   co-occur without ever sharing a cell. In pairSEQ's own framework a *cross-subject* α-β
   pair is the **definition of a false positive**. Two confounders in particular survive HLA
   stratification:

   - **Repertoire depth** inflates the lift by ``≈ 1+CV²(N)`` for rare clonotypes,
     independently of HLA, exposure or pairing (measured on the FMBA covid19 cohort:
     CV(N)=0.899 → θ_depth=1.81; a ≥1000-clonotype floor drops it to 1.50 *and* removes
     ~73% of the significant pairs). Filter near-empty samples before reading θ.
   - **Shared exposure**: two co-specific but unpaired clones stay associated within every
     stratum — indistinguishable from pairing by any cross-subject contingency table.

   Shared restriction by an allele of carrier frequency *f* cannot induce a lift above
   ``(1+CV²(N))/f`` (2.04 for HLA-A*02:01 before the depth term), so a θ far above that
   ceiling is not explicable by *that allele* — though it remains explicable by depth,
   ancestry, batch or exposure. De Witt et al. (*eLife* 2018) found shared HLA carriage
   explained the **majority** of strongly co-occurring TCRβ pairs, and restricted all
   downstream clustering to within-allele subsets. Read ``cooccurrence()`` output as
   co-occurrence and nothing more; physical pairing is established by within-subject
   designs (well-based subsampling or single-cell), not cross-subject tables.

Explore the whole screen interactively — condition (CMV / HLA-allele / CMH), test, match
scope, a live VDJdb overlay, and a co-occurrence panel — with
``marimo edit notebooks/biomarker_explorer.py`` (Emerson HIP via HuggingFace).

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
