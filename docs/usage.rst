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

   import numpy as np
   import polars as pl
   from vdjtools import io as vio
   from vdjtools.io import schema as S
   from vdjtools.model import load_bundled
   from vdjtools.model.generate import generate

   def demo_sample(locus="TRB", clones=2000, reads=20_000, a=1.0, seed=0):
       """A canonical clonotype frame: sequences from a bundled model, clone sizes from a rank law."""
       seqs = generate(load_bundled(locus, "olga"), clones, seed=seed, productive_only=True)
       rng = np.random.default_rng(seed)
       freq = np.arange(1, seqs.height + 1, dtype=float) ** -a      # f_i proportional to i**-a
       counts = rng.multinomial(reads, freq / freq.sum())
       out = seqs.with_columns(pl.Series(S.COUNT, counts)).filter(pl.col(S.COUNT) > 0)
       return vio.normalize(out, recompute_freq=True)

   sample = demo_sample(seed=1)
   sample.columns
   # ['v_call', 'd_call', 'j_call', 'c_call', 'junction_aa', 'junction_nt',
   #  'duplicate_count', 'frequency']

.. important::

   The clone sizes are what make this usable as a demo, and an earlier version of this page got it
   wrong. Drawing ``n`` sequences and counting them gives ``n`` clonotypes of size 1, on which every
   diversity statistic is degenerate rather than merely synthetic: Shannon equals the richness,
   normalised Shannon is exactly ``1.0``, inverse Simpson equals the richness, and Chao1 diverges
   (8.0e6 from 4,000 singletons). Those are not numbers any repertoire produces. Clone sizes drawn
   from the rank law ``f_i ∝ i**-a`` give 1,790 clonotypes over 20,000 reads with 371 singletons and
   a 12.1 percent top clone, which is the regime the estimators are built for.

   Note also that ``rng.zipf`` is **not** the right draw here: it samples integers *from* a Zipf
   distribution, whose mean is infinite at ``a = 1.5``, so normalising the result hands one clone
   almost all the mass and collapses the repertoire to a handful of survivors. The rank law above is
   the intended ``f_i ∝ i**-a``.

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

Adaptive/immunoSEQ gene names are resolved through a shipped, CDR-validated lookup table
(``resources/adaptive_imgt_map.tsv``) rather than by stripping zero-padding: whether the
trailing group of ``TCRAV01-01`` is an IMGT subgroup (``TRAV1-1``) or an allele
(``TCRAV22-01`` → ``TRAV22``) is a per-family fact no pattern can decide. Family calls
(``TCRBV12-X``), slash ties (``TCRBV03-01/03-02``) and co-locus names (``TCRAV38-02`` →
``TRAV38-2/DV8``) resolve too; tokens outside the table fall back to the legacy rewrite.

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

.. note::

   All three regions are cut from ``junction_aa``, so ``region="all"`` **includes** the Cys104 and
   Phe118 anchors and is two residues longer than the IMGT :term:`CDR3`. ``region="trimmed"`` cuts
   three residues from each end --- the anchor plus the two templated residues beside it, not one
   --- and ``region="center"`` takes the middle five, which on a typical junction is the
   non-templated part. Pick the region for the chemistry you mean, and say which one you used.

Somatic hypermutation (B cells)
-------------------------------

A T cell's V gene is germline for life; a B cell's is rewritten in the germinal centre, and the
*pattern* of that rewriting is the read-out. :mod:`vdjtools.stats.shm` summarises it.

``v_identity`` is **not** a canonical column — most repertoire formats do not carry it — so ask the
reader to keep it:

.. code-block:: python

   from vdjtools.io.read import read_airr
   from vdjtools.stats import shm_summary, shm_spectrum

   df = read_airr("sample.tsv", keep=("v_identity",))
   shm_summary(df)      # flat dict: mutation load, switching, the germinal-centre marks
   shm_spectrum(df)     # the mutation-level distribution, 20 bins

Without ``keep=`` every SHM field is ``nan``, which is the point: a repertoire whose aligner never
reported identity is not an unmutated repertoire, and the two must not produce the same number.

The motivating application is **tertiary lymphoid structure** detection in tumours — a TLS is an
ectopic germinal centre. A GC leaves three joint marks and none is sufficient alone:

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - mark
     - field
     - what scores falsely without it
   * - mutated V genes
     - ``frac_mutated``, ``mean_shm``
     - naive B-cell infiltrate scores 0
   * - class switching
     - ``frac_switched``
     - an IgM-only infiltrate scores 0
   * - *active* diversification
     - ``switch_shm_gap``, ``shm_entropy``
     - a **resident plasma-cell clone** scores high on both of the above

That last row is why these are reported separately rather than folded into one index. A clone that
switched and hypermutated somewhere else arrives finished — high load, high switching, no ongoing
diversification — and a composite score cannot tell it from an active germinal centre. Which of the
three carries the signal is usually the interesting part, so the composite is yours to make.

``weight="freq"`` answers "what fraction of the *cells* are mutated"; ``weight="unique"`` answers
"what fraction of the *lineages*". In a repertoire with one dominant plasma-cell clone these differ
enormously, and which you want depends on the question.

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
   preprocess.filter_productive(sample)             # AIRR: drop out-of-frame / stop-codon
   preprocess.filter_functional_genes(sample)       # IMGT: drop pseudogene / ORF V and J
   preprocess.filter_length(sample)                 # junction_aa 5..60 aa, inclusive
   preprocess.filter_frequency(sample, min_freq=1e-4)
   preprocess.correct(sample, max_mismatches=2)     # collapse likely sequencing errors

   # cross-batch V/J-usage bias
   usage = preprocess.correct_vj_usage(cohort, batch_col="batch", transform="sigmoid")
   fixed = preprocess.apply_vj_correction(sample, usage, sample_id="A0")

   # combine samples: pooled clonotype table, or an incidence/frequency join
   preprocess.pool_samples([a, b, c])
   preprocess.join_samples([a, b, c], min_samples=2)   # clonotypes seen in >=2 samples

.. seealso::

   :doc:`preprocessing` — the whole pre-processing surface, including why *productive* (the
   rearrangement, AIRR) and *functional* (the germline gene, IMGT) are different questions, and
   how ``frequency`` is handled at every step.

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

Longitudinal dynamics
---------------------

:mod:`vdjtools.dynamics` compares two samples of the **same donor** — which clonotypes changed
between timepoints. Depth is handled *per pair* via a two-step effective sample size (never by
normalising a whole cohort to a common depth, which is not defined here):

.. code-block:: python

   from vdjtools import dynamics, preprocess

   day0 = demo_sample(seed=1)
   day15 = preprocess.downsample(day0, day0["duplicate_count"].sum() // 2, seed=7)   # a later draw

   # paired within-donor test: every clonotype -> emergent / expanded / persistent /
   # contracted / vanishing (or `untested` when the pair is too shallow to call it)
   tracked = dynamics.test_pair(day0, day15)
   tracked["dynamics"].value_counts()

   # group near-identical CDR3s (a 1-Hamming ball, V/J-pinned) BEFORE testing, for power on
   # convergent expansions — scope "1,0,0,1" = 1 substitution, "1,1,1,1" = 1 edit (Levenshtein)
   grouped = dynamics.test_metaclonotypes(day0, day15, scope="1,0,0,1")

   # complementary edgeR NB-exact caller (TMM + qCML common dispersion) -> log2FC + p
   called = dynamics.expansion_test(day0, day15)

The **recapture model** (Pavlova, Zvyagin & Shugay, *Front Immunol* 2024) asks a different
question: do *particular* clonotypes reappear more than the background rate, and does that depend
on clone size? Clonotypes are binned by size (singleton / doubleton / tripleton / large), each
bin's recapture fraction gets a ``Beta`` credible interval, and a log-linear model tests the
**group** effect:

.. code-block:: python

   import polars as pl

   # tag the clonotypes of interest (e.g. antigen-specific) vs the rest, in the PRE sample
   specific = set(day0["junction_aa"].head(20))
   pre = day0.with_columns(
       pl.when(pl.col("junction_aa").is_in(specific)).then(pl.lit("specific"))
         .otherwise(pl.lit("background")).alias("group"))

   rates = dynamics.capture_rates(pre, day15, group_col="group")   # recapture per (group, size)
   dynamics.capture_test(rates)                                    # log-linear group effect + p

Across several donors, concatenate the per-donor ``capture_rates(..., donor=<id>)`` frames and use
:func:`~vdjtools.dynamics.capture_paired_test` for the per-size-class paired t-test. The full
worked example — sunken/alluvial tracking plots, a VDJdb overlay, and the capture ribbons — is
``marimo edit examples/vaccination_tracking.py``.

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

The **match scope** is set by ``key`` (``(junction_aa,)`` / ``+v`` / ``+v+j``) and ``match``:

``"exact"``
   The feature is the key itself.

``"fuzzy"``
   A single-mismatch **search**: ``incidence(c) = #subjects carrying ANY feature within `scope`
   of c``. The candidate **keeps its identity and gains incidence** — the mismatch is an
   incidence-*estimation* trick, so the output stays a list of individual biomarker clonotypes.
   Non-``junction_aa`` key columns (V/J) must match exactly, so ``key=(junction_aa, v_call)``
   pins the germline half of the contact surface while the CDR3 varies by one residue. This is
   the Vlasova-2026 operation, and it delegates the search to ``vdjmatch.cluster.overlap``.

``"1mm"``
   **Clustering** via :func:`~vdjtools.biomarker.metaclonotypes`: candidates are *merged* into
   groups and the group is tested. A legitimate operation, but a different one — it belongs
   *downstream* of a biomarker list (building a Hamming graph, a classifier), not to discovery.
   Use ``"fuzzy"`` to find biomarkers; use ``"1mm"`` to group them afterwards.

Candidate features are all public clonotypes (``min_incidence`` count or ``min_incidence_frac``
fraction) or an explicit ``candidates`` list (:func:`~vdjtools.biomarker.select_candidates`).
Note that under ``match="fuzzy"`` ``candidates`` is the **query** set only: the search universe
stays the whole cohort, because a candidate's neighbours are usually not candidates themselves.
``fisher_association`` is kept as the Emerson-2017 Fisher shortcut with the original schema.

.. admonition:: Choosing the unit — and, if you count rearrangements, the right null
   :class: important

   The sampling unit is the **subject**. Emerson 2017 tested template-weighted abundance
   head-to-head against presence/absence and abundance *lost*; weighting a contingency table by
   reads is pseudoreplication (Hurlbert 1984). ``association`` therefore tests subject incidence.

   If you instead count **rearrangements** (a unique nt row in a repertoire = one recombination
   event; a CDR3aa reached by three nt variants is three events), the counts are large (~10⁷ in a
   1,200-donor cohort), so an exact hypergeometric is the wrong tool — use a smooth test
   (conditional binomial / G-test), not factorials. **The null matters more than the test**: it
   must be the *subject* ratio ``n_pos/(n_pos+n_neg)``, not the *row* ratio
   ``N_pos/(N_pos+N_neg)``. Where sequencing depth differs between arms — in our FMBA cohort the
   controls carry 1.4–1.5× more rearrangements per donor — the two nulls differ by ~15–20%, and
   any clonotype whose row count does not scale with depth (a public one carried about once per
   donor regardless) picks up a spurious enrichment of exactly that size. At large counts a 1.4×
   bias is hyper-significant, so the wrong null does not add noise, it manufactures hits.

   Depth differences also bias the subject test, in the opposite direction: a deeper repertoire
   is more likely to carry any given clonotype. Where the design allows it — repeated samples of
   the same donor — down-sample each pair to a common read count
   (:func:`~vdjtools.preprocess.downsample`) before testing rather than trying to model it away.

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
``marimo edit examples/biomarker_explorer.py`` (Emerson HIP via HuggingFace).

The longitudinal and cohort workflows have their own interactive notebooks (``pip install
"vdjtools[examples]"``): ``marimo edit examples/vaccination_tracking.py`` (clonotype tracking +
the :mod:`vdjtools.dynamics` recapture model across yellow-fever / influenza / TBE vaccination
time courses), ``examples/aging.py`` (cohort-streaming diversity / clone-size / spectratype
across the Britanova "Cord Blood to Centenarians" cohort), and ``examples/ankspond_motif.py``
(the ankylosing-spondylitis TRBV9 "AS27" motif — disease vs HLA-B27 carriage). Each prefers a
local ``./data_dump/`` copy (gitignored; symlink your data there), else fetches from HuggingFace.

Single-cell
-----------

:mod:`vdjtools.sc` ingests 10x / AIRR-Cell / arda contigs into a flat ``cell_id``-keyed
frame, resolves and pairs chains with doublet / mispairing QC, and scores paired α/β
generation probability:

.. code-block:: python

   from vdjtools import sc

   cells = sc.read_airr_cell("outs/airr_rearrangement.tsv")  # -> cell_id-keyed frame
   cells = sc.resolve_chains(cells)                          # pick the productive chain per locus
   paired = sc.pair_chains(cells, locus_pair="TRA_TRB")      # one row per α/β cell

   sc.paired_pgen(paired)                 # adds pgen_alpha, pgen_beta, pgen_paired (= product)

   adata = sc.to_scirpy(cells)            # scverse-native AnnData (obsm["airr"])
   vdj   = sc.to_dandelion(cells)         # dandelion
   sc.write_screpertoire(cells, "airr_rearrangement.tsv")    # scRepertoire (R)

See :doc:`singlecell` for the full ingestion matrix, the interop round-trips, and the
``vdjtools sc`` command line.

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
   native.pgen_aa_batch(model, seqs, threads=0)   # many sequences, thread-parallel (~11x)
   generate(model, 1000)                          # sample a repertoire -> DataFrame

.. important::

   **These functions take a junction, not an IMGT CDR3**, even though the parameter is spelled
   ``cdr3_aa`` for historical reasons. The sequence must start with the conserved Cys104 and end
   with Phe118 (Trp118 for IGH) --- ``CASSLAPGATNEKLFF``, not ``ASSLAPGATNEKL``. The same applies
   to :func:`~vdjtools.model.infer.infer_nt`. Passing a trimmed CDR3 does not raise; it returns the
   Pgen of a different, shorter rearrangement. See :term:`junction`.

Pgen of a motif
^^^^^^^^^^^^^^^

:func:`~vdjtools.model.native.pgen_aa_degenerate` scores a **per-position set of permitted
residues** instead of one sequence: the total Pgen of every junction the motif matches. That is
the model-side estimate of how often a paratope motif — a V/J/length-pinned VDJdb cluster PWM
thresholded per position, or any wildcard/gapped pattern — is generated at all.

.. code-block:: python

   motif = list("CASSLAPGATNEKLFF")
   motif[4] = "ILV"        # position 4: any of I, L, V
   motif[5] = "X"          # position 5: any residue  ("" means the same)
   native.pgen_aa_degenerate(model, motif, v="TRBV5-1*01", j="TRBJ2-3*01")
   native.pgen_aa_degenerate_batch(model, [motif, ...], threads=0)   # thread-parallel

It is one transfer-matrix pass whatever the sets contain — the matching sequences are summed over,
never enumerated — and it is the DP that has always driven ``pgen_aa`` (every position pinned) and
``mismatches=1`` (one wildcard position at a time). Pinning every position reproduces
:func:`~vdjtools.model.native.pgen_aa` bitwise; wildcarding every position gives the model's length
marginal ``P(L, V, J)``.

.. warning::

   ``X`` is a wildcard **here only**. :func:`~vdjtools.model.native.pgen_aa` matches residues by
   exact character against the genetic code, so ``pgen_aa(model, "CASSLAPGATNEKLXF")`` has no
   codons at that position and returns ``0.0`` — a silent zero, not a degenerate query. Pass
   ``list(seq)`` to ``pgen_aa_degenerate`` when the sequence carries ``X``.

Reconstructing nucleotides and V/D/J markup from an amino-acid CDR3
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A VDJdb record carries ``(V, J, CDR3aa)`` and no nucleotides, so the boundary markup a repertoire
analysis wants is simply absent. :func:`~vdjtools.model.viterbi.infer_nt` reconstructs it under the
model: germline positions are pinned to their V/D/J segment, and each free N-region position takes
the nucleotide maximising ``P(nt_1) * prod P(nt_k | nt_{k-1})`` under the insertion model covering
it.

.. code-block:: python

   from vdjtools.model import infer_nt

   sc = infer_nt(model, "CASSLGQAYEQYF", v="TRBV5-1*01", j="TRBJ2-3*01")
   sc.cdr3_nt                       # the reconstructed nucleotide CDR3
   sc.v_end, sc.d_call, sc.d_start, sc.d_end, sc.j_start   # 0-based, half-open, CDR3-nt space
   sc.pgen, sc.margin               # exact Pgen of the result, and its lead over the runner-up

``v=``/``j=`` accept one allele, several (a list, or the comma-separated string an ambiguous AIRR
``v_call`` carries), or ``None`` — the search then marginalizes over every gene, which costs almost
nothing because the underlying DP sweeps V and J anyway. Reported timings: 2.5 ms per human TRB
CDR3, 0.5 ms per TRA, so all 80k VDJdb records take about three minutes.

.. warning::

   Report ``sc.margin``. With a long non-templated core many nucleotide sequences are near-equally
   likely, and a bare "most likely" claim is then close to meaningless — the margin is what says so.

.. note::

   Pass the :class:`~vdjtools.model.model.Model`, not a :func:`~vdjtools.model.pgen.prepare` -d one.
   A prepared model selects the pure-Python reference search, which is the implementation the native
   one is validated against and is roughly 600x slower on a VDJ locus.

Every plausible rearrangement, with its probability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:func:`~vdjtools.model.native.best_aa_scenarios` is the first stage alone: the top ``k`` scenarios
the same transfer matrix sums over, without reconstructing nucleotides. Use it for the three things
an aligner cannot produce — **the alternatives and their weights**, so boundary ambiguity is visible
instead of hidden behind one committed answer; **D geometry** (``idx5``, ``idx3``, ``pos``); and
**naming a missing V or J** by marginalising over the locus when the call is absent.

.. code-block:: python

   from vdjtools.model import best_aa_scenarios, best_aa_scenarios_batch

   for w, v, len_v, j, len_j, d, idx5, idx3, pos in best_aa_scenarios(
           model, "CASSLGQAYEQYF", v="TRBV9", j="TRBJ2-7", k=8):
       ...                          # v_end is len_v; j_start is 3*len(aa) - len_j

   df = best_aa_scenarios_batch(model, frame["junction_aa"].to_list(),
                                v=frame["v_call"].to_list(), j=frame["j_call"].to_list(),
                                k=8, threads=0)   # one row per scenario

A **gene**-level name (``TRBV4-3``, which is what CellRanger and most real repertoires carry)
resolves to that gene's representative allele; a name the model has no gene for raises and names it.
So an empty list means one thing only: the DP explains nothing, because a residue has no codon the
model can reach. The frame entry point returns ``row``, ``rank``, ``w``, ``v_call``, ``len_v``,
``j_call``, ``len_j``, ``d_call``, ``idx5``, ``idx3``, ``pos``, gives values identical to the
per-row loop at any thread count, and contributes no rows for a query it declines.

.. warning::

   This is the argmax of a probability model, not an aligner, and it is **worse than germline
   alignment at placing a boundary you already have the calls for**. Measured on 25,000 real human
   TRB clonotypes against the observed nucleotide markup, each method seeing only the amino acids
   and the V/J calls:

   .. list-table::
      :header-rows: 1

      * - method
        - ``v_end`` exact
        - ``j_start`` exact
        - declined
      * - germline alignment (:mod:`arda.cdr3fix`)
        - **73.3%**
        - **97.7%**
        - **0**
      * - ``best_aa_scenarios``, top-1 of ``k=8``
        - 59.8%
        - 90.3%
        - 642
      * - ``infer_nt``
        - 60.1%
        - 90.3%
        - 642

   Many junctions have several near-equally-probable boundaries, so the modal scenario is often not
   the one that happened. Reach for the scenario list when you want the distribution; reach for
   alignment when you want the boundary. The 642 declines above predate the gene-level resolution
   described here.

Building a model on your own germline library, fitting it to your own reads, checking it, comparing
two of them and asking how much diversity one describes are covered in
:doc:`Recombination model workshop <model>`. Explore any model's Bayes net interactively with
``marimo edit examples/model_explorer.py``, or work through the whole workshop with
``marimo edit examples/model_workshop.py``. See the :doc:`API reference <api>` for the full surface.

Command line
------------

Every workflow above has a CLI counterpart; inputs are auto-detected and results are written to
``-o`` — TSV, or **Parquet when the path ends in** ``.parquet`` / ``.pq`` — or to stdout:

.. code-block:: bash

   vdjtools models                                # list the bundled models
   vdjtools generate -m TRB -n 1000 -o gen.tsv    # (cf. olga-generate_sequences)
   vdjtools pgen seqs.tsv -m TRB -o pgen.tsv      # (cf. olga-compute_pgen)

   # data: convert any format to the canonical table, and preprocess
   vdjtools convert mixcr.txt.gz -o clones.parquet    # → canonical Parquet (or .tsv)
   vdjtools filter clones.parquet --coding --min-freq 1e-4 -o coding.tsv
   vdjtools downsample clones.parquet 100000 -o ds.tsv
   vdjtools pool s1.tsv s2.tsv --join --min-samples 2 -o joint.tsv

   # analytics (sample files, or -m metadata + --base-dir; parallel with -t, or stream --cohort)
   vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
   vdjtools overlap       *.tsv -o overlap.tsv
   vdjtools segment-usage *.tsv --segment v -o usage.tsv
   vdjtools spectratype   -m metadata.txt --base-dir samples/ -t 8 -o spectra.tsv
   vdjtools diversity     --cohort cohort_parquet/ -o div.tsv

   # longitudinal + enrichment
   vdjtools dynamics day0.tsv day15.tsv -o tracked.tsv    # paired within-donor expansion test
   vdjtools tcrnet   sample.tsv -o net.tsv                 # neighbourhood enrichment (control cohort)

Run ``vdjtools <command> --help`` for options.
