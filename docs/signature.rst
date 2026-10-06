Repertoire signatures
=====================

A *signature* is a fixed-order, name-addressed feature vector for one repertoire, on a scale a
downstream model can consume without fitting a scaler of its own. ``vdjtools`` emits the
**statistics** half, ``vsig``; the **geometry** half, ``rsig``, comes from `mirpy
<https://github.com/antigenomics/mirpy>`_ and shares this contract exactly.

It is built in three stages, and a **corpus** --- a large published reference collection of
repertoires --- fixes the last two:

.. code-block:: text

   raw features  ->  winsorization bounds  ->  rotation + per-PC scaling  ->  vsig
   (this sample                (all three from one corpus artifact)
    alone)                                                              +  channels, untouched

In plain terms, reading left to right:

#. **Raw features** are the ordinary measurements of one repertoire --- how diverse, how clonal,
   which V and J genes, what junction lengths, which isotypes. They come from this sample and
   nothing else.
#. :term:`Winsorization` clips each measurement to a range measured on the corpus, so that one
   pathological sample cannot dominate everything downstream. Clipping edits a real number, so how
   much of the row was clipped is always reported --- see ``winsor_frac`` below.
#. The :term:`rotation` re-expresses those correlated measurements as uncorrelated axes
   (:term:`principal component` s) ordered by how much variation each carries across the corpus,
   and each axis is put on a :term:`robust z-score` scale: how far this sample sits from the
   corpus's typical one, in robust deviations. That is what makes your matrix and a collaborator's
   directly comparable without either of you fitting a scaler.

:term:`Channel` s sit outside all of this. They are carried in their own units, never rotated and
never clipped, because a number whose job is to tell you whether to trust a row must not be mixed
into the row.

If you want the underlying measurements back rather than the rotated axes, ask for them:
``--named`` emits the raw groups the rotation is fitted on under their own names. See
:ref:`sig-named`.

Quickstart
----------

.. code-block:: bash

   # nine corpora are published; name one and it is fetched on first use
   vdjtools signature --corpus blood samples/*.tsv.gz -o sig.tsv

   # what exactly will I get?
   vdjtools signature --corpus blood --components 32 --describe

   # pre-warm the cache instead of fetching lazily -- the install-time step
   vdjtools corpus --fetch all

   # or build your own synthetic corpus, which needs no cohort at all
   vdjtools corpus --corpus synthetic-tissue -o /tmp/st.npz
   vdjtools signature --corpus /tmp/st.npz samples/*.tsv.gz -o sig.tsv

.. code-block:: python

   from vdjtools.signature import vsig, vsig_cohort
   from vdjtools.signature.corpus import Corpus, bundled_path

   corpus = Corpus.load(bundled_path("synthetic-blood"))
   row = vsig({"TRB": trb, "IGH": igh}, corpus)
   frame = vsig_cohort({"S1": sample1, "S2": sample2}, corpus, n_jobs=0)

Pick ``synthetic-blood`` or ``synthetic-tissue`` by which compartment your samples come from; both are
described below, with what they reproduce and what they do not.

.. _sig-choosing:

Choosing a corpus
-----------------

This is the first decision and the one that matters most, so start here. Match the corpus to how your
samples were produced, not to what you are testing for.

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - Your samples are
     - Use
     - Why
   * - Bulk RNA-seq of blood or PBMC
     - ``blood``
     - Fitted on 11,117 real blood samples across 947 study groups. The default choice for
       peripheral-blood work.
   * - Bulk RNA-seq of a solid tissue or tumour
     - ``tissue``
     - Fitted on 21,131 real non-blood samples across 1,934 study groups. Tissue repertoires are
       shallower and more skewed than blood, and a blood corpus clamps them harder.
   * - Targeted or amplicon deep TCR sequencing
     - ``deep-tcr``
     - Fitted on 3,936 amplicon samples. **TRA and TRB only** -- it has no other locus to offer.
   * - Any of the above, and you need the per-study cap's effect
     - ``blood-uncapped``, ``tissue-uncapped``
     - The same populations with no cap at 30 samples per study, so the difference is measurable
       rather than assumed.
   * - Something none of these describe
     - ``synthetic-blood``, ``synthetic-tissue``
     - Drawn across the measured depth and clone-size range of that compartment, using no cohort
       sample at all -- so they rebuild bit-for-bit anywhere and carry no population you have to
       accept.
   * - A reference point rather than a cohort
     - ``naive``, ``memory``
     - The pure regimes: what the recombination model emits, and the same under a Zipf
       rank-abundance clone-size law.

Two practical consequences:

**A real corpus fits real data better, and by a measurable margin.** On four TCGA tissue samples at
``--components 32``, the fraction of the row that the bounds clamped was 0.199 through ``tissue``
against 0.595 through ``synthetic-tissue`` -- roughly three times less clamping. ``deep-tcr`` on the
same samples clamps 0.483 with a median absolute z of 2.1, correctly far out: amplicon TCR is the
wrong reference for RNA-seq, and the numbers say so rather than quietly absorbing it.

**Check ``winsor_frac`` before you trust a matrix.** ``vsig:qc:-:winsor_frac`` is the fraction of the
row the corpus's bounds clamped. A high value means your samples sit outside the corpus's range, so
the standardisation is extrapolating -- pick a closer corpus rather than proceeding.

A corpus is required
--------------------

There is no default. A signature is comparable to another one **only** if both were rotated
through the same corpus, so the artifact has to be named, and the emitted row records which one it
was. Two matrices whose column names match but whose corpora differ are not comparable, and nothing
about the numbers says so — which is why the choice is not allowed to be implicit.

All nine corpora are published. Four are **synthetic** -- every receptor drawn from the bundled
recombination models, so anyone can rebuild the artifact bit-for-bit and no cohort sample is needed to
use one -- and five are **real**, fitted on repertoires. Both kinds exist on purpose: a synthetic
corpus spans a *measured* range by construction, while a real one carries the joint structure the
generative model does not produce (selection-shaped V/J usage, isotype and SHM structure, and the
cross-locus covariance of libraries prepared together).

.. list-table::
   :header-rows: 1
   :widths: 20 10 14 56

   * - corpus
     - kind
     - fitted on
     - what it is
   * - ``naive``
     - synthetic
     - 10,000 drawn
     - every clone size 1 -- what the recombination model emits
   * - ``memory``
     - synthetic
     - 10,000 drawn
     - Zipf rank-abundance, sampled by multinomial, zeros dropped
   * - ``synthetic-blood``
     - synthetic
     - 10,000 drawn
     - a naive/memory mixture across measured blood ladders (TRB richness 74-3,162, ``f1``
       0.215-0.949)
   * - ``synthetic-tissue``
     - synthetic
     - 10,000 drawn
     - the same across measured tissue ladders (TRB richness 30-2,977, IGH 40-10,352)
   * - ``blood``
     - real
     - 11,117 samples
     - public bulk RNA-seq blood, 947 study groups, capped at 30 per study
   * - ``blood-uncapped``
     - real
     - 22,441 samples
     - the same population with no cap, so the cap's effect is measurable
   * - ``tissue``
     - real
     - 21,131 samples
     - public bulk RNA-seq non-blood, 1,934 study groups, capped at 30 per study
   * - ``tissue-uncapped``
     - real
     - 33,874 samples
     - the same population with no cap
   * - ``deep-tcr``
     - real
     - 3,936 samples
     - targeted/amplicon deep TCR, **TRA and TRB only**

.. _sig-fetch:

Artifacts are fetched, not bundled
----------------------------------

At 256 components a ``vsig`` artifact is about 10 MB of float32 rotation, and the nine corpora across
both halves are roughly 110 MB. Making every ``pip install`` carry all of that in order to use one of
them is the wrong trade, so the wheel ships a few KB of **index** -- ``corpora.json``, naming every
corpus with its size and SHA-256 -- and the artifact itself is a GitHub release asset fetched on first
use into ``$VDJTOOLS_CORPUS_DIR`` (default ``~/.cache/vdjtools/signature``) and verified against that
digest.

The resolution order is **a local path, then the wheel, then the cache, then the release**. A local
path always wins, so a caller who fitted their own corpus and passes its path is never served a
download of the same name. A download whose digest does not match is deleted rather than cached: a
half-written rotation that loads is worse than one that is missing.

``vdjtools corpus --fetch all`` (and ``mir corpus --fetch all``) pre-warms the cache -- the
install-time step, since anything not pre-warmed is fetched lazily. The two halves are published by
the library that owns each and share one cache directory, because a ``vsig``/``rsig`` pair belongs in
one place.

**Which one?** Use ``blood`` or ``tissue`` when your samples are bulk RNA-seq of that compartment and
``deep-tcr`` for amplicon TCR: a real corpus is the closest reference to real data. Use
``synthetic-blood`` / ``synthetic-tissue`` when you need a reference that is reproducible from the
library alone, with no cohort behind it, or when your depths sit outside what the real corpora cover.
``naive`` and ``memory`` are the two pure regimes and neither is a bulk sample -- ``naive`` has no
clone-size structure at all, and ``memory`` at a fixed nominal size has a read count of exactly
``size × reads_per_clone`` in *every* sample -- so they are reference points for what selection does to
a rotation, not descriptions of a cohort.

.. _sig-cohort-corpus:

What comes out
--------------

Column names are four colon-separated parts throughout, so one :func:`~vdjtools.signature.parse`
reads every kind::

   vsig:pc:TRB:PC01          a rotated coordinate
   vsig:cov:TRB:cstar        a channel, carried through untouched
   vsig:qc:-:winsor_frac     a channel that is not per-locus

**Channels are never rotated and never clamped.** A provenance number that has been mixed with the
measurements it was supposed to qualify is no longer provenance. Read these three before trusting
any rotated column of a row:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - channel
     - what it tells you
   * - ``vsig:cov:<locus>:cstar``
     - The coverage this sample actually attained, in ``[0, 1]``. **Always** emitted, including
       when the coverage-standardised diversity features turn out not to be estimable — which is
       the case it matters most in.
   * - ``vsig:mask:<locus>:present`` / ``:estimable``
     - Why a column is a hole. "This locus is absent" and "this sample is too shallow to reach the
       target" are different facts that would otherwise both render as one ``nan``.
   * - ``vsig:qc:-:winsor_frac``
     - What fraction of this row the corpus's bounds clamped. See :ref:`sig-winsor`.

``vsig:qc:<locus>:v_fallback_frac`` and ``:j_fallback_frac`` are the other two worth reading: an
unrecognised V or J call raises nowhere, because the germline lookup falls back to the maximum
observed distance, so a cohort on Adaptive nomenclature or an older IMGT release yields a fully
populated, entirely plausible, systematically wrong result. These fractions are the only thing that
says so.

A hole is ``nan``, never ``0``.

.. _sig-winsor:

Winsorization: by percentile, and one-sided where the metric is
---------------------------------------------------------------

Bounds are estimated **once, on the corpus**, at a percentile — not as a multiple of a robust
standard deviation. A robust-SD bound is not a fixed amount of probability: it depends on the
distribution's shape, so an 8-SD bound trims nothing at all on a right-skewed count and a third of
the corpus on a near-degenerate column.

The side that gets trimmed comes from each feature's declared ``support``, and **never** from its
transform:

.. list-table::
   :header-rows: 1
   :widths: 16 20 64

   * - support
     - trimmed
     - because
   * - ``nonneg``
     - top only
     - ``[0, inf)`` — cannot run away downward. Counts, richness, Hill numbers, standard deviations.
   * - ``nonpos``
     - bottom only
     - ``(-inf, 0]`` — a log-probability is bounded above by 0.
   * - ``real``
     - both
     - ``(-inf, inf)`` — log-ratios, clr and logit coordinates, PC scores.
   * - ``unit``
     - neither
     - ``[a, b]`` — already bounded on both sides. Raw proportions, ``cstar``.

Three raw features declare ``transform="none"`` and have three different supports. Deriving the side
from the transform silently trims the wrong end of two of them, which is why ``support`` is a
declared field rather than an inference.

Two modes, and the mode fixes how the rotation was fitted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 14 38 48

   * - ``--winsorize``
     - at emit time
     - the corpus PCA was fitted by
   * - ``features``
     - clamp raw features to the corpus bounds, **then** rotate
     - winsorized corpus, median/MAD, plain SVD. Deterministic and bit-reproducible; the robustness
       sits in a step you can inspect.
   * - ``pcs``
     - rotate, **then** clamp the PC scores
     - ``RobustScaler`` and a full-solver PCA — there the corpus still has its tails at fit time,
       so the estimator has to resist them itself.
   * - ``none``
     - clamp nothing
     - whichever of the two the artifact records

Both bound sets and both percentiles (0.01 and 0.05) are stored, so ``--winsorize`` and
``--winsor-p`` are apply-time choices that need no refit.

**The clamp is never silent.** Clamping is many-to-one: it overwrites a real measurement of a real
sample, and the sample cannot tell afterwards. In the system this replaces, a reference whose centre
was wrong by 57 robust deviations went unnoticed because every affected sample was clamped to the
bound and looked like an ordinary tail case. So every row carries ``vsig:qc:-:winsor_frac``, and a
value near 1.0 means this sample does not belong to this corpus — not that it is unusual.

.. _sig-components:

How many components
-------------------

``--components`` / ``n_components=`` takes either kind of answer:

* an **integer** — a fixed count, identical at every locus. At fit time it is a cap: a block
  narrower than the request (the five cross-locus ratios against a request of 256) gets its full
  rotation, recorded per locus in the artifact.
* a **fraction** in ``(0, 1)`` — cumulative variance explained. The count then differs per locus
  *and per corpus*, because some corpora are genuinely easier to explain. That is biology, not a
  defect.

The artifact stores the rotation only up to the count it was fitted at, plus the **full eigenvalue
spectrum**. So truncating downward at apply time is free and **exact** — the components are ordered,
so the first *n* of a longer rotation are the same vectors. A **count** above the stored one is
capped per locus, exactly as the fit capped it; a **fraction** the stored spectrum cannot reach
raises, and the error quotes what it does reach. Neither ever pads: no column appears under a name
the rotation does not hold.

.. _sig-component-table:

The coverage target is a runtime argument
-----------------------------------------

Hill numbers are read at a coverage level, and **that level is not in the artifact**. A per-sample
quantity has no business in a corpus artifact: in the system this replaces it lived there, and one
shipped reference ended up standardising tissue samples to a level measured on blood, identical to
17 significant digits across all seven loci while its location and scale had genuinely been refitted.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - ``--cstar-target``
     - meaning
   * - ``min`` *(default)*
     - This cohort's own per-locus minimum attained coverage. Comparable across its samples, and a
       property of *this* cohort rather than of somebody's corpus. Costs a cheap extra pass over the
       count column — and resolves deferred samples twice, so pass a number for a one-pass run.
   * - a number
     - That level for every locus. One pass.
   * - ``own``
     - Each sample at its own attained coverage: the diversity features are then *observed* rather
       than standardised, and are not comparable across samples.

Whatever the target, ``vsig:cov:<locus>:cstar`` reports what the sample reached and
``vsig:mask:<locus>:estimable`` reports whether the target was reachable without extrapolating.

Raw feature groups
------------------

Every group at a locus is concatenated into **one** feature vector and rotated together, so the PCA
can find cross-group structure. One rotation per ``(half, locus)``, plus one for the cross-locus
group — *not* one joint rotation over all loci, because under a joint rotation one absent locus
silently moves every component. Per-locus blocking is what keeps a hole a hole.

.. list-table::
   :header-rows: 1
   :widths: 12 18 70

   * - group
     - width
     - what
   * - ``div``
     - 6
     - Coverage-standardised Hill numbers ``0D_c``/``1D_c``/``2D_c``, ``clonality``, ``0D_chao``,
       ``d50``
   * - ``depth``
     - 3
     - ``reads``, ``richness``, ``S_unseen`` (Chao's estimate of what was never drawn)
   * - ``clon``
     - 3
     - Clone-size composition ``f1``/``f2`` in clr coordinates, plus the top clone's share
   * - ``len``
     - 3
     - Weighted junction-length ``mean``, ``sd``, ``skew``
   * - ``vus`` / ``jus``
     - germline
     - V and J gene usage over the arda germline vocabulary, clr
   * - ``spec``
     - germline x 26
     - Spectratype: junction length resolved per V gene, clr
   * - ``kmer``
     - 400
     - Junction 2-mer composition, arcsine
   * - ``aa``
     - 20
     - Single-residue composition, arcsine
   * - ``pchem``
     - 30
     - Physicochemistry over two junction regions
   * - ``iso`` / ``shm``
     - 5 / 1
     - IGH only: isotype composition in clr, and mean V identity
   * - ``pair``
     - 5
     - Cross-locus: five log read-count ratios, depth divided out

Widths per locus run from 685 (TRD) to 3,744 (IGH), 13,483 raw features in total.

.. note::

   ``spec`` is by far the widest group and the sparsest. 121 IGH V genes x 26 lengths is 3,146 cells
   against a median 464 IGH clonotypes per real sample — about 0.15 observations per cell, mostly
   structural zeros. That sparsity is real and is not hidden: a corpus of 10,000 repertoires puts
   ~4.6M observations behind the same cells, so the rotation's directions are well determined even
   where one sample's scores along them are noisy. How noisy is readable per sample from
   ``vsig:cov:<locus>:cstar``.

``columns=`` selects output, it does not skip work
--------------------------------------------------

Because a locus is rotated jointly over all its groups, every group must be computed before any of
that locus's components exist. ``columns=`` therefore restricts the emitted frame; it does not avoid
the computation. **Declining a whole locus does skip it**, and is the way to make a run cheaper.

Building a corpus
-----------------

.. code-block:: bash

   vdjtools corpus --corpus synthetic-blood -o synthetic-blood.npz       # 10,000 repertoires
   vdjtools corpus --corpus naive  --out naive.npz
   vdjtools corpus --corpus memory --size n_eff --components 0.95 -o m.npz
   vdjtools corpus --smoke -o /tmp/smoke.npz                             # minutes, for tests

The synthetic corpora use **no samples from anybody's cohort**: every receptor is drawn from the
bundled recombination models.

The build runs across worker **processes**, one contiguous block of samples each; ``--jobs/-j``
sets how many and defaults to every core the process may use (``-j 1`` stays in-process). A sample
is a pure function of its index -- its generator is seeded from ``(seed, locus, index)`` -- so each
worker draws its own repertoires from memory-mapped pools and returns one row of numbers, and
**the artifact is identical at every** ``--jobs``. Measured on Aldan-3, seven loci at
``--samples 10000 --size 10000`` on 72 cores: 157 s to generate the receptor pools (4,850 s in one
process) and about 50 samples/s to featurise them.

.. note::

   Two settings, two layers. ``--jobs`` is worker **processes**; ``POLARS_MAX_THREADS`` and
   ``OMP_NUM_THREADS`` are **kernel threads**. The builder sets its workers to one kernel thread
   each, because ``jobs x cores`` threads is how a 15 s stage was once measured at 296 s.

The build is a deterministic function of ``(corpus, loci, samples, size, seed, source)`` and those
models, and is required to be byte-identical across processes, worker counts and thread counts:

.. code-block:: bash

   vdjtools corpus --corpus naive --loci TRB,TRG --samples 40 --size 150 -o /tmp/a.npz
   OMP_NUM_THREADS=1 POLARS_MAX_THREADS=1 \
     vdjtools corpus --corpus naive --loci TRB,TRG --samples 40 --size 150 -o /tmp/b.npz
   vdjtools corpus --corpus naive --loci TRB,TRG --samples 40 --size 150 -j 4 -o /tmp/c.npz
   cmp /tmp/a.npz /tmp/b.npz && cmp /tmp/a.npz /tmp/c.npz     # all three must be identical

.. important::

   Each repertoire's depth is **drawn**, log-uniformly over the measured real p05–p95 spread for its
   locus, not fixed. A corpus at one fixed depth has exactly zero variance in ``depth:reads``,
   ``depth:richness`` and all five ``pair:`` ratios, so those contribute nothing to the rotation and
   the cross-locus block cannot be fitted at all — and it fails quietly, as an artifact that simply
   omits them. Real repertoires span 4.1x (TRB) to 11.0x (IGH) between their 5th and 95th
   percentiles.

The artifact verifies itself on load: the raw column order is re-derived from the stored vocabulary
through the current layout and compared, and the recorded bundled-model version is checked against
the installed one. Both raise rather than producing numbers in a different coordinate system that
look perfectly reasonable — retraining the bundled models moves the synthetic corpora by 0.1–1.6%
per locus.

Parallelism
-----------

``--jobs`` / ``n_jobs`` selects worker processes: 1 (the default) runs in the
calling process, 0 requests all available cores, and a positive count is limited only
by the number of samples. No timing or sample-depth heuristic changes that count.

Spawned workers inherit one-thread limits for Polars, BLAS, OpenMP and Rayon before
importing numerical libraries. The parent's environment is restored when the pool
exits. Each worker loads its callable and frozen resources once, reads and reduces
one sample at a time, and returns a feature row. Workers take the next available
sample; output rows still follow input order. A failed pool raises.

The CLI passes filenames to workers. For the Python API, use picklable zero-argument
readers for the same bounded-memory path. Already-loaded frames remain resident in
the parent. With ``n_jobs=1``, the caller's kernel thread settings remain in effect.
The default cohort coverage target requires two passes through the input; an explicit
``cstar_target`` avoids the coverage pass without changing worker allocation.

.. _sig-named:

Getting the named statistics back
---------------------------------

A signature is a rotation, and a rotated coordinate has no name a domain reader can use. But the
rotation is fitted on features that *do*: diversity, depth, clone-size fractions, junction length,
isotype composition, SHM, cross-locus yield. Those are computed for every sample either way.

``named=`` returns them:

.. tab-set::

   .. tab-item:: Python

      .. code-block:: python

         from vdjtools.signature import vsig, vsig_cohort
         from vdjtools.signature.corpus import Corpus, bundled_path

         corpus = Corpus.load(bundled_path("blood"))

         vsig(sample, corpus, n_components=32, named=True)
         vsig(sample, corpus, n_components=32, named=("div", "depth"))   # or pick blocks
         vsig_cohort(samples, corpus, n_components=32, named=True)

   .. tab-item:: Command line

      .. code-block:: bash

         vdjtools signature --corpus blood --components 32 --named all -o sig.tsv
         vdjtools signature --corpus blood --components 32 --named div,depth -o sig.tsv

``named=()``, the default, emits exactly the rotated columns and the channels -- the output is
byte-identical to a run without the argument. The reportable blocks are:

.. list-table::
   :header-rows: 1
   :widths: 14 20 66

   * - half
     - block
     - what it holds
   * - ``vsig``
     - ``div``
     - Coverage-standardised Hill numbers ``1D_c``/``0D_c``/``2D_c``, ``0D_chao``, clonality,
       ``d50``
   * - ``vsig``
     - ``depth``
     - ``reads``, ``richness``, ``S_unseen``
   * - ``vsig``
     - ``clon``
     - Singleton, doubleton and top-clone fractions
   * - ``vsig``
     - ``len``
     - Junction-length mean, sd, skew
   * - ``vsig``
     - ``iso`` / ``shm``
     - IGH isotype composition, and mean V identity
   * - ``vsig``
     - ``pair``
     - Cross-locus yield log-ratios
   * - ``rsig``
     - ``depth`` / ``band`` / ``band_igh``
     - ``n_eff`` and observed mass; clone-size and isotype shares of ``Phi``

.. warning::

   **These carry their declared transform, not their natural scale.** A ``log10`` diversity comes
   back as ``log10`` and a ``clr`` composition as a log-ratio with no unique inverse. Do not read
   ``vsig:div:TRB:1D_c = 2.386`` as a clone count -- it is ``log10`` of one.
   :func:`~vdjtools.signature.layout.channel_table` reports the transform for every feature, and
   ``--describe`` prints it per emitted column.

The reason this exists rather than being a convenience: the diversity floor of any study using this
library is made of exactly these blocks. Without a documented route to them there is no floor, and
a signature that beats nothing gets reported as if it beat something.

.. _sig-fit-cohort:

Fitting a corpus on your own cohort
------------------------------------

The published corpora are one route; fitting the rotation on your own training half is the other,
and it is the first thing anybody with a cohort asks for.

.. code-block:: python

   from vdjtools.signature.corpus import Corpus, fit_cohort

   corpus = fit_cohort(train_samples, sig="vsig", name="my-cohort",
                       loci=("TRA", "TRB", "IGH"), n_components=64, n_jobs=0)
   corpus.save("my-cohort.npz")

   # and it is the same artifact type the shipped corpora are
   mine = Corpus.load("my-cohort.npz")
   vsig(held_out_sample, mine, n_components=64)

Samples are ``{sample_id: {locus: frame}}``, or an iterable of frames, or **picklable** zero-argument
callables that defer the read into the worker. ``mir.signature`` must be imported before
``sig="rsig"`` will resolve; nothing in vdjtools imports ``mir``, so the featuriser reaches the
corpus through a registry rather than an import.

.. warning::

   **A cross-validated score computed on a rotation fitted inside the same cohort is not evidence
   that the rotation generalises.**

   Measured on raw repertoire embeddings, 5,376 columns over seven loci, one logistic head, with a
   balanced train/test split of one cohort and a second cohort held out entirely. At matched widths:

   .. list-table::
      :header-rows: 1
      :widths: 40 20 20 20

      * - read-out
        - in-cohort leads at
        - median, in-cohort
        - median, corpus
      * - training out-of-fold (what selects)
        - 4 of 5 widths
        - **0.5039**
        - 0.4968
      * - held-out half
        - 2 of 5
        - 0.5222
        - **0.5635**
      * - external cohort
        - **0 of 5**
        - 0.4757
        - **0.5270**

   Fitting the rotation in-cohort improves the number the configuration is *selected* on, and
   neither held-out read-out. A rotation fitted on 612 samples of one trial learns that trial's
   covariance, which is what a within-cohort cross-validation rewards and what does not travel.

   Both routes are legitimate and they answer different questions. Fit here when you want a basis
   for *this* cohort; name a shipped corpus when the number has to travel.


Further reading
---------------

:doc:`signature-methods` is the measurement record: how each corpus's population was selected, what
the quantile ladders are, how far a corpus carries across sequencing depths, and which alternatives
were tried and rejected. :doc:`channels` documents the columns that are never rotated and never
clamped.

API
---

.. automodule:: vdjtools.signature
   :members:
   :imported-members:

.. automodule:: vdjtools.signature.layout
   :members:

.. automodule:: vdjtools.signature.features
   :members:

.. automodule:: vdjtools.signature.corpus
   :members:

.. automodule:: vdjtools.signature.signature
   :members:

.. automodule:: vdjtools.signature.transform
   :members:
