Repertoire signatures
=====================

A *signature* is a fixed-order, name-addressed feature vector for one repertoire, on a scale a
downstream model can consume without fitting a scaler of its own. ``vdjtools`` emits the
**statistics** half, ``vsig``; the **geometry** half, ``rsig``, comes from `mirpy
<https://github.com/antigenomics/mirpy>`_ and shares this contract exactly.

It is built in three stages, and a **corpus** fixes the last two:

.. code-block:: text

   raw features  ->  winsorization bounds  ->  rotation + per-PC scaling  ->  vsig
   (this sample                (all three from one corpus artifact)
    alone)                                                              +  channels, untouched

Quickstart
----------

.. code-block:: bash

   # build a corpus once -- uses no samples from anybody's cohort
   vdjtools corpus --corpus naive --smoke -o /tmp/naive.npz

   # then score repertoires through it
   vdjtools signature --corpus /tmp/naive.npz samples/*.tsv.gz -o sig.tsv

   # what exactly will I get?
   vdjtools signature --corpus /tmp/naive.npz --components 32 --describe

.. code-block:: python

   from vdjtools.signature import vsig, vsig_cohort
   from vdjtools.signature.corpus import Corpus

   corpus = Corpus.load("/tmp/naive.npz")
   row = vsig({"TRB": trb, "IGH": igh}, corpus)
   frame = vsig_cohort({"S1": sample1, "S2": sample2}, corpus, n_jobs=0)

A corpus is required
--------------------

There is no default. A signature is comparable to another one **only** if both were rotated
through the same corpus, so the artifact has to be named, and the emitted row records which one it
was. Two matrices whose column names match but whose corpora differ are not comparable, and nothing
about the numbers says so — which is why the choice is not allowed to be implicit.

Five corpora are planned; the two synthetic ones need no cohort at all and are reproducible by
anyone who installs the library.

.. list-table::
   :header-rows: 1
   :widths: 14 20 66

   * - corpus
     - what it is
     - clone sizes
   * - ``naive``
     - synthetic, unselected
     - every clone size 1 — what the recombination model emits
   * - ``memory``
     - synthetic, selected
     - Zipf rank-abundance, sampled by multinomial, zeros dropped
   * - ``deep-tcr``
     - real, targeted/amplicon TCR
     - 7 cohorts, 4,080 samples
   * - ``blood``
     - real, bulk RNA-seq blood
     - 23,234 samples across 947 study groups
   * - ``tissue``
     - real, bulk RNA-seq non-blood
     - 13,577 samples across 1,024 study groups

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
  narrower than the request (the five cross-locus ratios against a request of 128) gets its full
  rotation, recorded per locus in the artifact.
* a **fraction** in ``(0, 1)`` — cumulative variance explained. The count then differs per locus
  *and per corpus*, because some corpora are genuinely easier to explain. That is biology, not a
  defect.

The artifact stores the rotation only up to the count it was fitted at, plus the **full eigenvalue
spectrum**. So truncating downward at apply time is free and **exact** — the components are ordered,
so the first *n* of a longer rotation are the same vectors — while asking for more than was fitted
raises, and the error quotes what the spectrum does reach.

.. warning::

   **Do not choose the count by explained variance alone.** A truncated SVD keeps the directions of
   greatest variance *in the fitting corpus*, which for repertoires are depth, V usage and batch. A
   motif carried by a handful of clonotypes in a handful of donors is not one of those.

   Measured on ankylosing spondylitis against healthy, both HLA-B27+, basis fitted on a disjoint
   cohort: the sum of the 17 columns the published motif occupies reads AUC 0.769 at permutation
   *p* = 0.031, while the best of 64 SVD components reads AUC 0.841 at *p* = 0.20. The larger number
   is the meaningless one — a maximum over 64 components reaches \|AUC − 0.5\| = 0.34 under a label
   permutation null.

   For rare, discriminative signal do not rotate at all: take
   :func:`~vdjtools.signature.raw_and_channels` and an L1 model. For a broad compositional shift,
   rotate, and keep components that survive a study-disjoint refit rather than components that reach
   a variance threshold.

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

   vdjtools corpus --corpus naive  --out naive.npz                       # 10,000 repertoires
   vdjtools corpus --corpus memory --size n_eff --components 0.95 -o m.npz
   vdjtools corpus --smoke -o /tmp/smoke.npz                             # minutes, for tests

The synthetic corpora use **no samples from anybody's cohort**: every receptor is drawn from the
bundled recombination models. The build is a deterministic function of
``(corpus, loci, samples, size, seed, source)`` and those models, and is required to be
byte-identical across processes and thread counts:

.. code-block:: bash

   vdjtools corpus --corpus naive --loci TRB,TRG --samples 40 --size 150 -o /tmp/a.npz
   OMP_NUM_THREADS=1 POLARS_MAX_THREADS=1 \
     vdjtools corpus --corpus naive --loci TRB,TRG --samples 40 --size 150 -o /tmp/b.npz
   cmp /tmp/a.npz /tmp/b.npz          # must be identical

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

``--jobs``/``n_jobs`` is worker **processes**. There is exactly one concurrency knob, and its help
text says which layer it reaches; a flag named ``--threads`` that started processes is a mistake
this subsystem has made before. Workers are spawned, not forked — polars cannot be combined with
``fork`` — and a pool that cannot start **raises** rather than falling back to one process, because
a correctness-preserving fallback is what hid a 20x slowdown here once.

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
