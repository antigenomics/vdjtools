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

   # four corpora ship with the wheel; name one, no build and no cohort needed
   vdjtools signature --corpus synthetic-blood samples/*.tsv.gz -o sig.tsv

   # what exactly will I get?
   vdjtools signature --corpus synthetic-blood --components 32 --describe

   # or build your own -- still uses no samples from anybody's cohort
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

A corpus is required
--------------------

There is no default. A signature is comparable to another one **only** if both were rotated
through the same corpus, so the artifact has to be named, and the emitted row records which one it
was. Two matrices whose column names match but whose corpora differ are not comparable, and nothing
about the numbers says so — which is why the choice is not allowed to be implicit.

Seven corpora are planned and **four ship today**. Every one of the four is synthetic: each receptor
is drawn from the bundled recombination models, so anyone who installs the library can rebuild the
artifact bit-for-bit, and no sample from anybody's cohort is needed to use one.

.. list-table::
   :header-rows: 1
   :widths: 18 10 20 52

   * - corpus
     - ships
     - what it is
     - richness and clone sizes
   * - ``naive``
     - yes
     - synthetic, unselected
     - every clone size 1 — what the recombination model emits
   * - ``memory``
     - yes
     - synthetic, selected
     - Zipf rank-abundance, sampled by multinomial, zeros dropped
   * - ``synthetic-blood``
     - yes
     - synthetic, blood-like
     - a drawn naive/memory mixture across measured blood ladders (TRB richness 74–3,162,
       ``f1`` 0.215–0.949; n = 34,365 samples)
   * - ``synthetic-tissue``
     - yes
     - synthetic, tissue-like
     - the same across measured tissue ladders (TRB richness 30–2,977, IGH 40–10,352;
       n = 22,298 and 47,031)
   * - ``deep-tcr``
     - not yet
     - real, targeted/amplicon TCR
     - 7 cohorts, 4,080 samples
   * - ``blood``
     - not yet
     - real, bulk RNA-seq blood
     - 23,234 samples across 947 study groups
   * - ``tissue``
     - not yet
     - real, bulk RNA-seq non-blood
     - 13,577 samples across 1,024 study groups

**Use ``synthetic-blood`` or ``synthetic-tissue`` unless you have a reason not to.** ``naive`` and
``memory`` are the two pure regimes, and neither is a bulk sample: ``naive`` has no clone-size
structure at all, and ``memory`` at a fixed nominal size has a read count of exactly
``size × reads_per_clone`` in *every* sample, so its sequencing depth does not vary either. They are
reference points for what selection does to a rotation, not descriptions of a cohort.

.. _sig-cohort-corpus:

The two cohort corpora: three measured ladders per locus
--------------------------------------------------------

A ``synthetic-<cohort>`` corpus draws each repertoire from three quantile ladders measured on the
cohort it is named after — one streaming pass over the harmonized AIRR store, every sample with at
least 100 reads in that locus, five quantiles per locus:

* **richness**, the clonotype count;
* **reads per expanded clone**, ``(reads − singletons) / (richness − singletons)``;
* the **singleton fraction** ``f1``, which in a reasonably deep library is what stands in for the
  naive compartment, and which falls with donor age (`Britanova et al. 2014
  <https://doi.org/10.4049/jimmunol.1302064>`_, `2016
  <https://doi.org/10.4049/jimmunol.1600005>`_).

V and J usage need no ladder: they come out of the rearrangement model.

The mixture is then **constructed** to hit its drawn richness, read count and singleton fraction
exactly — ``round(f1 × N)`` clonotypes at one read, the rest sharing the remainder with Zipf
rank-abundance frequencies on a floor of 2 — rather than sampled and measured afterwards. The three
are drawn through a Gaussian copula at the cohort's measured rank correlations, which is for the
**rotation** rather than for any one marginal: a corpus's rotation is its covariance structure, and
``f1`` against expansion size is −0.503 on blood TRB and −0.846 on blood IGK (n = 43,678), so
independent draws would hand the PCA a correlation the cohort does not have.

Two design notes worth knowing before reading a manifest:

* **Reads per expanded clone, not reads per clonotype.** A repertoire with ``f1`` singletons whose
  other clones all carry at least 2 reads has at least ``2 − f1`` reads per clonotype. Drawing reads
  per *clonotype* independently of ``f1`` therefore lands in a forbidden region about half the time
  on blood TRB; drawing the read count itself does so for 21.1% of 20,000 draws. Reads per
  *expanded* clone is ≥ 2 whatever ``f1`` is, so the trio has no forbidden region.
* **A corpus describes the cohort's middle 90%.** The draw is uniform on ``[p05, p95]`` of each
  ladder, so the cohort's p05 is the corpus's minimum and its median is the corpus's median. A
  sample outside that band is winsorized, and ``vsig:qc:-:winsor_frac`` says how much.

What that costs, measured on 500 drawn samples per cohort against the cohort's own percentiles
(medians; reads is the strict check, being a product of all three drawn quantities):

.. list-table::
   :header-rows: 1
   :widths: 12 8 12 12 12 12 10 10

   * - cohort
     - locus
     - samples
     - richness drawn
     - richness real
     - reads drawn
     - reads real
     - ``f1`` drawn / real
   * - blood
     - TRB
     - 34,365
     - 430
     - 457
     - 780
     - 763
     - 0.734 / 0.759
   * - blood
     - TRA
     - 32,732
     - 260
     - 276
     - 432
     - 456
     - 0.727 / 0.745
   * - blood
     - IGH
     - 33,245
     - 765
     - 697
     - 1,327
     - 1,298
     - 0.668 / 0.697
   * - blood
     - IGK
     - 43,678
     - 554
     - 549
     - 1,280
     - 1,316
     - 0.575 / 0.555
   * - tissue
     - TRB
     - 22,298
     - 150
     - 158
     - 464
     - 308
     - 0.657 / 0.688
   * - tissue
     - IGH
     - 47,031
     - 594
     - 602
     - 2,172
     - 2,194
     - 0.539 / 0.561
   * - tissue
     - IGK
     - 58,706
     - 551
     - 519
     - 2,974
     - 2,634
     - 0.431 / 0.410

The singleton fraction lands within 0.03 everywhere and richness within 10% on 12 of the 14
(cohort, locus) pairs. The read count is within 15% on 11 of 14; the exceptions are the shallow
tissue TR loci, where it runs 23–51% high (tissue TRB 464 against 308) because a product of three
heavy-tailed factors has a median above the product of their medians unless the joint tails match
exactly, and no marginals-plus-copula draw does that. The drawn p95 of the read count is also
*lower* than the cohort's — blood TRB 3,964 against 5,274 — so the fitted depth ceiling is tighter
than the cohort's own p95 by roughly a quarter.

**Why this matters more than it sounds.** ``naive`` and ``memory`` draw their depth across
:data:`~vdjtools.signature.corpus.DEPTH_SPREAD`, 2.4× on TRD to 11.0× on IGH, measured on 1,168 deep
blood samples. Against the full reference cohort the real axis is an order of magnitude wider: blood
TRB richness spans 43× and tissue IGH 259×. Bounds, centre and per-PC scaling are all estimated from
the draw and none of them extrapolates, so a corpus drawn across a tenth of the axis cannot
standardise the rest of it.

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
so the first *n* of a longer rotation are the same vectors. A **count** above the stored one is
capped per locus, exactly as the fit capped it; a **fraction** the stored spectrum cannot reach
raises, and the error quotes what it does reach. Neither ever pads: no column appears under a name
the rotation does not hold.

.. _sig-component-table:

The two forms are not interchangeable — measured
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read from the shipped artifacts' own stored spectra, so this is a measurement rather than a
preference. ``p_raw`` is the locus's raw feature count, ``k`` what the shipped rotation stores, and
``k@0.90`` / ``k@0.95`` the counts those variance fractions would need.

.. list-table::
   :header-rows: 1
   :widths: 8 10 8 10 8 10 10 10

   * - half
     - corpus
     - locus
     - p_raw
     - k
     - k@0.90
     - k@0.95
     - var@k
   * - vsig
     - naive
     - TRB
     - 2261
     - 128
     - 571
     - 772
     - 0.4236
   * - vsig
     - naive
     - IGH
     - 3744
     - 128
     - 676
     - 839
     - 0.4116
   * - vsig
     - naive
     - TRG
     - 767
     - 128
     - 282
     - 365
     - 0.7841
   * - vsig
     - memory
     - TRB
     - 2261
     - 128
     - 501
     - 648
     - 0.4825
   * - vsig
     - memory
     - IGH
     - 3744
     - 128
     - 840
     - 1008
     - 0.3771
   * - vsig
     - memory
     - TRG
     - 767
     - 128
     - 220
     - 340
     - 0.8538
   * - rsig
     - naive
     - TRB
     - 772
     - 128
     - 17
     - 32
     - 0.9923
   * - rsig
     - naive
     - IGH
     - 775
     - 128
     - 12
     - 35
     - 0.9887
   * - rsig
     - naive
     - TRD
     - 772
     - 128
     - 5
     - 11
     - 0.9966
   * - rsig
     - memory
     - TRB
     - 772
     - 128
     - 11
     - 19
     - 0.9957
   * - rsig
     - memory
     - IGH
     - 775
     - 128
     - 8
     - 26
     - 0.9911
   * - rsig
     - memory
     - TRD
     - 772
     - 128
     - 3
     - 4
     - 0.9987
   * - vsig
     - synthetic-blood
     - TRB
     - 2261
     - 128
     - 633
     - 828
     - 0.4629
   * - vsig
     - synthetic-blood
     - IGH
     - 3744
     - 128
     - 836
     - 1078
     - 0.4601
   * - vsig
     - synthetic-blood
     - TRG
     - 767
     - 128
     - 160
     - 280
     - 0.8836
   * - vsig
     - synthetic-tissue
     - TRB
     - 2261
     - 128
     - 489
     - 637
     - 0.5381
   * - vsig
     - synthetic-tissue
     - IGH
     - 3744
     - 128
     - 883
     - 1132
     - 0.5136
   * - vsig
     - synthetic-tissue
     - TRG
     - 767
     - 128
     - 127
     - 246
     - 0.9011
   * - rsig
     - synthetic-blood
     - TRB
     - 772
     - 128
     - 20
     - 36
     - 0.9910
   * - rsig
     - synthetic-blood
     - IGH
     - 775
     - 128
     - 13
     - 38
     - 0.9880
   * - rsig
     - synthetic-blood
     - TRD
     - 772
     - 128
     - 6
     - 17
     - 0.9951
   * - rsig
     - synthetic-tissue
     - TRB
     - 772
     - 128
     - 20
     - 37
     - 0.9909
   * - rsig
     - synthetic-tissue
     - IGH
     - 775
     - 128
     - 13
     - 37
     - 0.9882
   * - rsig
     - synthetic-tissue
     - TRD
     - 772
     - 128
     - 7
     - 17
     - 0.9950

The two halves sit on opposite sides of 128, and by a wide margin:

* **On the statistics half, 128 components are far short of 0.90 nearly everywhere.** Reaching it
  needs 127 (TRG, ``synthetic-tissue``) to 883 (IGH, ``synthetic-tissue``) components, so
  ``--components 0.90`` and anything above it **raises on every shipped vsig corpus and locus except
  that one** rather than returning a narrower matrix. That is the intended refusal — the rotation
  genuinely stops at 128 — and the message names the fraction actually reached. Use a count there, or
  refit at a larger ``n_components``.
* **On the geometry half, 128 is generous.** 0.90 needs 3 to 28 components, because its raw features
  are a few hundred embedding coordinates rather than a few thousand sparse usage shares, so a
  fraction is the natural knob and yields a far narrower matrix.
* **A corpus changes the answer.** IGH needs *more* components under ``memory`` than ``naive`` (840
  against 676) while every other locus needs fewer — clonal expansion concentrates most loci and
  spreads IGH, whose isotype and SHM blocks only vary once clones are selected.
* **Drawing across a real depth range concentrates the TR loci and spreads IGH further.** TRB at 0.90
  needs 489 components under ``synthetic-tissue`` and 633 under ``synthetic-blood``, against 571 under
  ``naive``; TRG drops to 127 and 160 from ``naive``'s 228. IGH goes the other way, 836 and 883
  against 676. Variance at the shipped 128 rises on every locus of both cohort corpora — TRB 0.4629
  and 0.5381 against ``naive``'s 0.4236 — because a corpus drawn across two decades of depth has a
  genuine dominant direction (depth) that a fixed-depth corpus does not.

Full per-locus numbers for all seven loci and both corpora are in each artifact's ``manifest.json``
under ``variance_at_k`` beside the stored spectrum, so any threshold can be read off without a refit.

.. _sig-depth-sweep:

A corpus is depth-portable on the geometry half and not on the statistics half
------------------------------------------------------------------------------

Every bound, centre and scale in a corpus is fitted at the depth the corpus was drawn at, so the
question a user actually has is: *can I score a deep cohort through a corpus drawn shallow?* Measured
rather than assumed — twelve extra corpora at N=1,000, drawn at the per-locus **p05 / median / p95**
of real ``n_eff`` (from 1,168 real bulk blood samples), both halves, both regimes. ``centre shift`` is
the median per-column distance from the p05 centre to the p95 centre, **in p05 robust-SD units**;
``pc scale`` is the median ratio of per-PC scale, p95 over p05.

.. list-table::
   :header-rows: 1
   :widths: 8 8 8 8 8 10 10 12

   * - half
     - corpus
     - locus
     - p05
     - p95
     - bound width
     - pc scale
     - centre shift
   * - vsig
     - naive
     - TRB
     - 195
     - 790
     - 1.48
     - 0.628
     - 3.53
   * - vsig
     - naive
     - IGH
     - 111
     - 1219
     - 2.67
     - 0.311
     - 6.48
   * - vsig
     - naive
     - TRD
     - 14
     - 34
     - 1.10
     - 0.803
     - 1.69
   * - vsig
     - memory
     - TRB
     - 195
     - 790
     - 1.51
     - 0.589
     - 3.57
   * - vsig
     - memory
     - IGH
     - 111
     - 1219
     - 2.64
     - 0.285
     - 5.54
   * - rsig
     - naive
     - TRB
     - 195
     - 790
     - 0.97
     - 1.012
     - 0.06
   * - rsig
     - naive
     - IGH
     - 111
     - 1219
     - 0.96
     - 0.977
     - 0.07
   * - rsig
     - memory
     - TRB
     - 195
     - 790
     - 7.31
     - 1.024
     - 0.05
   * - rsig
     - memory
     - IGH
     - 111
     - 1219
     - 4.62
     - 1.010
     - 0.04

Across all seven loci and both regimes:

* **On the statistics half, depth moves the fit.** The centre shifts **1.67 to 6.48** p05-robust-SD
  between the shallow and deep point, and per-PC scale falls to **0.285–0.904** of its shallow value.
  So scoring a deep cohort through a shallow corpus puts samples many robust deviations out for a
  reason that is sequencing depth, not biology. Match the corpus's depth to the cohort's, or expect
  the offset. The effect tracks each locus's own depth spread — 11.0x on IGH against 2.4x on TRD —
  which is why it is largest on IGH and IGL and smallest on TRD.
* **On the geometry half it does not.** The centre shifts **0.03–0.08** SD and per-PC scale stays
  within **±5%**, at every locus and both regimes. That follows from what the coordinates are: an
  ``rsig`` coordinate is a weighted mean of fixed embedding vectors, so depth changes its *variance*
  (:math:`\approx \sigma^2/n_\text{eff}`) and not its value, while a richness or Hill number is a
  function of depth directly.
* **Bounds are the exception on ``rsig`` ``memory``**, widening **4.6–7.8x** where ``naive`` stays at
  0.94–0.98. Bounds are percentiles, so clonal expansion fattens the tails while leaving the centre
  and scale where they were — the winsorization bound is the depth-sensitive part of the geometry
  half, and the standardisation is not.

The sweep artifacts are not shipped; they are a measurement of the shipped ones. Reproduce with
``vdjtools corpus --corpus naive --size p05 --samples 1000`` (``p05`` / ``n_eff`` / ``p95``).

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

   Two knobs, two layers. ``--jobs`` is worker **processes**; ``POLARS_MAX_THREADS`` and
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
