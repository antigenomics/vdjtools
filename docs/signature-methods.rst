How signatures were fitted
==========================

This page is the measurement record behind :doc:`signature`: what was measured, on how much data,
and which alternatives were tried and rejected. You do not need it to use a signature -- the
:doc:`guide <signature>` is self-contained -- but you do need it to defend a number in a manuscript,
or to decide whether a corpus applies to a cohort it was not fitted on.

Everything here is measured. Where a design decision looks arbitrary, the section says what was
measured to settle it, and what the rejected alternative cost.

.. contents::
   :local:
   :depth: 1

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
     - 452
     - 457
     - 812
     - 763
     - 0.774 / 0.759
   * - blood
     - TRA
     - 32,732
     - 253
     - 276
     - 421
     - 456
     - 0.766 / 0.745
   * - blood
     - IGH
     - 33,245
     - 715
     - 697
     - 1,187
     - 1,298
     - 0.687 / 0.697
   * - blood
     - IGK
     - 43,678
     - 578
     - 549
     - 1,260
     - 1,316
     - 0.569 / 0.555
   * - tissue
     - TRB
     - 22,298
     - 156
     - 158
     - 467
     - 308
     - 0.704 / 0.688
   * - tissue
     - IGH
     - 47,031
     - 626
     - 602
     - 1,934
     - 2,194
     - 0.541 / 0.561
   * - tissue
     - IGK
     - 58,706
     - 566
     - 519
     - 2,596
     - 2,634
     - 0.417 / 0.410

The singleton fraction lands within 0.03, and richness within 10%, on **all 14** (cohort, locus)
pairs. The read count is within 15% on 10 of 14; the exceptions are the four shallow tissue TR loci,
where it runs 20–52% high (tissue TRB 467 against 308) because a product of three heavy-tailed factors
has a median above the product of their medians unless the joint tails match exactly, and no
marginals-plus-copula draw does that. The drawn p95 of the read count is also *lower* than the
cohort's — blood TRB 3,755 against 5,274 — so the fitted depth ceiling is roughly 29% tighter than the
cohort's own p95.

**Why this matters more than it sounds.** ``naive`` and ``memory`` draw their depth across
:data:`~vdjtools.signature.corpus.DEPTH_SPREAD`, 2.4× on TRD to 11.0× on IGH, measured on 1,168 deep
blood samples. Against the full reference cohort the real axis is an order of magnitude wider: blood
TRB richness spans 43× and tissue IGH 259×. Bounds, centre and per-PC scaling are all estimated from
the draw and none of them extrapolates, so a corpus drawn across a tenth of the axis cannot
standardise the rest of it.

.. _sig-real-corpus:

The three real populations: how each was selected
-------------------------------------------------

Five of the nine published corpora are real, drawn from three populations --- ``blood`` and
``tissue`` each ship in a capped and an uncapped variant of the same population, and ``deep-tcr``
in one. A real corpus is fitted on repertoires, so what it means depends entirely on which
repertoires. The selection is a rule, and each artifact records it:

* **blood** -- public bulk RNA-seq, blood compartment. **tissue** -- the same, non-blood compartments.
  Both **task-disjoint** by study group *and* by sample against the evaluation panel, because a
  reference fitted on studies it will later be scored against makes the transfer claim circular.
* **deep-tcr** -- targeted/amplicon deep TCR sequencing, TRA and TRB only.
* A per-locus **read floor**, applied to reads *measured off the raw data*.

**Why measured reads and not the metadata flag.** The previous reference slice selected on
``has_TRB | has_IGH``, and that field carries two incompatible definitions: ``>= 100 reads`` in the
base metadata table, and "any row present" in roughly 24,600 appended rows. A population built on it
admits thousands of ultra-shallow repertoires. The harmonized sample table has no per-locus read
column at all, so the floor is applied against a read count measured in one pass over the raw
clonotypes. Blood comes out at 22,441 task-disjoint samples where the old slice recorded 23,234 — that
difference is the defect being removed, not a change of scope.

**TRG and TRD are a stratum with a lower floor, and the number is measured.** They cannot support 100
reads at scale:

.. list-table::
   :header-rows: 1
   :widths: 14 14 18 18 18 18

   * - cohort
     - locus
     - samples ≥100
     - study_ids ≥100
     - samples ≥30
     - study_ids ≥30
   * - blood
     - TRG
     - 5,686
     - 442
     - 13,381
     - 705
   * - blood
     - TRD
     - 3,563
     - 331
     - 9,525
     - 580
   * - tissue
     - TRG
     - 367
     - 78
     - 1,574
     - 317
   * - tissue
     - TRD
     - 378
     - 42
     - 1,059
     - 189

So the five main loci keep the 100-read floor and TRG/TRD use 30. Because each locus is fitted only on
the rows that observed it, this is a stratum rather than a filter on the whole corpus, and
``n_obs`` per column records exactly how many samples each column rests on. At a 30-read floor a TRG
repertoire has a few dozen clonotypes: its richness, V usage and junction-length features are
computable while the coverage-standardised diversity columns stay holes, which is the honest outcome
rather than an extrapolated one.

**Both cap variants ship.** Repertoire corpora are dominated by a handful of large submissions — the
top 10 hold 17.7% of SRA blood, the top 50 hold 40.9% — so the predicate caps each study at 30
samples, applied *before* any other row limit. The uncapped fit ships beside it so the cap's effect is
measurable rather than asserted; emission is the expensive half and is done once, so each extra fit
costs seconds.

**What the manifest does not contain.** No cohort label, no dataset name, no accession, no sample id.
A rotation is a function of features and components; a cohort label is not one of its inputs. What is
recorded is the population as rule text, the per-locus floors, aggregate sample and study counts, the
measured read band per locus, and the germline fingerprint the V/J and spectratype columns are indexed
by — which is the load gate even for a real corpus, since nothing was drawn from those models but the
columns are still named by them.

The two forms are not interchangeable — measured
------------------------------------------------

Read from the shipped artifacts' own stored spectra, so this is a measurement rather than a
preference. Every corpus is fitted at **k = 256** and ships every component; ``k@0.90`` is what that
variance fraction would need on the same spectrum.

.. list-table::
   :header-rows: 1
   :widths: 22 12 12 12 12 12 12

   * - corpus
     - vsig TRB k@0.90
     - vsig TRB var@256
     - vsig IGH k@0.90
     - vsig IGH var@256
     - rsig TRB k@0.90
     - rsig TRB var@256
   * - ``naive``
     - 571
     - 0.6283
     - 676
     - 0.5768
     - 14
     - 0.9996
   * - ``memory``
     - 501
     - 0.7305
     - 840
     - 0.5411
     - 14
     - 0.9996
   * - ``synthetic-blood``
     - 633
     - 0.6392
     - 844
     - 0.5926
     - 14
     - 0.9996
   * - ``synthetic-tissue``
     - 483
     - 0.7320
     - 885
     - 0.6263
     - 18
     - 0.9995
   * - ``blood``
     - 464
     - 0.7525
     - 951
     - 0.5596
     - 14
     - 0.9996
   * - ``blood-uncapped``
     - 485
     - 0.7386
     - 973
     - 0.5464
     - 14
     - 0.9995
   * - ``tissue``
     - 301
     - 0.8541
     - 1,069
     - 0.5627
     - 18
     - 0.9995
   * - ``tissue-uncapped``
     - 311
     - 0.8438
     - 1,148
     - 0.5487
     - 18
     - 0.9994
   * - ``deep-tcr``
     - 611
     - 0.7501
     - —
     - —
     - 7
     - 1.0000

``p_raw`` is 2,261 on TRB and 3,744 on IGH for every corpus (the layout's width does not depend on
what it was fitted on); ``rsig`` carries 772 and 775.

The two halves sit on opposite sides of 256, and by a wide margin:

* **On the statistics half, even 256 components are short of 0.90 everywhere.** Reaching it needs 301
  (TRB, ``tissue``) to 1,148 (IGH, ``tissue-uncapped``), so ``--components 0.90`` **raises** rather
  than returning a narrower matrix. That is the intended refusal — the rotation genuinely stops at 256
  — and the message names the fraction actually reached. Use a count, or refit wider.
* **On the geometry half, 256 is far past the point of diminishing returns.** 0.90 needs **7 to 18**
  components and 0.95 needs 10 to 32, because its raw features are a few hundred embedding coordinates
  rather than a few thousand sparse usage shares. ``--components 32`` there gives a far narrower matrix
  at essentially no cost in explained variance, and is usually what you want.
* **A real corpus concentrates TRB and spreads IGH.** ``tissue`` reaches 0.90 on TRB in 301 components
  against ``naive``'s 571, while its IGH needs 1,069 against 676. Real IGH carries isotype and SHM
  structure, and real libraries carry between-study heterogeneity, that a generative draw does not
  produce; real TRB in tissue is dominated by a depth axis that concentrates it.
* **The cap barely moves the rotation.** Capping blood at 30 samples per study changes TRB's
  ``k@0.90`` from 485 to 464 and its variance at 256 from 0.7386 to 0.7525 — the cap is a guard
  against a handful of submissions dominating the *centre and scale*, not a lever on the components.

**Widening a fit is a strict superset — measured, not assumed.** Refitting ``synthetic-blood`` from
k=128 to k=256 reproduces all 128 original components with ``|cos| = 1.000000`` on TRB, IGH and TRG and
identical eigenvalues to ``0.000e+00``. That is why shipping 256 costs existing users nothing: their
numbers at 128 are bit-for-bit the same, and truncating downward at apply time needs no refit.

Full per-locus numbers for all nine corpora are in each artifact's ``manifest.json`` under
``variance_at_k`` beside the stored spectrum, so any threshold can be read off without a refit.

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


Three alternatives were tried and rejected
------------------------------------------

Each was measured before being dropped, and each is recorded so that it is not reinstated as an
apparent simplification.

**Reads per clonotype as the third drawn quantity.** Because
``reads >= singletons + 2*(richness - singletons)`` is an identity, mean reads per clonotype is at
least ``2 - f1``, where ``f1`` is the singleton fraction. An independently drawn pair therefore
violates the constraint about half the time on blood TRB, and drawing the read count itself violates
it for 21.1 percent of 20,000 draws. Clipping back into the feasible region silently rewrites the
marginal it was drawn from. Reads per **expanded** clone is at least 2 whatever ``f1`` is, so it has
no forbidden region, and its Spearman correlation with log richness is -0.036 on blood TRB against
-0.103 for reads per clonotype.

**A three-point quantile ladder.** The derived read count is a product of all three drawn quantities
and so is the sharpest check on the draw. Against blood TRB's measured median of 763 reads, a
five-point ladder gives 762 and a three-point ladder gives 846; blood IGH's measured 1,298 comes out
1,276 against 1,490. Resolution, not parameterisation, is what made the corpora land.

**Independent draws instead of a copula.** The copula was kept, and for the rotation rather than for
any marginal: a corpus's rotation *is* its covariance structure, and the singleton fraction against
expansion size correlates -0.503 on blood TRB and -0.846 on blood IGK. On the read-count median alone
independence makes no measurable difference (735 against the copula's 762, with 763 measured), so the justification is the
covariance, not the margin.

Two implementation notes that follow from reproducibility rather than from statistics:

The copula's correlation root is a closed-form 3x3 Cholesky factor, never ``numpy.linalg.eigh``.
Eigenvector **signs** are a LAPACK convention rather than mathematics: negating a column leaves the
covariance, and therefore every marginal and every rank correlation, exactly as it was, while
changing the realised sample. An artifact would then be statistically identical and byte-different
between two machines whose LAPACK disagreed, which is the one thing a shipped corpus may not be.
Measured: flipping one column moved the realised rank correlation on blood TRB from -0.052 to -0.028
against a target of -0.036 -- nothing detectable in the distribution, and a completely different
draw.

The mixture is **constructed rather than sampled**: ``round(f1*N)`` clonotypes at one read, and the
rest sharing the remaining reads with Zipf rank frequencies on a floor of 2. So a drawn sample hits
its richness, read count and singleton fraction exactly. ``memory`` instead drops a multinomial's
zeros, which is why its richness and singleton fraction come out of the sampling rather than being
asked for -- and therefore why a ``memory`` corpus cannot be pointed at a measured cohort.
