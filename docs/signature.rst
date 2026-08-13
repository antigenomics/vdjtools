Repertoire signature — the statistics half
==========================================

:mod:`vdjtools.signature` turns one AIRR sample into a **fixed, named, positional** vector of
repertoire statistics — ``vsig``. Its companion is :mod:`mir.signature` (``rsig``), which covers the
embedding geometry; the two are namespaced so they concatenate on ``sample_id`` without colliding,
and the shared contract machinery (column layout, transform registry, frozen-reference rescaling)
lives here, in vdjtools, because mirpy depends on vdjtools and not the reverse.

.. code-block:: python

   from vdjtools.signature import vsig, vsig_cohort
   from vdjtools.signature import layout

   v = vsig({"TRB": df})                       # {column: value}, in layout order
   F = vsig_cohort(samples, tier="standard")   # one row per sample, positional
   layout.columns("full")                      # the contract, computed from no data

Columns are named ``vsig:<block>:<locus>:<feature>``. Tiers ``core`` / ``standard`` / ``full`` are
**exact index subsets** of one frozen layout — a narrower tier is a slice of a wider one, never a
differently-computed number.

.. contents::
   :local:
   :depth: 1

Every column is transformed, and the transform is denominator-aware
-------------------------------------------------------------------

A learner cannot be handed a read count, an amino-acid fraction and a Hill number in one matrix and
be expected to weight them sensibly. Each block therefore declares a transform:

.. list-table::
   :header-rows: 1
   :widths: 16 22 62

   * - transform
     - blocks
     - why this one
   * - ``arcsine``
     - ``aa`` (×20)
     - Anscombe :math:`\arcsin\sqrt{(xm + 3/8)/(m + 3/4)}` — a proportion with a known denominator
   * - ``logit``
     - ``qc`` ``shm`` ``div:clonality`` ``pgen:frac_atypical`` ``clon:top``
     - Haldane–Anscombe, so ``0/3`` and ``0/500`` are different numbers
   * - ``clr``
     - ``iso`` ``clon`` (ships :math:`k-1` parts)
     - a genuine closed composition; all :math:`k` parts are linearly dependent
   * - ``log10`` / ``log1p``
     - ``depth`` ``div``
     - non-negative, right-skewed, unbounded above
   * - ``none``
     - ``len`` ``pchem`` ``pair`` ``pgen:mean_log10``
     - already a moment, a log-ratio or a log-probability
   * - ``none``, exempt
     - ``mask``
     - a 0/1 indicator; rescaling it against a reference would destroy it

Why arcsine and not log1p with winsorization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This comes up every time, so here are the measurements (``benchmark_transform_arcsine.py`` in
`2026-mirpy-analysis <https://github.com/antigenomics>`_ regenerates them).

**Depth invariance.** Hold the true proportion :math:`p` fixed and vary the denominator from
:math:`m = 91` (the corpus median junction-residue count for a shallow RNA-seq sample) to
:math:`m = 50{,}000`:

.. list-table::
   :header-rows: 1
   :widths: 26 24 24 26

   * - transform
     - value at m = 91
     - value at m = 50,000
     - drift
   * - ``log1p(count)``
     - 0.919
     - 6.908
     - **7.5×**
   * - ``arcsine(p, m)``
     - 0.1476
     - 0.1419
     - **1.04×**

The same biological composition sequenced twice as deep is a *different feature value* under
``log1p`` of counts. That alone decides it: depth is the largest nuisance axis in a public-data
cohort, and a transform that encodes depth into every column hands the learner the batch label.

**Variance stabilisation.** Ratio of the largest to the smallest per-bin standard deviation across
the 20 amino-acid columns:

.. list-table::
   :header-rows: 1
   :widths: 20 26 26 28

   * - denominator
     - ``log1p``
     - ``arcsine``
     - reading
   * - m = 91
     - 3.1–3.9×
     - 3.1–3.9×
     - indistinguishable when counts are small
   * - m = 5,000
     - 21.6×
     - **1.0×**
     - only arcsine survives depth

**Zeros carry information.** ``arcsine(0, m)`` depends on :math:`m` — a zero from 91 residues and a
zero from 50,000 residues are different evidence, and the transform says so. ``log1p(0) = 0``
always, which asserts that never observing a residue in a shallow sample and never observing it in
a deep one are the same observation.

**Winsorization is a fitted parameter.** A 5% clip has to be estimated from a corpus, so it inherits
that corpus's depth distribution — and at the top of the distribution it clips the dominant V genes
and the most abundant residues, which is where the biology usually is. Arcsine has no free
parameter: the :math:`3/8` is Anscombe's, not ours.

**Can you compose them?** ``arcsine(log1p(x))`` is undefined for :math:`x \ge 2` — ``log1p(2) = 1.10``
and :math:`\arcsin` needs :math:`[0,1]`. The useful composition is the one the pipeline already
does: log the *counts* into a clone weight, normalise to a proportion, then arcsine that proportion
against the true count denominator. Under heavy clonal expansion the ``log2p1`` clone weight cuts
the column standard deviation about 4× (0.0319 → 0.0080) and restores arcsine's 1.14× stabilisation,
which raw frequency weighting loses.

Holes are ``nan``, never 0
--------------------------

A locus that was not sequenced, or a statistic the sample is too shallow to estimate, yields ``nan``
plus a ``vsig:mask:`` column. A model that reads "absent" as "zero" reads an unsequenced chain as a
biological finding.

Ambiguous V calls resolve to the first gene
-------------------------------------------

Real data does not hand you one V gene. Adaptive/immunoSEQ material realigned with MiXCR against the
IMGT reference — the better reference, but junction plus short flanks rather than full primered
reads — routinely calls a comma-separated tie. In one store slice that produced **1,296 distinct
"genes"**, of which **1,235 were comma-strings** rather than genes at all, which silently shatters
any V-keyed feature space into singleton columns.

:func:`vdjtools.io.schema.resolve_gene` takes the **first** gene in the list and strips the allele —
MiXCR orders by alignment score, so the first is the best-supported call:

.. code-block:: python

   from vdjtools.io.schema import resolve_gene, strip_allele

   df.with_columns(resolve_gene(pl.col("v_call")))   # "TRBV5-1*01,TRBV5-5*01" -> "TRBV5-1"
   df.with_columns(strip_allele(pl.col("v_call")))   # keeps the tie, sorted, for reporting

The two are deliberately different functions. ``strip_allele`` **sorts** the tie so that reporting is
order-insensitive; composing it with "take the first" would therefore give you the *alphabetically*
first gene (``TRBV19`` where the aligner said ``TRBV5``), which is a plausible-looking wrong answer.
Use ``resolve_gene`` for anything keyed on V.

The V + k-mer space
-------------------

:mod:`vdjtools.features.kmer_space` is the heavy feature: a junction k-mer profile keyed jointly on
the V gene, TF-IDF scaled, and projected onto a frozen truncated-SVD basis. The k-mer counting is
C++ (``src/kmer.cpp``); a Python implementation is not viable at corpus scale.

.. code-block:: python

   from vdjtools.features.kmer_space import fit_kmer_space, save_kmer_spaces
   from vdjtools.signature.kmer import register_kmer

   sp = fit_kmer_space(frames, pattern="xxxx", n_groups=8, flank=4,
                       min_df=0.02, max_df=0.80, n_components=32)
   register_kmer({"TRB": sp})            # adds vsig:kmer:TRB:PC01.. to the layout
   sp.transform(df, weight="freq")       # project one repertoire

**Patterns.** ``"xxxx"`` is a contiguous 4-mer; ``"xx.x"`` is gapped (the ``.`` is a skipped
position). Gapped patterns exist because a conserved motif with one variable position is otherwise
split across every substitution at that position.

**Alphabets.** ``n_groups=20`` is the plain amino-acid alphabet. Smaller values cluster residues by
BLOSUM62: :func:`seqtree.SubstitutionMatrix.blosum62().penalty` is
:math:`s(a,a) + s(b,b) - 2s(a,b)`, which is a *squared* distance in the Gram sense, so classical MDS
recovers the exact embedding and Ward linkage on it gives the groups. No hand-picked chemistry
classes.

**The sparsity trade, measured.** On TRB with V-keying and a 4-mer:

.. list-table::
   :header-rows: 1
   :widths: 18 26 28 28

   * - alphabet
     - code space
     - surviving columns
     - effect
   * - A = 20
     - 3,360,000
     - 1,022
     - most codes never occur; the survivors are common motifs
   * - A = 8
     - 86,016
     - 25,402
     - dense, but distinct motifs collapse into one column

Neither is universally right. ``max_df`` matters as much as the alphabet: at ``A=8, max_df=0.80``
a known disease motif was cut from the vocabulary entirely (0/4 windows present), while at
``A=20, max_df=0.99`` all 4 survived. If you are looking for a *rare discriminative* motif, keep the
document-frequency ceiling high and do not project — see below.

Do not select components by variance
-------------------------------------

A truncated SVD keeps the directions of greatest variance **in the fitting corpus** — which, for a
repertoire matrix, are sequencing depth, V-gene usage and batch. A disease motif carried by a
handful of clonotypes in a handful of donors is not one of those, so projecting onto the leading
components discards it by construction.

Measured on ankylosing spondylitis vs healthy, both HLA-B27+, with the space fitted on a disjoint
TRB cohort (``benchmark_kmer_ankspond.py``):

.. list-table::
   :header-rows: 1
   :widths: 44 16 16 24

   * - read-out
     - AUC
     - perm. p
     - verdict
   * - sum of the 17 columns the published motif occupies (A = 20)
     - 0.769
     - **0.031**
     - the representation carries it
   * - same, A = 8
     - 0.731
     - 0.061
     - the reduced alphabet blurs it
   * - best of 64 SVD components
     - 0.841
     - 0.20
     - **chance** — null median 0.297, 95th 0.385
   * - best single vocabulary column
     - 0.813
     - —
     - a label-selected maximum, not an estimate

The two "best-of-N" rows are the point. They are the larger numbers and they are the ones that mean
nothing: a maximum taken over 64 components reaches :math:`|\mathrm{AUC} - 0.5| = 0.34` under a
label permutation null. **Any best-of-N read-out must be nulled or not quoted.**

How many components, then? Measured on the emitted signature matrix (14,553 samples × 1,369
columns, 182 studies), mean AUC over the four largest tasks with study-disjoint folds and the
rotation refit inside every fold:

.. list-table::
   :header-rows: 1
   :widths: 16 14 14 14 14 14 14

   * - components
     - 8
     - 16
     - 64
     - 256
     - 512
     - all 1,369
   * - mean AUC
     - 0.575
     - 0.581
     - 0.589
     - **0.591**
     - 0.584
     - 0.555

Flat from 16 to 256, then a cliff — **the full matrix scores worse than 16 components**. Use
**16–64**: 16 buys 98% of the achievable AUC at a quarter of the width, 64 is the plateau. Keep all
columns only for a regularised learner (L1, boosting) or a rare sparse signal.

The practical consequence for feature selection:

- **rare, discriminative signal** — keep the un-projected TF-IDF columns and use an L1 model. The
  SVD is the wrong tool; it optimises for the wrong thing.
- **broad compositional shift** — project, and keep only components that survive a study-disjoint
  split-half refit (per-component score correlation ≥ 0.95), not components that reach a
  cumulative-variance threshold.
- **either way** — report the permutation null alongside the number.

Regenerating these numbers
--------------------------

Every table on this page is produced by a script in the companion
`2026-mirpy-analysis <https://github.com/antigenomics>`_ repo, so it can be re-measured rather than
trusted:

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - script
     - what it re-measures
   * - ``benchmark_transform_arcsine.py``
     - the arcsine vs log1p depth-drift and variance-stabilisation tables
   * - ``benchmark_kmer_ankspond.py``
     - the k-mer transfer AUCs and their permutation nulls
   * - ``benchmark_signature_dimension.py``
     - how many components reproduce across a study-disjoint refit
   * - ``fit_kmer_spaces.py``
     - the fitted per-locus spaces and their vocabulary sizes

*Measurements on this page were last taken 2026-08-13.*

API
---

.. automodule:: vdjtools.signature
   :members:
   :undoc-members:

.. automodule:: vdjtools.features.kmer_space
   :members:
   :undoc-members:
