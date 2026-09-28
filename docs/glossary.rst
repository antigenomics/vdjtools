Glossary
========

Terms as this documentation uses them. Where a word is used inconsistently across the field, the
entry says which convention vdjtools follows, because several of these distinctions change results
rather than only wording.

.. glossary::
   :sorted:

   AIRR
      The Adaptive Immune Receptor Repertoire Community's data
      `schema <https://docs.airr-community.org/>`_. vdjtools standardises on AIRR Rearrangement
      column names, so ``v_call``, ``j_call``, ``junction_aa``, ``duplicate_count`` mean what the
      standard says they mean. See :mod:`vdjtools.io.schema`.

   allele
      A specific sequence variant of a germline :term:`gene`, written with a star suffix:
      ``TRBV9*01``. Recombination models are keyed by allele. Real repertoire files usually carry
      gene-level calls instead, which is a common source of silent zeros -- see
      :term:`gene-level call`.

   anchor
      The conserved residues bounding the :term:`junction`: Cys104 at the V end and Phe118 (or
      Trp118 for IGH) at the J end. Their positions in each germline sequence are what let a V or J
      segment be cut to its CDR3-region contribution. vdjtools resolves anchors from arda.

   ALICE
      Neighbourhood enrichment measured against a V(D)J :term:`generation model` rather than against
      a control cohort: a clonotype with more similar neighbours than its generation probability
      predicts is evidence of antigen-driven expansion. Contrast :term:`TCRnet`. Command:
      ``vdjtools alice``.

   C-star
      Written ``cstar``. The sample :term:`coverage` a repertoire attained, in ``[0, 1]`` -- the
      estimated fraction of the underlying population represented by the observed reads. Used as the
      common level at which :term:`Hill number` s are read, so that two samples sequenced to
      different depths can be compared. Always reported in a :term:`signature`, even when the
      diversity columns it gates are holes.

   CDR3
      The third complementarity-determining region, in the IMGT definition: the hypervariable loop
      that contacts antigen, **excluding** the Cys104 and Phe118 :term:`anchor` residues. Two
      residues shorter than the :term:`junction` for the same receptor.

   clonotype
      One distinct receptor sequence in a :term:`repertoire`, and one row of the canonical frame.
      What counts as distinct is a choice: vdjtools keys on the nucleotide junction plus V and J
      calls by default, and amino-acid keying is available where the analysis calls for it.

   clone size
      How many :term:`read` s support a clonotype -- the ``duplicate_count`` column. The
      distribution of clone sizes across a repertoire is heavily skewed, which is why unnormalised
      richness is not comparable between samples.

   cohort
      A set of samples analysed together, usually described by a metadata table with one row per
      sample. Cohort commands either parallelise over sample files or stream a pre-ingested Parquet
      store.

   corpus
      A fitted reference: a large collection of repertoires from which winsorization bounds, a
      per-locus rotation, and per-component centre and scale were all estimated in one pass. A
      :term:`signature` is only comparable to another one rotated through the same corpus. Nine are
      published; see :doc:`signature`.

   coverage
      The fraction of the underlying receptor population that the observed sample represents,
      estimated from the abundance spectrum. The basis for coverage-standardised diversity, and the
      reason :term:`rarefaction` beats fixed-depth downsampling. See :term:`C-star`.

   d50
      A dominance index, in the legacy ``getDxxIndex`` definition vdjtools keeps: rank clonotypes by
      descending count, take the smallest number ``k`` whose cumulative read fraction reaches 50
      percent, and report ``1 - k/Sobs`` -- the fraction of clonotypes **not** needed. Higher means
      more dominated by few clones. A value of 0.98 means 2 percent of clonotypes carry half the
      reads.

   downsampling
      Randomly reducing a sample to a fixed number of reads or clonotypes so that depth-sensitive
      statistics can be compared. Discards data; :term:`rarefaction` estimates the same quantity
      without doing so.

   expanded clone
      A clonotype supported by more than one read. ``(reads - singletons) / (richness - singletons)``
      -- reads per expanded clone -- is at least 2 by construction, which is why it, and not reads
      per clonotype, can be drawn independently of the :term:`singleton` fraction.

   gene
      A germline V, D, J or C segment, named without an allele suffix: ``TRBV9``. See
      :term:`allele`.

   gene-level call
      A ``v_call`` or ``j_call`` naming a :term:`gene` rather than an :term:`allele`. Recombination
      models are keyed by allele, so passing a gene-level name to a Pgen function raises rather than
      silently marginalising over every allele -- an earlier version returned a V/J-agnostic value
      2.38 times too high with no error. Resolve to a representative allele explicitly.

   generation model
      A probabilistic model of V(D)J :term:`recombination`: which segments pair, how many
      nucleotides are deleted from each, and what is inserted between them. Yields :term:`Pgen`.
      Models for all seven human loci ship in the wheel. See :doc:`model`.

   germline
      The unrearranged V, D, J and C segment sequences a receptor is assembled from. vdjtools
      resolves all germline through arda, which is the single source of truth: mixing germline
      sources within one model produces wrong answers that look plausible.

   Hill number
      A one-parameter family of diversity measures indexed by order ``q``, all in units of an
      effective number of clonotypes. ``q = 0`` is richness, ``q = 1`` is the exponential of Shannon
      entropy, ``q = 2`` is inverse Simpson. Higher ``q`` weights abundant clones more.

   holes
      A column that could not be computed, represented as ``nan`` and never as ``0``. "This locus
      was not sequenced", "this sample is too shallow to estimate this", and "this field was absent
      from the input" are different facts, and the ``mask`` channels record which applies. See
      :doc:`channels`.

   incidence
      Whether a clonotype is present in a sample at all, as opposed to how abundant it is.
      Incidence-based association tests ask whether presence co-varies with a phenotype across
      donors. See :mod:`vdjtools.biomarker`.

   junction
      The nucleotide or amino-acid sequence spanning the V-D-J join **including** the Cys104 and
      Phe118 :term:`anchor` residues -- AIRR's ``junction_nt`` and ``junction_aa``. Two residues
      longer than the :term:`CDR3`. vdjtools uses junction columns throughout; confirm which
      convention any external dataset uses before matching against it.

   locus
      Which receptor chain a clonotype belongs to: TRA, TRB, TRG, TRD, IGH, IGK, IGL. Note that TRA
      and TRD share V segments, so a locus must never be re-derived per row from a gene name in a
      frame that is already keyed by locus.

   MAD
      Median absolute deviation, scaled by 1.4826 to match the standard deviation of a Gaussian. A
      robust scale estimate: a :term:`signature`'s per-component scaling uses the median and MAD of
      the corpus's component scores rather than mean and standard deviation, so a few extreme
      samples cannot set the scale.

   metaclonotype
      A group of similar clonotypes treated as one unit -- typically a ball of a given edit distance
      around a CDR3, optionally requiring matching V and J. Recovers signal that exact matching
      misses, because two donors rarely share a receptor exactly but often share a motif.

   N region
      The non-templated nucleotides inserted between V and D, and between D and J, during
      :term:`recombination`. Together with the deletions at each segment end, this is what makes
      the junction hypervariable and what a :term:`generation model` has to model.

   Pgen
      The probability that a V(D)J :term:`recombination` event produces a given receptor sequence,
      marginalised over every recombination scenario consistent with it. Low Pgen means a sequence
      unlikely to arise twice independently, so convergent observations of it are informative.
      Available over nucleotides, amino acids, and the Hamming-1 ball.

   productive
      A property of the **rearrangement**: in frame, with no stop codon, so it can encode a chain.
      Distinct from :term:`functional gene`, and filtering on one says nothing about the other -- a
      productive rearrangement can use a pseudogene V segment. See :doc:`preprocessing`.

   functional gene
      A property of the **germline gene**, in IMGT's classification: F (functional), ORF (open
      reading frame) or P (pseudogene). Orthogonal to :term:`productive`, which is why vdjtools
      keeps the two filters on separate flags.

   public clonotype
      A receptor observed in many unrelated donors. Usually a consequence of high :term:`Pgen`
      rather than of shared antigen exposure, which is why association testing conditions on
      generation probability.

   rarefaction
      Estimating what a diversity statistic would have been at a smaller sample size or a lower
      :term:`coverage`, from the abundance spectrum, without discarding reads. The basis of
      coverage-standardised comparison. Extrapolation is the same machinery run the other way. See
      :func:`vdjtools.stats.inext`.

   read
      One sequencing observation supporting a clonotype, counted in ``duplicate_count``. Where the
      protocol uses molecular barcodes the unit is a UMI rather than a raw read, which matters for
      depth-dependent statistics because UMI counts are not inflated by PCR.

   recombination
      The somatic process assembling a receptor gene: one V, optionally one or two D, and one J
      segment are joined, with nucleotides deleted from each segment end and non-templated
      nucleotides inserted between them. See :term:`N region`, :term:`generation model`.

   repertoire
      The full set of receptor sequences carried by one individual, in one sample, at one time.
      Represented as a table of :term:`clonotype` s with :term:`clone size` s.

   richness
      The number of distinct clonotypes observed. Strongly depth-dependent, so raw richness is not
      comparable across samples; Chao1 and Efron-Thisted estimate what was missed, and
      :term:`rarefaction` puts two samples on a common footing.

   rotation
      The per-locus linear map from raw features to principal components inside a :term:`corpus`
      artifact. A rotation is a function of features and components only: it carries no sample data
      and no cohort labels. Truncating it to fewer components is exact, so a wide artifact serves
      every narrower request.

   signature
      One fixed-width, named, already-standardised feature vector per repertoire. Split across two
      packages: ``vsig`` is the statistics half (vdjtools) and ``rsig`` the geometry half
      (`mirpy <https://github.com/antigenomics/mirpy>`_); they join on ``sample_id``. See
      :doc:`signature`.

   singleton
      A clonotype supported by exactly one read. The singleton fraction drives every coverage and
      richness estimator, and in a synthetic draw it stands in for the naive share of the
      repertoire.

   spectratype
      The distribution of junction lengths across a repertoire, optionally resolved per V gene.
      Skew away from the germline-determined shape indicates selection or clonal expansion.

   tandem D
      A rearrangement incorporating two D segments. Real but rare, and neither OLGA nor IGoR models
      it; vdjtools supports it exactly. Costs roughly 2.5 times a single-D evaluation, and no exact
      skip exists, because essentially no read has zero tandem-D contribution.

   TCRnet
      Neighbourhood enrichment measured against a **control repertoire**: a clonotype with more
      similar neighbours than the control predicts is a candidate for antigen-driven expansion.
      Contrast :term:`ALICE`, which uses a generation model instead. Command: ``vdjtools tcrnet``.

   winsorization
      Clamping a value to a percentile bound estimated on the :term:`corpus`, so that one extreme
      sample cannot dominate a rotation. Which side is clamped follows the feature's declared
      support, not its transform: a non-negative count is clamped at the top only, a
      log-probability at the bottom only, an unbounded score at both. Clamping edits a real
      measurement, so the fraction clamped is always reported in ``qc:-:winsor_frac``.
