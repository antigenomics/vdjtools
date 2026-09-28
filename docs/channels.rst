Channels
========

A **channel** is a named family of columns emitted in its own units, never passing through the
rotation. That is the whole definition, and it is narrower than it was before the 4.0 rewrite: a
channel is not "a block of the signature" but specifically *the part that is not standardised
against a corpus*.

Why channels exist
------------------

A provenance number that has been mixed with the measurements it was supposed to qualify is no
longer provenance. If ``v_fallback_frac`` went through a rotation, the coordinate carrying it would
also carry diversity, usage and length, and a caller could no longer ask "is this row comparable to
ours at all?" -- which is the only question that column exists to answer.

So the split is by *role*, not by cost:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - role
     - treatment
   * - a **measurement** of the repertoire
     - raw feature: transformed, winsorized against the corpus, rotated, scaled
   * - a statement about **the row itself**
     - channel: emitted as-is, in ``[0, 1]`` or a raw count

.. _channels-mask:

``mask`` is a feature block, not diagnostics
--------------------------------------------

.. important::

   **Do not drop the mask columns before modelling.** They are grouped separately from ``qc`` on
   this page for exactly this reason: ``qc`` describes how much to trust the row, and ``mask``
   describes the donor.

   Which loci a sample resolved at all is biology. In a cohort of 874 labelled bulk RNA-seq
   samples, 720 of 7,399 (sample, locus) pairs sit below a five-clonotype floor, almost all of them
   TRD and TRG -- and *which* donors those are tracks lymphocyte content rather than noise.
   Measured on that cohort with a second cohort held out entirely, adding the five presence flags
   to the recommended feature set moved the **external** ROC-AUC from **0.6243 to 0.6676** for 17
   extra columns, and the held-out one from 0.6371 to 0.6419.

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - channel
     - loci
     - what it says
   * - ``vsig:mask:<locus>:present``
     - 7
     - Whether the locus had any productive clonotype at all.
   * - ``vsig:mask:<locus>:estimable``
     - 7
     - Whether the coverage target was reachable without extrapolating. ``present=1`` with
       ``estimable=0`` means the sample is too shallow, which is a different fact from an absent
       locus and used to render as the same ``nan``.
   * - ``vsig:mask:IGH:c_call`` / ``:shm``
     - IGH
     - Whether the input carried a constant-gene call, and a usable ``v_identity``. Both are
       optional AIRR fields, so their absence is a property of the file rather than of the donor.
   * - ``rsig:mask:<locus>:present`` / ``:estimable``
     - 7
     - The same convention on the geometry half. ``estimable=0`` means the locus was present but
       held fewer than ``min_clonotypes`` clonotypes, so its embedding-diversity family is a hole.

``cov`` -- the coverage this sample attained
---------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - channel
     - loci
     - what it says
   * - ``vsig:cov:<locus>:cstar``
     - 7
     - The Chao coverage this sample attained at its observed depth, in ``[0, 1]``. **Always
       emitted**, including when the diversity columns it gates are holes. This is the number that
       makes a diversity hole readable, and the one to check before reading any component of a
       shallow locus.

``div`` and ``disp`` -- embedding diversity, on the geometry half
------------------------------------------------------------------

Rao quadratic entropy is the metric-space analogue of Simpson diversity, and it sees something no
Hill number can: that two clonotypes are one substitution apart. To a counting index they are two
species, as distinct from each other as from anything else.

Until 4.3 Rao was the only such read-out the geometry half gave up; everything else it knew was
inside a rotated component, and a component has no name a domain reader can use. These are the rest
of the analogues, and all of them come out of the same chunked pass that already computes ``Phi``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - channel
     - what it is
   * - ``rsig:div:<L>:rao``
     - Rao's ``Q`` over the whole embedding. Simpson, with distance.
   * - ``rsig:div:<L>:q_v`` / ``q_j`` / ``q_c``
     - The same restricted to the V, J and junction strides of ``Phi``. The strides are literal
       column offsets, so "how much of this diversity is V-driven" is answerable with no
       attribution model. The three sum to ``rao``.
   * - ``rsig:div:<L>:q_frac_v`` / ``q_frac_j`` / ``q_frac_c``
     - Each stride's share of the total dispersion, in clr coordinates -- a composition of *where*
       the diversity sits. No clonotype-counting index has an analogue.
   * - ``rsig:div:<L>:evenness``
     - ``Q(w) / Q(uniform)``: the effect of the clone-size weighting with the composition held
       fixed. Bounded, and much less depth-fragile than richness.
   * - ``rsig:div:<L>:eff_dim``
     - ``exp(H(lambda))`` over the weighted covariance spectrum -- **richness**, as the number of
       independent directions of receptor space the repertoire occupies. A thousand clones inside
       one convergent cluster occupy few; a thousand unrelated clones occupy many; a clonotype
       count cannot tell those apart.
   * - ``rsig:div:<L>:eff_dim_pr``
     - ``(sum lambda)^2 / sum lambda^2``, the order-2 version, dominated by the leading directions.
   * - ``rsig:div:<L>:q_top`` / ``q_singleton``
     - Rao **inside** a clone-size compartment, weights renormalised within it -- the diversity
       *of* the expanded compartment, rather than its share, which ``band`` already carries.
   * - ``rsig:div:<L>:q_ratio_top``
     - ``q_top / q_singleton``, the two compared directly.
   * - ``rsig:disp:<L>:top_singleton`` / ``cos_top_singleton``
     - The distance and the cosine between the expanded and singleton **centroids**. A displacement
       is a quantity no index has: where a compartment sits, rather than how large it is.
   * - ``rsig:disp:IGH:IgG_IgM`` / ``IgA_IgM`` and their cosines
     - The same between isotype compartments -- class-switch geometry.
   * - ``rsig:disp:<L>:norm``
     - ``||Phi||``.

.. warning::

   Every one of these is **within one locus**. Each locus has its own panel of ``K`` prototype
   receptors, so ``Phi(TRA)`` and ``Phi(TRB)`` are vectors in different spaces and a distance
   between them is arithmetic without a meaning. Cross-locus comparison belongs in ``vsig:pair``,
   which is a ratio of scalars.

   A locus below ``min_clonotypes`` holes this entire family and sets ``mask:<L>:estimable = 0``.
   It is not zeroed: Rao of a one-clonotype locus is arithmetically ``0.0`` and is not a diversity
   measurement, it is the presence mask in different units -- and read as a measurement it sat
   about 29 robust deviations below the 1st percentile of the real values. The geometry itself
   (``phiv``/``phij``/``phic``) is still computed there, because ``Phi`` is perfectly measurable
   from three clonotypes even when its dispersion is not.

``qc`` -- diagnostics
---------------------

These describe the **row**, not the donor. Read them before trusting a row; do not model on them.

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - channel
     - loci
     - what it says
   * - ``vsig:qc:<locus>:v_fallback_frac`` / ``:j_fallback_frac``
     - 7
     - Weight fraction whose V or J call is outside the arda germline vocabulary. An unrecognised
       call raises nowhere -- the germline lookup falls back to the maximum observed distance -- so
       a cohort on Adaptive nomenclature or an older IMGT release yields a fully populated, entirely
       plausible, systematically wrong result. This is the only column that says so.
   * - ``vsig:qc:<locus>:nonstd_aa_frac``
     - 7
     - Weight fraction dropped as non-productive or unparseable. ``nan`` under
       ``prefiltered=True``, because a caller who already filtered would otherwise get a confident
       floor of zero.
   * - ``vsig:qc:-:n_loci_present``
     - --
     - How many of the seven loci this sample had.
   * - ``vsig:qc:-:winsor_frac`` / ``rsig:qc:-:winsor_frac``
     - --
     - What fraction of this row's finite values the corpus's bounds clamped. See
       :ref:`sig-winsor`. A value near 1.0 means this sample does not belong to this corpus.

Reading a row
-------------

.. code-block:: python

   from vdjtools.signature import vsig
   from vdjtools.signature.corpus import Corpus, bundled_path

   corpus = Corpus.load(bundled_path("blood"))   # a published name; fetched on first use
   row = vsig(sample, corpus)

   # diagnostics -- read these before trusting the row, do not model on them
   row["vsig:qc:-:winsor_frac"]          # did the corpus's bounds edit this sample?
   row["vsig:qc:TRB:v_fallback_frac"]    # is our germline vocabulary this cohort's vocabulary?

   # coverage -- how completely was this locus sampled?
   row["vsig:cov:TRB:cstar"]

   # mask -- a FEATURE. Which loci this donor resolved is biology, not bookkeeping.
   row["vsig:mask:TRD:present"]
   row["vsig:mask:TRB:estimable"]        # is vsig:pc:TRB:* resting on a real diversity estimate?

What can I get, and in what units
----------------------------------

One call, rather than three lookups and a reading of the source:

.. code-block:: python

   import polars as pl
   from vdjtools.signature.layout import channel_table

   pl.DataFrame(channel_table("vsig"))

Each row is one ``(block, feature)`` with its ``transform``, its ``support``, the loci it is
emitted for, and a ``kind``:

.. list-table::
   :header-rows: 1
   :widths: 16 84

   * - kind
     - meaning
   * - ``rotated``
     - A rotation input. You do not get this column; you get ``<sig>:pc:<locus>:PCnn``.
   * - ``named``
     - A rotation input that is **also** reportable on its own, via ``named=`` -- see
       :ref:`sig-named`. Same number, carrying its declared transform.
   * - ``channel``
     - Carried through untouched: never winsorized, never rotated.

``vdjtools signature --describe`` prints the same three kinds for the exact columns an invocation
will emit.

Channels are declared, not detected
-----------------------------------

:func:`~vdjtools.signature.layout.channels` is the registry and
:func:`~vdjtools.signature.channel_columns` enumerates it. A channel declares a ``support`` like a
raw feature does, but purely for documentation and a range assertion -- no trimming is applied
either way.

``mir.signature`` registers its own channels into the same registry with
:func:`~vdjtools.signature.register_channel`. Nothing in vdjtools imports ``mir``.

.. code-block:: python

   from vdjtools.signature import channel_columns, channels, support_of

   # One entry per DECLARATION, not per name: `mask` and `qc` are each declared twice, once for
   # all seven loci and once for the IGH-only or cross-locus features.
   [c.name for c in channels("vsig")]     # ['cov', 'mask', 'mask', 'qc', 'qc']
   len(channel_columns("vsig"))           # 46
   support_of("vsig:cov:TRB:cstar")       # 'unit'

.. note::

   **A channel is pass-through, so adding one does not invalidate a fitted corpus.** The rotation
   is indexed by :func:`~vdjtools.signature.layout.raw_columns` and nothing else, and
   :func:`~vdjtools.signature.corpus.apply` fills every registered channel from the sample. The
   ``div``, ``disp`` and ``mask`` families added to ``rsig`` in 4.3 therefore appear in the output
   of every artifact already published, with no refit and no new download.
