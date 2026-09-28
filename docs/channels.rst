Channels
========

A **channel** is a named family of columns that is emitted in its own units and never passes through
the rotation. That is the whole definition, and it is a narrower thing than it was before the 4.0
rewrite: a channel is not "a block of the signature" but specifically *the part that is not
standardised against a corpus*.

Why channels exist
------------------

A provenance number that has been mixed with the measurements it was supposed to qualify is no
longer provenance. If ``v_fallback_frac`` went through a rotation, the coordinate carrying it would
also carry diversity, usage and length, and a caller could no longer ask "is this row comparable to
ours at all?" — which is the only question that column exists to answer.

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

The seven vsig channels
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 30 10 60

   * - channel
     - loci
     - what it says
   * - ``vsig:cov:<locus>:cstar``
     - 7
     - The Chao coverage this sample attained at its observed depth, in ``[0, 1]``. Always emitted.
       This is the number that makes a diversity hole readable, and the one to check before reading
       any component of a shallow locus.
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
   * - ``vsig:qc:<locus>:v_fallback_frac`` / ``:j_fallback_frac``
     - 7
     - Weight fraction whose V or J call is outside the arda germline vocabulary. An unrecognised
       call raises nowhere — the germline lookup falls back to the maximum observed distance — so a
       cohort on Adaptive nomenclature or an older IMGT release yields a fully populated, entirely
       plausible, systematically wrong result. This is the only column that says so.
   * - ``vsig:qc:<locus>:nonstd_aa_frac``
     - 7
     - Weight fraction dropped as non-productive or unparseable. ``nan`` under ``prefiltered=True``,
       because a caller who already filtered would otherwise get a confident floor of zero.
   * - ``vsig:qc:-:n_loci_present``
     - —
     - How many of the seven loci this sample had.
   * - ``vsig:qc:-:winsor_frac``
     - —
     - What fraction of this row's finite values the corpus's bounds clamped. See
       :ref:`sig-winsor`. A value near 1.0 means this sample does not belong to this corpus.

The geometry half adds ``rsig:div:<locus>:rao`` — Rao quadratic entropy in embedding coordinates,
carried through rather than rotated — and its own ``rsig:qc:-:winsor_frac``.

Reading a row
-------------

.. code-block:: python

   from vdjtools.signature import vsig
   from vdjtools.signature.corpus import Corpus, bundled_path

   corpus = Corpus.load(bundled_path("blood"))   # a published name; fetched on first use
   row = vsig(sample, corpus)

   # before trusting any vsig:pc:* value of this row
   row["vsig:qc:-:winsor_frac"]          # did the corpus's bounds edit this sample?
   row["vsig:cov:TRB:cstar"]             # how completely was TRB sampled?
   row["vsig:mask:TRB:estimable"]        # is vsig:pc:TRB:* resting on a real diversity estimate?
   row["vsig:qc:TRB:v_fallback_frac"]    # is our germline vocabulary this cohort's vocabulary?

Channels are declared, not detected
-----------------------------------

:data:`~vdjtools.signature.layout.channels` is the registry, and
:func:`~vdjtools.signature.channel_columns` enumerates it. A channel declares a ``support`` like a
raw feature does, but purely for documentation and a range assertion — no trimming is applied either
way.

``mir.signature`` registers its own channels into the same registry with
:func:`~vdjtools.signature.register_channel`. Nothing in vdjtools imports ``mir``.

.. code-block:: python

   from vdjtools.signature import channel_columns, channels, support_of

   [c.name for c in channels("vsig")]     # cov, mask, qc
   len(channel_columns("vsig"))           # 46
   support_of("vsig:cov:TRB:cstar")       # 'unit'
