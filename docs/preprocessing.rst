Data pre-processing
===================

Everything that happens between a file on disk and a clonotype frame you can compute on: reading
and format conversion, the three filtering axes, error correction, depth normalisation, and
frequency handling.

All of it is free functions over the canonical clonotype frame (:mod:`vdjtools.io.schema`), and
all of it is mirrored on the command line.

.. contents::
   :local:
   :depth: 2


Reading and format conversion
-----------------------------

:func:`vdjtools.io.read` sniffs the format and returns a normalised frame. The canonical columns
are ``v_call, d_call, j_call, c_call, junction_aa, junction_nt, duplicate_count, frequency``.

.. code-block:: python

   from vdjtools import io

   df = io.read("sample.tsv")                  # format sniffed
   df = io.read("sample.tsv", fmt="mixcr")     # or stated

Readers ship for AIRR, vdjtools' own format, parquet, MiXCR, MiGEC, MiTCR, immunoSEQ/Adaptive,
IMGT, Vidjil, RTCR, TRUST4 and arda. Adaptive material is remapped to IMGT names on the way in,
because Adaptive's gene nomenclature is not IMGT's and comparing the two without remapping silently
compares different genes.

Keeping a non-canonical column
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``keep=`` carries extra columns through the normalisation:

.. code-block:: python

   df = io.read("sample.tsv", keep=("v_identity",))

``v_identity`` is the usual case — the signature's somatic-hypermutation block cannot be computed
without it, and it is not part of the canonical schema.

Optional AIRR annotation columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three AIRR Rearrangement fields are *optional* and not part of ``SCHEMA``, but are used when
present: ``productive``, ``stop_codon`` and ``vj_in_frame``. See
:ref:`filtering-productive` — where a file states its own productivity, that statement wins over
anything re-derived from ``junction_aa``.

Whole cohorts
~~~~~~~~~~~~~

:func:`vdjtools.io.ingest_cohort` and :func:`vdjtools.io.scan_cohort` read a metadata sheet plus a
directory of samples; :func:`vdjtools.io.iter_samples` and :func:`vdjtools.io.map_samples` stream
them one at a time when the cohort does not fit in memory.


.. _filtering-three-axes:

Filtering — three axes, and they are not the same axis
------------------------------------------------------

This is the part that most often goes wrong, because one English word has been doing three jobs.
Keep them apart:

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - axis
     - the question
     - the standard's word
   * - the **file**
     - is ``junction_aa`` parseable at all?
     - — (not a biological question)
   * - the **rearrangement**
     - does this sequence encode a chain?
     - AIRR ``productive``
   * - the **germline gene**
     - is the V/D/J/C gene it uses real?
     - IMGT functionality: **F** / **ORF** / **P**

`AIRR Rearrangement <https://docs.airr-community.org/en/latest/datarep/rearrangements.html>`_
defines ``productive`` as an open reading frame with no defect in the start codon, splicing sites
or regulatory elements, no internal stop codon, and an in-frame junction.
`IMGT functionality <https://www.imgt.org/IMGTindex/functionality.php>`_ classifies a *germline
gene* as **F** (functional), **ORF** (has a reading frame but a defect in splicing, regulatory
elements or conserved-residue hydropathy — *not* functional) or **P** (pseudogene).

They are orthogonal. A rearrangement can be perfectly in frame with no stop codon while using a
pseudogene V; a functional V gene can rearrange out of frame. **Filtering one tells you nothing
about the other.**

.. _filtering-productive:

Productive — a property of the rearrangement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from vdjtools.preprocess import filter_productive, productive_mask

   df = filter_productive(df)                              # keep productive
   df = filter_productive(df, keep="nonproductive")        # the complement
   df = filter_productive(df, recompute_frequencies=False) # leave frequencies alone

**The file's own annotation wins.** If ``productive`` is present it is used; failing that,
``stop_codon`` and ``vj_in_frame`` together; and only if none are present is productivity derived
from ``junction_aa``, where a stop codon is ``*`` and an out-of-frame junction carries one of the
legacy markers ``[atgc#~_?]``. :func:`productive_mask` returns the predicate *and* the evidence it
rests on, so you can always ask which was used:

.. code-block:: python

   mask, source = productive_mask(df)
   print(source)        # 'productive' | 'stop_codon+vj_in_frame' | 'junction_aa'

The derived fallback is a proxy, and worth knowing the limits of: it cannot see a defect in a
splicing site or a regulatory element, which AIRR's ``productive`` can.

.. note::

   ``filter_functional(keep="coding")`` is a deprecated alias and warns. "Functional" is IMGT's
   word for a gene; this filters rearrangements. Map ``keep="coding"`` to ``keep="productive"``.

Functional genes — a property of the germline gene
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from vdjtools.preprocess import filter_functional_genes

   df = filter_functional_genes(df)                          # keep IMGT F only
   df = filter_functional_genes(df, keep=("F", "ORF"))       # the common looser choice
   df = filter_functional_genes(df, segments=("V", "J", "D"))

A call that cannot be resolved against the germline reference is **kept**, not dropped. An
unrecognised gene name means the reference is incomplete or the caller uses a different
nomenclature — discarding those rows would report a vocabulary bug as biology.

Length
~~~~~~

.. code-block:: python

   from vdjtools.preprocess import filter_length

   df = filter_length(df)                              # 5..60 aa, inclusive
   df = filter_length(df, min_len=8, max_len=30)
   df = filter_length(df, keep="outside")              # inspect what a bound would discard

Both bounds are **inclusive**: the defaults keep a 5-mer and keep a 60-mer. They are a data-sanity
bound, not a biological claim — below 5 aa a junction cannot span the Cys104..Phe118 anchors with
any diversity between them, and above 60 aa it is beyond what the germline can produce. Real
junctions sit far inside both.

.. warning::

   These bounds are on ``junction_aa``, which **includes** the Cys104 and Phe118 anchors and is
   therefore two residues longer than the IMGT CDR3. Subtract 2 if you are reasoning in CDR3
   lengths.

   Nothing in this package imposes a length bound by default, so switching this on **will** change
   counts on a corpus that carries junk. Run it with ``keep="outside"`` first.

Frequency, segment and cross-sample filters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from vdjtools.preprocess import filter_frequency, filter_segment, filter_by_sample

   df = filter_frequency(df, min_freq=1e-5)
   df = filter_frequency(df, top_quantile=0.5)
   df = filter_segment(df, v=["TRBV20-1"], j=["TRBJ2-2"])
   df = filter_segment(df, v=["TRBV20-1"], keep=False)      # remove instead
   df = filter_by_sample(df, other)                          # incidence against another sample


Frequencies
-----------

**By default `frequency` is derived from `duplicate_count`, at read time and again after every
filter.** That is the historical behaviour and it is right for the common case, where a file's
frequency is just ``count/total`` anyway.

It is wrong when it is not. A UMI-corrected frequency, or one already normalised against something
other than the row counts, is silently replaced by ``count/total`` and cannot be recovered. To keep
what the file actually said:

.. code-block:: python

   df = io.read("sample.tsv", recompute_frequencies=False)

Honoured by the ``vdjtools``, ``airr`` and ``parquet`` readers. A file carrying no frequency column
derives one either way — there is nothing to preserve. Where ``read_airr`` collapses rows to
clonotypes, a preserved frequency is **summed** across the collapsed rows, since two rows becoming
one clonotype contribute additively to its share.

.. note::

   Before 3.10.0 this was not possible at all: ``_AIRR_ALIASES`` had no ``frequency`` entry and
   every reader passed ``recompute_freq=True`` unconditionally, so "use the frequencies as in the
   file" could not be honoured at any layer above.

Downstream, a frequency is recomputed whenever a filter removes rows, and the filters that take a
switch expose it the same way:

.. code-block:: python

   filter_productive(df, recompute_frequencies=True)    # default: renormalise over survivors
   filter_productive(df, recompute_frequencies=False)   # leave the file's frequencies alone

``recompute_frequencies=True`` is the default because it is what almost every caller wants — after
dropping rows, the remaining frequencies should sum to 1. Pass ``False`` when the frequencies
themselves are the quantity of interest and must stay comparable to the unfiltered file.

.. warning::

   **``recompute_frequencies=False`` is undone by the next filter in a chain.** Only
   :func:`filter_productive`, :func:`filter_length` and :func:`filter_functional_genes` take the
   switch. :func:`filter_frequency`, :func:`filter_segment`, :func:`filter_by_sample` and
   :func:`downsample` renormalise unconditionally, and :func:`select_top` spells the same idea
   ``renormalize=True``. So this preserves nothing::

       df = filter_productive(df, recompute_frequencies=False)
       df = filter_frequency(df, min_freq=1e-5)      # <- renormalises anyway

   If the file's frequencies must survive a chain, do the frequency-preserving filter **last**.

To renormalise explicitly at any point:

.. code-block:: python

   from vdjtools.io import recompute_frequency

   df = recompute_frequency(df)      # frequency = duplicate_count / sum(duplicate_count)

.. note::

   :func:`vdjtools.signature.blocks.work_frame` also writes the ``frequency`` column, but it is
   **not** the file's frequency — it is the signature's internal clone weight
   ``log2(1+count)/Σ``, which merely borrows the column name. It is applied inside the signature
   and does not affect anything documented on this page.


Error correction and decontamination
------------------------------------

.. code-block:: python

   from vdjtools.preprocess import correct, decontaminate

   df = correct(df, max_mismatches=2, ratio=0.05)
   df = decontaminate(df, others=[sample_b, sample_c], ratio=20.0)

:func:`correct` collapses a low-count clonotype into a near neighbour when the neighbour is
sufficiently more abundant — the standard PCR/sequencing-error model. :func:`decontaminate` removes
clonotypes far more abundant in another sample than in this one, the cross-contamination case.


Depth normalisation
-------------------

.. code-block:: python

   from vdjtools.preprocess import downsample, select_top

   df = downsample(df, size=10_000)                 # by reads
   df = downsample(df, size=5_000, by="clonotypes")
   df = select_top(df, n=1_000)

Both recompute ``frequency`` over what survives.

.. warning::

   Order matters. Anything that recomputes ``frequency`` — ``downsample``, ``select_top``, or a
   filter with ``recompute_frequencies=True`` — must run **before** the signature's
   :func:`~vdjtools.signature.blocks.work_frame`, never after, or it silently restores read
   weighting on top of the clone weight.


Pooling, joining and batch correction
-------------------------------------

.. code-block:: python

   from vdjtools.preprocess import pool_samples, join_samples, correct_vj_usage

   pooled = pool_samples([a, b, c], key="aa")            # sum counts
   joined = join_samples([a, b, c], key="aaVJ")          # incidence
   corrected = correct_vj_usage(samples, batch_col="batch")

Match keys are ``strict | nt | ntV | ntVJ | aa | aaV | aaVJ``.
:func:`correct_vj_usage` standardises V/J usage within a batch — for a *named technical* batch
variable, not for a study identifier that is collinear with the biology you are measuring.


.. note::

   **"Non-functional" means the opposite thing one subpackage over.**
   :mod:`vdjtools.model` trains the bundled recombination models on
   ``LABELS = ("functional", "nonfunctional")`` reads, where *nonfunctional* is exactly the
   out-of-frame population :func:`filter_productive` removes — and it wants them, because a
   rearrangement that never met selection is the cleanest read of the recombination process
   itself.

   That is a third sense of the word, on top of IMGT's, and it is why this page uses AIRR's
   *productive* for the rearrangement axis. ``model``'s ``nonfunctional`` is this page's
   ``keep="nonproductive"``.


On the command line
-------------------

.. code-block:: bash

   vdjtools convert sample.tsv -o sample.airr.tsv

   # the two filtering axes are separate flags, deliberately
   vdjtools filter sample.tsv --productive          -o productive.tsv
   vdjtools filter sample.tsv --nonproductive       -o nonproductive.tsv
   vdjtools filter sample.tsv --functional-genes    -o functional_v.tsv
   vdjtools filter sample.tsv --functional-genes --keep-orf -o f_and_orf.tsv

   # length, inclusive bounds
   vdjtools filter sample.tsv --min-len 5 --max-len 60 -o sane.tsv

   # keep the file's own frequencies instead of renormalising
   vdjtools filter sample.tsv --productive --keep-frequencies -o kept.tsv

   # combine
   vdjtools filter sample.tsv --productive --min-len 5 --max-len 60 --min-freq 1e-5 -o clean.tsv

``--coding`` / ``--noncoding`` still work as hidden deprecated aliases and print a notice.
``vdjtools filter --productive`` reports which evidence it used on stderr.


How mirpy differs
-----------------

**In vdjtools, productive filtering is optional.** Non-productive rearrangements are real data —
they carry the second allele's rearrangement, and their share is a measurable per-sample quantity —
so this package will happily compute statistics on them, and ``--nonproductive`` exists precisely
to isolate them.

**In mirpy it is mandatory and cannot be turned off.** mirpy embeds receptors into a geometry, and
a stop codon is *in* the alphabet its distance code uses: an unfiltered frame does not crash, it
returns a finite, meaningless number and contaminates the geometry silently. mirpy therefore
filters on every read and raises if you ask it not to.

If you want the non-productive fraction, that is a vdjtools question.
