Channels: reading a signature
=============================

Twenty names that turn "the model found something" into a sentence with a noun in it.

A signature is a wide vector — 152, 688 or 1,403 columns depending on the tier. Wide vectors score
well and explain badly. **Channels** are the interpretive layer over them: a small, fixed,
purely structural vocabulary, the same for every sample anyone emits, in which a finding can be
stated and compared across labs.

What a channel is
-----------------

A signature column is named ``<sig>:<channel>:<locus>:<feature>``. The second field is the
**channel**: the named group of columns that measures one thing.

.. code-block:: text

   vsig : div : TRB : 1D_c
   ──┬──  ─┬─   ─┬─   ──┬──
     │     │     │      └── feature — the individual number
     │     │     └───────── locus
     │     └─────────────── channel — WHAT IS BEING MEASURED
     └───────────────────── half: statistics (vsig) or geometry (rsig)

Channels are **disjoint and exhaustive**: every emitted column belongs to exactly one, so a set of
per-channel shares sums over the whole vector with nothing left over. That is what makes them
usable as an accounting unit rather than a filing convenience.

The vocabulary is a property of the contract, not of any corpus or any cohort. The same twenty
names describe every sample anyone emits, at every tier, in every locus. Two labs therefore compare
not just numbers but *findings*.

Why the vocabulary exists
-------------------------

Classical repertoire analysis answers "which summary statistic separates my groups?" by running a
fixed menu of named statistics side by side — diversity, clonality, junction length, V usage — one
test per statistic, and reading off which one moved. The menu **is** the explanation: every number
has a name, so a result is a sentence.

A wide feature vector scores better and explains worse. ``signature()`` returns 688 anonymous
columns; a model trained on them reports "it separates the groups", which has no noun in it. The
channel vocabulary restores the noun without giving up the vector:

.. code-block:: text

   without channels   the classifier reaches AUC 0.78
   with channels      the classifier reaches AUC 0.78, carried by IGH diversity and
                      IGH isotype composition; junction length adds nothing

The second sentence is a hypothesis someone can test with a different assay. The first is not.

The twenty channels
-------------------

Printable at any time, from the library or the command line, reading no input at all:

.. code-block:: bash

   vdjtools signature --channels            # one row per channel, and what it measures
   vdjtools signature --channels --tier full   # sized for the tier you will emit

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - Channel
     - Half
     - What it measures
   * - ``mask``
     - ``vsig``
     - Which loci and which statistics this sample can support at all.
   * - ``qc``
     - ``vsig``
     - How far the annotation had to reach: unrecognised V/J calls, non-standard residues.
   * - ``depth``
     - ``vsig``
     - How much was sequenced — reads, observed richness, unseen-species mass.
   * - ``div``
     - ``vsig``
     - Diversity at a fixed coverage level: Hill numbers, clonality, d50.
   * - ``clon``
     - ``vsig``
     - Clonal dominance — the size of the largest clones as a share of the sample.
   * - ``len``
     - ``vsig``
     - Junction length distribution: mean, spread, asymmetry.
   * - ``pair``
     - ``vsig``
     - Relative yield between loci from one library — the α/β, γ/δ and B/T ratios.
   * - ``iso``
     - ``vsig``
     - Isotype composition of the IGH repertoire (class-switch state).
   * - ``shm``
     - ``vsig``
     - Somatic hypermutation load — mean V identity to germline.
   * - ``pgen``
     - ``vsig``
     - How typical the repertoire's rearrangements are under the V(D)J model.
   * - ``kmer``
     - ``vsig``
     - V-gene-and-junction k-mer composition, projected onto a frozen basis.
   * - ``aa``
     - ``vsig``
     - Amino-acid composition of the junction.
   * - ``pchem``
     - ``vsig``
     - Physicochemical profile of the junction (charge, hydrophobicity, bulk).
   * - ``depth``
     - ``rsig``
     - Depth as the geometry sees it: effective clone count and observed mass.
   * - ``div``
     - ``rsig``
     - Sequence-aware dispersion — Rao entropy, which a Hill number cannot express.
   * - ``band``
     - ``rsig``
     - Shares of the repertoire held by clone-size bands and by isotype.
   * - ``contrast``
     - ``rsig``
     - Signed deviation from unselected V(D)J output — the selection imprint.
   * - ``phiv`` / ``phij`` / ``phic``
     - ``rsig``
     - Where the repertoire sits in V-gene, J-gene and junction coordinates.

Two channel names appear on both halves — ``depth`` and ``div`` — and they are **different
measurements of the same idea**, not duplicates. ``vsig:div`` is a Hill number of the clone-size
vector; ``rsig:div`` is Rao entropy in embedding coordinates, which can see that two clonotypes are
one substitution apart. A channel key therefore always carries its half: ``vsig:div``, never ``div``.

Attributability: which channels can name clonotypes
----------------------------------------------------

Every channel declares, at build time, whether it has a **clonotype pre-image** — whether "which
clones drive this" is a well-posed question.

- ``rsig:contrast``, ``rsig:phiv``, ``rsig:phij``, ``rsig:phic`` **are attributable.** They are
  linear functionals of a kernel-mean over clonotypes, so the question has an answer and the
  library will compute it.
- Everything else **is not.** A Hill number is a summary of a distribution and a read fraction is a
  ratio of totals; neither is a sum over clonotypes, so asking which clones drive it is a category
  error rather than an unanswered question.

This is declared, never inferred from the name. Asking an unattributable channel for its drivers
raises instead of returning a plausible-looking list.

Per-locus channels
------------------

A channel spans loci — ``vsig:div`` is the diversity channel of all seven. Ask for the finer
grouping when the locus is part of the finding, which it usually is:

.. code-block:: python

   from vdjtools.signature import channels

   channels("standard")                      # {"vsig:div": [21, 22, ...], ...}
   channels("standard", per_locus=True)      # {"vsig:div:TRB": [...], ...}
   channels(columns=frame.columns[1:])       # index the frame you actually have

"IGH diversity moved and TRB diversity did not" is a different claim from "diversity moved", and
usually the more useful one.

Using a channel map
-------------------

``channels()`` hands back the name → column-index map that a bare matrix does not carry, keyed to
whatever column list you give it:

.. code-block:: python

   import polars as pl
   from vdjtools.signature import channel, channel_table, channels

   F = pl.read_csv("sig.tsv", separator="\t")
   idx = channels(columns=F.columns[1:])           # sample_id dropped
   F[:, [i + 1 for i in idx["vsig:div"]]]          # just the diversity channel

   channel("vsig:div:TRB:1D_c")                    # -> "vsig:div"
   channel_table("standard")                       # the vocabulary, sized for the tier

To ask *which* channel carries a signal, use :func:`mir.signature.channel_spec` and
:func:`mir.explain.channel_report` in mirpy — the ablation machinery lives there, next to the
geometry half. The map above is what it consumes.

Read the two deltas together:

``delta_in``
   Does this channel carry signal **on its own**? Marginal, so it is inflated by correlation
   between channels — two redundant channels both look important.
``delta_out``
   Is this channel's signal **anywhere else**? Conditional, and deflated by the same correlation.

High in / high out means irreplaceable. High in / near-zero out means **redundant** — the signal is
duplicated elsewhere in the vector, which is itself a finding.

Channels and the reference
--------------------------

Channels are orthogonal to which scale reference you use. Choosing ``blood`` over ``deep-tcr``
changes the location and scale a column is standardised against; it does not move a column into a
different channel, add a channel, or remove one. The vocabulary is fixed by the contract, and the
contract is frozen — see :doc:`signature` for what is fitted on data and what is not.
