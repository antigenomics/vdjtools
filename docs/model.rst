Recombination model workshop
============================

A model in :mod:`vdjtools.model` is three things: a **manifest** declaring the recombination Bayes
net, a set of **tidy polars marginal tables** (one per event), and the **germline frames** those
tables are keyed against. Because the probabilities are ordinary DataFrames, everything below is
table-in / table-out — nothing is hidden behind a binary format.

This page covers building a model from your own reference, fitting it to your own sequences,
checking it, comparing two of them, and asking how much diversity one actually describes. For
generation probability and sampling see :doc:`usage`; for the complete surface see :doc:`api`.

.. contents::
   :local:
   :depth: 1


Published model downloads
-------------------------

The `Hugging Face model mirror <https://huggingface.co/datasets/isalgo/airr_feature_signatures>`_
contains the public recombination models, corpus bundles and example sources. Its
dataset card provides selective-download examples. Ordinary bundled-model loading
continues to use the models installed with vdjtools.

The working-model APIs and CLIs default to ``source="isalgo"``.
vdjtools owns model fitting, probability tables, generation and Pgen; ARDA supplies
the germline reference scaffolds, without OLGA probability-parameter seeding.
The replacement working-model catalog targets twelve contexts: all seven human
loci and mouse TRA, TRB, IGH, IGK and IGL. Each manifest describes the parameter
source and reference. Delivery and validation of these twelve parameter sets are
release gates; API changes alone do not establish completed training or publication.
Human bundle keys remain bare loci; mouse keys are prefixed ``mouse_``. Missing
assets fail explicitly and never fall back to a human model.

No mouse TRG or TRD model is bundled. Acquiring suitable 10X or other public
training data for these two loci is a separate future project, not a blocker for
this twelve-context release. Their absence from the bundle does not establish
that suitable public data do not exist.
``source="olga-legacy"`` and CLI ``--source olga-legacy`` retain the frozen human
bootstrap models for oracle checks and historical reproduction, including
``collapse=False`` for allele-level comparisons. ``olga`` is an alias of
``olga-legacy``. Compatibility names ``learned`` and ``arda`` both select fresh
``isalgo`` fits; reproducing their previous parameter families requires a pinned
previous release. Old bundled weights are never substituted for a missing isalgo fit.
Existing frozen corpus manifests retain their recorded source and fingerprints;
the source omitted by older manifests still means OLGA. A custom corpus fitted on
previous ``learned``/``arda`` weights reports a fingerprint mismatch against the
new aliases; use its matching pinned release for reproduction.


Custom germline libraries
-------------------------

:func:`~vdjtools.model.io.from_germline` builds a model scaffold from **any** V(D)J library — your
own FASTA, a population-specific reference, or arda's. ``from_arda`` is a thin wrapper over it.

``from_arda`` defaults to one major allele per gene: a usable functional ``*01``
when present, otherwise the lowest-numbered usable functional allele. Explicit
``gene_level=False`` retains the allele-level reference for diagnostics. Alternative
genes remain distinct candidates. Inference masks and generation-probability calls
resolve allele suffixes to the model's major representative, including a fallback
such as ``*02``. Repeated aliases do not add probability mass; unknown genes raise
instead of silently turning into an unrestricted search. Generated calls use the
same major names as the model's germline tables.

``load_germline`` preserves every locus membership declared by arda's anchor table.
Shared ``TRAV…/DV…`` alleles therefore occur in both TRA and TRD libraries; an
allele's single compatibility ``Anchor.locus`` label does not restrict membership.

.. code-block:: python

   from vdjtools.model import read_germline_fasta, validate_germline
   from vdjtools.model.io import from_germline

   germline = read_germline_fasta("V.fasta", "J.fasta", "D.fasta")   # D optional -> VJ vs VDJ
   validate_germline(germline)          # tidy issue frame; empty means clean
   template = from_germline(germline, locus="TRB", organism="human")

The frame needs ``allele``, ``segment`` (``V``/``D``/``J``) and ``sequence``; ``gene``,
``functional``, ``cdr3_anchor`` and ``full_germline`` are filled in if absent. ``sequence`` is the
**CDR3-region** germline — for V from the conserved Cys104 codon to the 3' end, for J from the 5'
end through the [FW]118 codon — or the full germline for D. If your FASTAs are full-length, pass
``anchors=`` (OLGA's ``*_gene_CDR3_anchors.csv`` format) and they are sliced for you.

.. note::

   The single most damaging mistake in a custom library is a **misplaced CDR3 anchor**: it shifts
   every deletion profile by a constant, and nothing downstream complains. That is why
   :func:`~vdjtools.model.reference.validate_germline` checks that each V starts on a Cys codon and
   each J ends on Phe/Trp, and warns when one does not.

The marginals of a fresh template are placeholders. Their support ranges **bound what EM can later
learn**, which is why ``ins_max`` is a parameter: a model whose insertion support stops at 40 can
never learn a 45-nt N-region.

New germline scaffolds use model schema 2.1 and explicit ``*_first_nt`` tables:
four rows ``nt, p`` for A/C/G/T coded 0/1/2/3. EM learns the first nucleotide of
each nonempty insertion separately from subsequent Markov transitions. A one-base
insertion therefore contributes first-base evidence even though it has no
transition. VD/VJ/DD chains start at the V-proximal side; DJ starts at the
J-proximal side. This independent first-base distribution differs from IGoR's
boundary-conditioned first transition. Unobserved transition conditionals retain
their normalized initial rows so combining insertion lengths cannot lose mass.

Joint single/tandem D histories
--------------------------------

``vdjtools.model.dd.to_history(template)`` explicitly converts a VDJ scaffold to
schema 3.0; inference applies this automatically to new ARDA VDJ scaffolds.
Its ``d_history`` table stores ``j_allele, n_d, d_allele,
d2_allele, p`` and normalizes the complete history distribution within each J.
Single-D rows use an empty string for the absent second D. A J with no physical
tandem path has only single-D mass; D1–D2–J2 remains permitted in TRB.

Human TRD forward coding histories respect the mapped genomic order
TRDD1–TRDD2–TRDD3. All three genes are forward-oriented in
`IMGT's coordinates <https://www.imgt.org/IMGTrepertoire/index.php?group=TRD&repertoire=GenePositions&section=LocusGenes&species=human>`_.
Nonclassical joining of the opposite RSS pair deletes the coding D pair into an
episomal circle rather than retaining reversed D templates in chromosomal VDDJ
(`Cieslak et al., Fig.8 <https://pubmed.ncbi.nlm.nih.gov/25135298/>`_).
The shared support rule applies this curated order only to recognized human
TRDD genes; it does not infer unknown species or custom gene order from names.

Multiple-D incorporation is a normal, abundant feature of TRD; there is no universal
10% upper bound shared with TRB or IGH. See the
`Chudakov-last-author 2023 conference abstract <https://doi.org/10.1055/s-0043-1768549>`_
and the `cross-species primary study <https://pubmed.ncbi.nlm.nih.gov/42358008/>`_.
Validation checks history normalization, genomic support and generation/Pgen
consistency instead of forcing the fitted TRD tandem probability below that ceiling.

Single-D deletion profiles allow a fully deleted D. Separate ``d_tandem_del``
and ``d2_del`` profiles retain at least one nucleotide from each tandem D and
normalize on that support. This partitions latent histories; different histories
can still emit the same nucleotide sequence. Pgen sums all of them, and inference
assigns posterior weights instead of fixing the longest preliminary alignment.
The native scenario search includes both branches and returns ``d2_call``,
``d2_start`` and ``d2_end`` for joint-history models. Its scenario tuple API appends
``d2_call, idx5_2, idx3_2, pos2`` after the existing nine fields; the batch API adds
the corresponding columns. Legacy models retain their existing output layout.

Entropy contributions weight each branch-specific event by its actual activation
probability, summed over J and the complete history table. The global D-count
probabilities are derived summaries, not extra independent factors. Existing
legacy model tables and assets keep their original representation.
Scoring, both samplers and saved models share this distribution. Models without
these events retain the legacy stationary distribution derived from their Markov
matrix; existing OLGA assets and seeded streams are unchanged.


Learning from your own sequences
--------------------------------

:func:`~vdjtools.model.infer.infer_frame` takes a clonotype frame, finds the junction column, and
builds the per-read V/J masks that make EM tractable on a D-bearing locus.
With ``error_rate`` enabled it uses fresh batched ARDA proposals instead of the
frame's preliminary gene assignments. Proposed genes expand to every functional
model allele; missing V/J proposals leave that segment unrestricted. D and D2
remain unrestricted by annotation, including when no second D was called.

.. code-block:: python

   from vdjtools.model.infer import infer_frame, training_frame

   model, report = infer_frame(template, clones, max_iter=15, gene_prior=1.0)
   training_frame(model)       # run, iter, loglik, n_scoreable, rel_change

Train on **non-productive** reads — out-of-frame *or* stop-codon, since both escaped selection.
Keeping only the out-of-frame half conditions the training set on junction length modulo 3, which
the insertion-length model would then happily learn.

``gene_prior`` is a Dirichlet pseudocount over the germline's functional alleles. ``P(V) = 0`` is an
**absorbing state** of this EM — the E-step weights scenarios by ``P(V)``, so a zeroed allele can
never be re-attributed — and without the prior one unlucky iteration deletes a real gene for good.

Fine-tuning is the same call with a warm start:

.. code-block:: python

   tuned, report = infer_frame(model, new_clones, init="template", max_iter=5)

Watching, and surviving, a long fit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``infer_native(..., error_rate=0.01)`` enables a fitted substitution likelihood
on retained germline bases. Each match contributes ``1-error_rate`` and each
mismatch ``error_rate/3``; insertion bases retain their observed Markov emissions.
Alternative trims and insertion lengths contribute posterior weights even when
inserted bases happen to match the germline. This training likelihood is separate
from exact generation Pgen, which remains error-free. The default ``None`` retains
the exact-match inference mode.
In error mode the default alignment initializer uses uniform parameters so an
exact-match vote cannot give a mismatch-only gene absorbing zero choice mass.
An explicit ``init="template"`` keeps the supplied model's support.

The error-rate M-step divides posterior expected mismatches by posterior expected
retained germline bases. ``report.error_rate`` records the updated rate after each
iteration, and saved training runs declare ``error_model``. ``resume`` restores
the last fitted rate unless explicitly overridden. ``report.n_scoreable`` counts
positive-likelihood inputs; subtract it from ``report.n_sequences`` to count
unsupported inputs. ``report.loglik`` is the mean over scoreable inputs, not a
finite full-cohort likelihood when unsupported observations remain.

Native inference uses ``prune_ratio=1e-5`` by default. Set it to ``0`` for exact
inference; generation Pgen is always exact. Max-product bounds identify the best
complete scenario and discard disjoint single-D cells or tandem-D bridges only
when every scenario they contain falls below that ratio. VJ models prune complete
cells. Retained cells can still contain weaker component histories. Likelihood
and sufficient statistics use the same retained set; no fixed insertion-length
or five-base reassignment cutoff is imposed.

This follows the relative-likelihood principle in the
`IGoR implementation <https://github.com/qmarcou/IGoR/blob/ea42ac8958645103d1c77e5b266a9a969b66cf91/igor_src/Genmodel.cpp>`_,
while using a global maximum and factorized groups rather than IGoR's dynamic
traversal. The ratio is **not** a bound on total omitted probability.
``report.pruned_mass_mean``, ``report.pruned_mass_max`` and
``report.n_pruned_observations`` record actual omitted mass and affected input
counts per iteration. These columns also appear in pruned training logs.
``resume`` restores the recorded ratio unless overridden. Bundled builds expose
``EM_PRUNE_RATIO``; ``model learn`` and ``model build`` expose ``--prune-ratio``.

EM on a D-bearing locus runs for a long time — IGH enumerates roughly 1,225 D pairs per read
against TRB's 9 — and the training log only becomes readable once the fit *returns*, so a slow run
and a stuck one look identical. Two things fix that:

.. code-block:: python

   from vdjtools.model.infer import infer_native, print_progress, resume

   model, rep = infer_native(template, seqs, progress=print_progress(), checkpoint="ckpt/IGH")
   # ...interrupted...
   model, rep = resume("ckpt/IGH", seqs, max_iter=10)

``progress`` reports the log-likelihood and its relative change after every iteration — the exact
quantity compared against ``tol``, so you can watch it approach. ``checkpoint`` saves the model
after each iteration (swapped into place, so a kill mid-write leaves the previous checkpoint
loadable), and :func:`~vdjtools.model.infer.resume` picks it up.

Resuming is **exact**: three iterations plus a resumed four reach the same log-likelihood, and the
same tables, as an uninterrupted seven. The checkpoint carries its training log and the resumed run
appends to it, so ``training_frame`` spans every attempt.

From the CLI: ``vdjtools model learn ... -v --checkpoint DIR``, then
``vdjtools model learn ... --resume DIR``.

Each fit appends a run to ``model.training["runs"]``, which
:func:`~vdjtools.model.io.save_model` writes beside the model as ``training.json``. Models that were
never fitted here — every bundled one — simply have ``training is None``.


Checking a model
----------------

:func:`~vdjtools.model.check.check_model` returns a tidy issue frame rather than raising, so every
problem in a model is visible at once and the result sorts, filters and writes like any other table.

.. code-block:: python

   from vdjtools.model import check_model

   issues = check_model(model)
   issues.filter(pl.col("severity") == "error")     # empty == sound

``severity`` is the contract: ``error`` means the model will crash or score wrongly, ``warn`` means
suspicious but usable, ``info`` is a note. The checks exist because each one has produced a silently
wrong answer at some point — a functional gene left at ``P = 0`` makes Pgen exactly zero for every
clonotype using it; deletion mass past a germline's length is unreachable and quietly lost; an
allele present in a marginal but not in the germline crashes the native packer far from the cause.

Unreachable deletion mass is reported as a **fraction per allele**, ranked, rather than one row per
cell: a shared deletion-bin grid across alleles of different lengths always strands a little mass on
the short ones, so what matters is how much.

.. note::

   On a model imported from OLGA this particular finding is **inherited**, and is reported at
   ``warn`` rather than ``error``. OLGA's own marginals carry the same mass — ``IGHV4-30-4*01``
   has ``Pgen`` identically zero in OLGA too — and vdjtools reproduces those arrays bit-faithfully
   on purpose, so correcting it would break the exact-OLGA-Pgen invariant. The gene-collapsed
   models (the default) are clean.

From the command line, ``vdjtools model check`` exits 1 on any error-severity issue, which makes it
usable as a build gate.


Information content and diversity
---------------------------------

Two different questions, reported side by side because they differ by orders of magnitude:

.. code-block:: python

   from vdjtools.model.analyze import total_entropy
   from vdjtools.model.score import diversity

   total_entropy(model)          # per-event contribution, in bits
   diversity(model, n=5000)      # both entropies and both Hill numbers

**Scenario entropy** is the information in one recombination *event*, summed over the Bayes net. It
is an upper bound on the sequence entropy, because different scenarios can produce the same
junction. **Sequence entropy** is the entropy of the junction distribution itself, estimated by
Monte Carlo: sequences drawn from the model *are* distributed as ``Pgen``, so ``E[-log2 Pgen]`` over
generated sequences is an unbiased estimator and the standard error comes free.

From those follow two Hill numbers — ``2^H`` (the usual "~10^x distinct sequences" figure) and
``1 / E[Pgen]`` (how many draws before two coincide, exact because ``Σ Pgen² = E[Pgen]``). Human
TRB comes out at roughly 52 bits per scenario, 45 bits per sequence, and ~3·10¹³ effective
sequences.


Comparing two models
--------------------

.. code-block:: python

   from vdjtools.model.analyze import compare_models, compare_net_dot, compare_usage, render_dot

   compare_models(olga, learned, by="gene")     # per-event tv, tv_max, jsd_bits, support diffs
   compare_usage(olga, learned, "v")            # gene usage side by side
   render_dot(compare_net_dot(olga, learned), "diff.pdf", fmt="pdf")

In 5.0, ``compare_models`` accepts only ``by`` as a keyword argument. Remove the
formerly unused ``labels`` keyword from table comparisons. Graph titles and labels
remain available through ``compare_net_dot(..., labels=("olga", "learned"))``.

Tables are aligned on the **union** of their realization keys with zero fill, so a gene one model
knows and the other does not contributes to the distance instead of vanishing. Conditioned events
are averaged over parent groups **weighted by the parent's marginal**, so a rarely-used V's deletion
profile cannot dominate; ``tv_max`` reports the worst single group, which is what finds the one
broken gene an average hides.

Jensen-Shannon is the headline metric: symmetric, bounded by one bit, and finite when the supports
differ — which is exactly the case here. KL is deliberately not reported, being infinite whenever
one model assigns zero to something the other does not.

Use ``by="gene"`` to compare models built on different germline vintages or sources; an
OLGA-namespace model and an arda-namespace one only line up at gene level.


.. _release-5-4-validation:

Release 5.4 validation
----------------------

All twelve ``isalgo`` models were freshly fitted to their complete nonfunctional
inputs and converged. The seven frozen OLGA bundles retain their original bytes.
Independent raw-OLGA engine comparisons agree on all 3,378 selected queries,
including zeros, at relative tolerance ``1e-9`` and absolute tolerance zero.
The table separates error-free, V/J-marginalized nucleotide Pgen support from
observation likelihood with each fitted model's stored substitution rate.

.. list-table:: Final fitted models and software-validation counts
   :header-rows: 1
   :widths: 18 16 12 16 19 19

   * - Context
     - Fitted clonotypes
     - EM iterations
     - Selected queries
     - Positive exact Pgen
     - Positive observation likelihood
   * - human TRA
     - 36984
     - 16
     - 300
     - 299
     - 300
   * - human TRB
     - 133622
     - 12
     - 300
     - 300
     - 300
   * - human TRG
     - 14780
     - 11
     - 300
     - 299
     - 300
   * - human TRD
     - 11433
     - 17
     - 300
     - 296
     - 300
   * - human IGH
     - 143477
     - 18
     - 300
     - 294
     - 300
   * - human IGK
     - 295371
     - 27
     - 300
     - 296
     - 300
   * - human IGL
     - 39184
     - 35
     - 300
     - 293
     - 300
   * - mouse TRA
     - 52234
     - 9
     - 300
     - 300
     - 300
   * - mouse TRB
     - 30865
     - 9
     - 300
     - 300
     - 300
   * - mouse IGH
     - 156338
     - 27
     - 300
     - 300
     - 300
   * - mouse IGK
     - 2357
     - 22
     - 300
     - 275
     - 300
   * - mouse IGL
     - 128
     - 12
     - 78
     - 6
     - 78

**Scope and interpretation.** Queries are the first 300 eligible unique junctions
from fresh annotations of the first 20,000 reads in each independent
functional-source input; mouse IGL supplies only 78 eligible queries. This source
bucket does not guarantee productivity. All selected queries, including stop or
noncanonical junctions, are retained. Observation scoring performs no M-step and
does not estimate predictive superiority or unbiased nonfunctional held-out
likelihood. Exact Pgen zeros remain visible: retained-germline substitutions can
have positive observation likelihood while error-free Pgen is zero. The fitted
substitution rate does not identify instrument error, polymorphism or SHM.
The ``1e-5`` relative scenario cutoff does not bound total omitted mass;
validation records its actual mean, maximum and affected-query count.

Likelihood and BIC
------------------

.. code-block:: python

   from vdjtools.model.score import compare_pgen, model_fit, pgen_summary

   model_fit(model, held_out_junctions)            # loglik, k, AIC, BIC
   pgen_summary(compare_pgen(a, b, sequences))     # how differently two models score one set

Two conventions matter. **Likelihoods use nucleotide Pgen**, because ``Σ Pgen_nt`` over all nt junctions
is 1, so ``log Pgen_nt`` is a proper log-likelihood. ``Pgen_aa`` sums only the in-frame, stop-free
nucleotide fiber of a translation, so an amino-acid log-likelihood is unnormalized and its missing
constant *differs between models* — it is a relative score on one fixed sequence set, never an
absolute one. And a sequence the model cannot generate gets ``pgen = 0`` with a null
log-probability, never ``-inf``; aggregates run over the scoreable subset and always report
``n_scoreable`` beside ``n``, so a flattering log-likelihood earned on 10% of the data is visible.

The free-parameter count behind BIC is **support-based**: per normalization group,
``max(occupied cells - 1, 0)``, dropping groups that are undefined (all-zero) or unreachable (their
parent has zero marginal). Counting rows instead would put human TRB's ``v_3_del`` at ~3,600
parameters when ~700 are real. It is a lower bound — a support count cannot distinguish a structural
zero from a parameter EM drove to zero — so BIC is comparable only between models counted the same
way, which is the case for any two compared through this function.

.. warning::

   A log-likelihood computed on the sequences a model was **trained** on is that model's own EM
   objective, which EM increases by construction. It validates nothing. Score a held-out set.


.. _olga-caveats:

Known quirks of the OLGA models
-------------------------------

The bundled ``olga`` model set preserves OLGA's published probability tables.
Importing or loading a model also applies documented germline anchor repairs.
For example, human TRB ``TRBV4-3*02`` moves from anchor 267 to 270;
its repaired germline can contribute probability absent from original OLGA output.
Engine oracle tests therefore compare identical raw germlines and probability
tables, using relative tolerances with zero absolute tolerance. Anchor repairs
are tested separately. Equality of probability tables alone does not establish
equality of Pgen after a germline repair.
``appendix/compare_models.py`` keeps this raw engine oracle separate from the
frozen loaded baseline used for functional-source model comparisons. This source
bucket does not guarantee productivity of every selected junction. The validator
reports both gene-call-conditioned and marginalized error-free Pgen, including
zeros, separately from the fitted model's observation likelihood with its stored
substitution rate. Observation scoring uses one unrestricted native E-step with
the declared relative cutoff and records actual omitted mass; it does not refit
parameters or establish an unbiased nonfunctional held-out likelihood.

Insertion-workspace boundary in the oracle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OLGA 1.3.0's codon-group iteration can omit the final nonzero insertion-length
bin for some support lengths. A fixed 27-residue synthetic TRD regression with
insertion support ending at 40 nucleotides exposes this: native Pgen agrees with independent nucleotide scenario
enumeration, while unmodified OLGA gives a smaller probability. Adding three zero
bins to the test-only OLGA workspace restores strict relative agreement without
changing any probability table. The regression checks this explicitly; production
scoring and bundled model parameters are unchanged.
The comparison validator uses this zero-only workspace padding for both VJ and
VDJ oracle calculations. Independent nucleotide enumeration also confirms the
correction on the mouse TRA queries that expose this boundary.

Deletion mass on trims that cannot be reached
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OLGA stores ``P(delV | V)`` on one deletion-bin grid per locus, sized for the longest allele. An
allele whose CDR3-region germline is shorter than that grid can carry probability on trims longer
than it has nucleotides. The Pgen DP never visits those, so the probability is not redistributed —
it is simply absent from every Pgen through that allele.

Measured on the shipped models (fraction of each allele's own deletion distribution that is
unreachable; identical in OLGA's arrays and in ours):

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - allele
     - unreachable
     - consequence
   * - ``IGHV4-30-4*01``
     - 100 %
     - ``Pgen`` is **identically zero** — in OLGA too
   * - ``IGKJ4*02``
     - 80.9 %
     - Pgen through this J is roughly 5x too low
   * - ``TRAV20*03``
     - 54.7 %
     - Pgen through this V is roughly 2x too low
   * - ``IGHV3-30-3*01``
     - 52.0 %
     -
   * - ``TRAV36/DV7*03``
     - 14.0 %
     -
   * - ``IGLJ3*01``
     - 10.6 %
     -

Across the whole set this affects a minority of alleles per locus (IGH V is the worst: 60 of 62
alleles carry some, most of it small). :func:`~vdjtools.model.check.check_model` reports it as
``deletion_unreachable`` at ``warn`` severity for an OLGA-sourced model, with the fraction, so you
can see whether a gene you care about is affected.

**What to do about it.** Use the gene-collapsed models — the default — where the effect is absorbed
by the collapse. The replacement ``isalgo`` models are fitted on ARDA germline
scaffolds rather than the OLGA deletion grid. ``learned`` and ``arda`` are aliases
of the new isalgo family; they do not select the old fitted weights. If you need allele resolution *and* an affected gene, be aware Pgen through it is
an underestimate.

Genomically impossible D–J pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the deletional D–J arrangements modeled here, a D can only be joined to a J
lying 3' of it. Almost every locus puts all D genes 5' of all J genes, so
this constrains nothing. **TRB interleaves**:

.. code-block:: text

   5' — TRBV… — TRBD1 — TRBJ1 cluster — TRBD2 — TRBJ2 cluster — TRBC2 — 3'

so a TRBD2–TRBJ1 join is physically impossible. OLGA's model does **not** encode that: it gives
``P(TRBD2*01 | TRBJ1-6*01) = 0.333``, a third of that J's D mass on a pair that cannot occur. The
bundled ``olga`` set reproduces it faithfully and :func:`~vdjtools.model.check.check_model` reports
it as ``impossible_dj_pair`` at ``warn``.

New ``isalgo`` models enforce genomic order on complete histories during
inference, generation, Pgen and normalization. TRBD1–TRBJ2 and
TRBD1–TRBD2–TRBJ2 remain allowed; missing preliminary D2 evidence does not veto
the tandem history.

Earlier single-D fits applied a compatibility mask inside the M-step. Their
historical diagnostic showed why genomic support matters:
a D match is only 10–18 nt — weak evidence on short reads — and an earlier build of the human TRB
model had reached ``P(TRBD2 | TRBJ1-1) = 0.091`` before the constraint existed. Enforcing it also
*improved* the fit, converging in 9 iterations instead of 11 to a slightly better log-likelihood,
because mass that had been leaking onto impossible pairs went to the D that actually generated the
read. These were training-objective observations, not independent validation of
the replacement fits.

:func:`~vdjtools.model.infer.enforce_dj_order` repairs a model fitted before the constraint existed.

Genes with no CDR3-region germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OLGA leaves the CDR3-region germline empty for a handful of ORF alleles while still giving them
usage. They can be drawn by the generator but score ``Pgen = 0``, and human ``TRBV23-1`` — 8.6 % of
a real TRB repertoire — is one of them. ``check_model`` reports these as ``unscoreable_gene_mass``.
:func:`~vdjtools.model.io.from_olga` accepts ``derive_orf=True`` to reconstruct the missing germline
from the full germline and the anchor. Historical OLGA-scaffold learned models used
that option; the replacement isalgo models use ARDA reference scaffolds directly.
The explicit ``olga`` set keeps this extra derivation off. Its documented anchor
repairs still apply; an exact engine oracle additionally needs identical raw
germlines, as described above.

V/J usage is protocol-specific
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not a defect, but the most common way to get a wrong answer with these models. OLGA's were fit to
DNA-multiplex data; the historical ``learned`` family used 5'RACE reads. The two amplify different V genes at very
different rates. An earlier comparison measured TRBV19 usage at 3.1 % in OLGA
and 37 % in its 5'RACE population; those are historical population-specific
observations, not validation results for the replacement fits. Neither marginal
is automatically right for *your* library. The junction model (trims, insertions, dinucleotides) is the
shared, transferable part. Use :func:`~vdjtools.model.rescale.rescale_usage` before scoring.

Out-of-frame input
~~~~~~~~~~~~~~~~~~~

``olga``'s own ``compute_nt_CDR3_pgen`` rejects an out-of-frame junction outright ("Invalid
nucleotide CDR3 sequence"). vdjtools scores it, because the generation model is defined before
selection and out-of-frame rearrangements are exactly what it is trained on. This is a deliberate
difference in input validation, not in the probability: on any sequence OLGA accepts, the two agree
exactly.


Extending the allele library
----------------------------

.. code-block:: python

   from vdjtools.model.infer import extend_alleles

   bigger = extend_alleles(model, load_germline("TRB", "human"))
   bigger, report = infer_frame(bigger, clones, init="template", max_iter=5)

A new allele of a gene the model already has is seeded from a gene-mate; a brand-new gene is seeded
from the germline-nearest existing allele at a floor mass, because there is no evidence at all for
how often it is used. Deletion rows copied from a donor are clipped to the new allele's own germline
length, so an extension can never introduce unreachable mass.

Each pre-existing **gene** keeps its total usage: alleles of one gene are alternative versions of
the same gene — a diploid carries at most two — so a richer library must split a gene's mass more
finely, never multiply it. Existing alleles are never modified, including their germline.

This seeds; it does not estimate. Follow it with a warm-start fit.


Rescaling V/J usage
-------------------

.. code-block:: python

   from vdjtools.model import rescale_usage

   scored_with = rescale_usage(load_bundled("TRB", "isalgo"), my_sample)

V and J usage is **protocol-dependent** and the junction model is not: 5'RACE and DNA-multiplex
amplify different V genes at very different rates, so neither usage marginal is right in general.
Learn the junction model once, then set ``P(V)`` from the repertoire you are actually about to
score. Pass an out-of-frame sample only if that is what you will score — a pseudogene's
rearrangements are never productive and are therefore enriched out of frame.


Exporting and importing tables
------------------------------

.. code-block:: python

   from vdjtools.model.io import load_model, marginals_frame, save_model, set_marginals

   marginals_frame(model)                        # every marginal as one long frame
   set_marginals(model, edited_frame)            # ...and back again
   save_model(model, "mymodel", fmt="tsv")       # a hand-editable model directory
   load_model("mymodel")                         # format auto-detected, dtypes restored


Building from the read corpus
-----------------------------

Default raw-read builds select one usable reference allele per V/D/J gene
**before EM**: ``*01`` when usable, otherwise the lowest-numbered usable allele.
Empty, ambiguous and misanchored germlines are excluded.
Allele-only duplicate clonotypes merge with exact read-count sums;
comma-separated gene ambiguity and distinct duplicated genes are retained.
An explicit template keeps its declared resolution and event graph.

Prediction accepts a bare gene or any allele suffix when that gene has one model
representative: ``TRBV9``, ``TRBV9*01`` and ``TRBV9*99`` condition on the same entry.
Unknown genes raise in native Pgen rather than silently marginalizing over all genes.
Explicit multi-allele models retain their allele-specific probabilities.
The training inputs are the original nonfunctional read sets (out of frame or
carrying a stop codon); a junction-only stop/frame filter does not reproduce
whole-rearrangement nonfunctionality.

The EM-fitted model sets shipped in the wheel use real 5'RACE reads; ``olga`` is an import.
The fitting pipeline —
fetch FASTQ, map with arda, collapse to unique clonotypes, run EM — is
:func:`~vdjtools.model.data.build_all`, parallel across chains:

.. code-block:: bash

   vdjtools model build --chains TRB,TRA,IGH --workers 4 -o models/

It needs HuggingFace access to the source dataset and arda's mmseqs2.

For examples and tests, two arda-mapped clonotype sets (human TRA and TRB, out-of-frame) ship in
the source tree as gzipped FASTA and load with no network, no arda and no mmseqs2:

.. code-block:: python

   from vdjtools.model.data import load_prepared

   clones = load_prepared("human", "TRB")     # junction, v_call, j_call, d_call, d2_call, count
   model, report = infer_frame("TRB", clones, max_iter=10)

The V/J/D calls ride in the FASTA header (``>{i}|{v_call}|{j_call}|{d_call}|{d2_call}|{count}``)
because EM needs them for its per-read masks. These files live under ``tests/`` and are not
packaged into the wheel; from an installed vdjtools use
:func:`~vdjtools.model.data.prepare` or pass ``path=``.

.. note::

   Sequencing nucleotide inputs use A/C/G/T only. Standard DNA IUPAC ambiguity
   symbols are replaced with ``A`` before alignment, clonotype grouping, training
   and nucleotide scoring. This is a mandatory input rule: ambiguity is not a
   latent observation alphabet, and no completion enumeration or correction model
   is used. The compatibility ``ambiguous`` argument accepts only ``"A"``.
   Other invalid DNA symbols are input errors. Amino-acid inputs and immutable
   germline references retain their own alphabets. Original raw files and quality
   strings are preserved; normalization applies to the observations consumed by
   the library.


Command line
------------

Every operation above has a ``vdjtools model`` counterpart. A model is named either as a directory
or as ``LOCUS[:source[:organism]]``:

.. code-block:: bash

   vdjtools model check TRB:isalgo                       # exits 1 on any error
   vdjtools model template --locus TRB -o template/
   vdjtools model learn clones.tsv -t template/ -o fitted/
   vdjtools model log fitted/                              # loglik per iteration
   vdjtools model entropy TRB:olga-legacy --table total
   vdjtools model diversity TRB:olga-legacy -n 5000
   vdjtools model compare TRB:olga-legacy TRB:isalgo --by gene --dot diff.pdf
   vdjtools model compare-pgen TRB:olga-legacy TRB:isalgo seqs.tsv --summary
   vdjtools model loglik seqs.tsv TRB:isalgo              # loglik, k, AIC, BIC
   vdjtools model extend fitted/ --locus TRB -o extended/
   vdjtools model rescale TRB:isalgo my_sample.tsv -o rescaled/
   vdjtools model export TRB:olga-legacy --format tsv -o trb_tsv/
   vdjtools model net TRB:olga-legacy --format pdf -o bn.pdf

Options that change a model
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every subcommand takes ``--organism`` (``human`` by default). The rest are grouped by what they
affect; ``vdjtools model <sub> --help`` remains authoritative for defaults.

.. list-table::
   :header-rows: 1
   :widths: 22 24 54

   * - Subcommand
     - Option
     - What it changes
   * - ``template``, ``extend``
     - ``--germline-v``, ``--germline-j``, ``--germline-d``, ``--anchors``
     - Use your own FASTA segment libraries and anchor table instead of arda's. Supply all of them
       together: mixing germline sources within one model produces wrong answers that look
       plausible.
   * - ``template``
     - ``--ins-max``
     - The largest non-templated insertion the scaffold allows, per junction. Raising it enlarges
       every insertion table and slows every later evaluation.
   * - ``learn``
     - ``--template``
     - The scaffold to start from, if not the bundled one for that locus.
   * - ``learn``
     - ``--max-iter``, ``--tol``
     - When expectation-maximisation stops: an iteration cap, and the log-likelihood change below
       which it is converged.
   * - ``learn``
     - ``--init``
     - How the marginals are seeded before the first E-step.
   * - ``learn``
     - ``--gene-prior``, ``--nd-prior``
     - Pseudocounts on gene usage and on the number of D segments. Non-zero keeps a rarely-observed
       gene from being driven to exactly zero, which the E-step treats as absorbing.
   * - ``learn``
     - ``--single-d``
     - Fit single-D only, skipping tandem D. Roughly 2.5 times cheaper and biases ``P(n_D = 2)``
       downward by a measurable amount --- see :term:`tandem D`.
   * - ``learn``
     - ``--no-calls``
     - Ignore the input's own V/D/J calls and re-align every read.
   * - ``learn``
     - ``--checkpoint-every``
     - How often a long fit writes its state, so ``--resume`` can pick it up.
   * - ``build``
     - ``--chains``, ``--groups``, ``--cap``, ``--max-iter``
     - Which loci to build, which read groups to draw from, how many sequences per chain, and the
       EM iteration cap.
   * - ``rescale``
     - ``--no-v``, ``--no-j``, ``--aggregate``
     - Which usage tables to replace from your sample, and whether to aggregate alleles to the gene
       before doing it.
   * - ``check``
     - ``--germline``
     - The reference library to audit the model's germline against.
   * - ``loglik``
     - ``--per-sequence``, ``--weights-col``
     - Report one row per sequence rather than a total, and weight each by a count column.
   * - ``compare``
     - ``--dot-format``
     - The image format for ``--dot``.
   * - ``export``
     - ``--long``
     - One long frame of every probability instead of a directory of per-event tables.

Resource budgets and annotation fields
--------------------------------------

``build_model(..., threads=N)`` forwards the budget to annotation, batched ARDA
candidate search and the native E-step. It fits a retained-germline substitution
rate by default (initially 0.01), and retains junctions regardless of preliminary
V/J calls. ``candidate_masks`` returns enumeration proposals plus explicit search
limits, missing-proposal fallback counts and a sequence/mask identity hash. These
are recorded in training reports and checkpoints. Resuming a candidate-restricted
fit requires matching sequences and masks; declaring ``candidate_search={}``
explicitly records a deliberate support change. Candidate omissions can affect
the fitted model and are a computational approximation.

When full AIRR reads are available, their coding-strand sequences provide the gene
proposals. All full-read alternatives for a clonotype are retained and unioned;
they do not replace the observed junction or fix its deletion/insertion boundaries.
Queries run in bounded batches, with a separate digest and full-read/fallback
counts. Junction-only inputs explicitly fall back to junction queries. Short
junctions can lack enough V/J sequence for a proposal; a missing proposal leaves
that segment unrestricted. Clonotype deduplication preserves read support and
uses a deterministic order for tied counts.

New ARDA scaffolds exclude records marked ``status=truncated`` or ``partial``; an
incomplete germline tail is not a complete short allele template.
``infer_native(..., threads=N)`` controls the E-step directly; zero uses the available CPU
allocation. ``build_all(..., workers=N)`` parallelizes builds with one kernel thread each by
default. An explicit ``threads`` override is honored, so account for both layers when choosing it.
Packed-model state is bounded; successive EM iterations do not retain an unbounded model history.

``annotate_junctions`` also returns the stage-one ``v_fix``, ``j_fix``, ``v_canonical`` and
``j_canonical`` fields. Consumers can reuse the repair classifications without rerunning markup.
The repaired amino-acid junction remains authoritative when applying nucleotide flanks.

Native sampling from the command line
-------------------------------------

For large new sequence draws, select the existing native sampler explicitly::

   vdjtools generate --model TRB --number 10000 --seed 7 --engine native --threads 0 -o draws.tsv

``--threads 0`` selects the allocated CPU budget; a positive count selects a smaller
budget. Native draws are invariant to thread count. The default reference sampler
preserves historical single-D seeded draws and runs serially. Both samplers reject
whole tandem proposals lacking a supported second D or a positive contribution from
both Ds. This fixes the tandem seeded stream to match Pgen scenario weights;
explicit OLGA single-D streams remain unchanged. Native and reference engines
use different random streams, so the same seed does not produce identical sequences
across engines. Neither option fits or rebuilds a model.


Explicit selection-mass scaling
-------------------------------

``vdjtools.model.selection_mass(pgen, selection_q=1)`` returns the explicit
product Q times generation mass, without clipping or normalization. It validates
finite nonnegative input and a positive finite multiplier. The shared
``ALICE_TRB_SELECTION_Q`` constant is 9.41, the historical ALICE average on the
OLGA TRB scale; it is not an established TRA factor. ALICE uses this helper for
its Poisson expectation. Empirical TCRnet controls already reflect selection
and must not receive this generative multiplier a second time.


``model.fit_selection_factor(counts, pgen, total_events=N)`` fits an effective
selection multiplier from receptor-event multiplicities using a zero-truncated
Poisson likelihood. Conditioning on positive counts is necessary because the
input contains observed sequences only. Counts must not be sequencing reads.
The unconditioned mean is N times Q times Pgen; fitted Q is model-, species-
and locus-specific. Unsupported zero-mass sequences require an explicit coverage
audit. All-singleton data do not identify a positive estimate. This effective
constant does not establish a sequence-specific thymic survival probability.


ALICE exposure convention
-------------------------

Modern ALICE multiplies joint V/J neighbourhood generation mass by the total
eligible nucleotide-clonotype count in the locus, then by Q. ``n_group`` remains
the descriptive V/J subgroup count; ``n_exposure`` records the actual multiplier.
``legacy=True`` (CLI ``--legacy``) restores the historical subgroup-sized
exposure with joint probability. ``probability_convention`` and ``selection_q``
make the choice explicit in output. Model source remains a separate choice.
Rows without V or J do not contribute to either the degree or exposure population.

Gene-aware joint Pgen uses ``model.native.pgen_aa_gene_batch(model, sequences, v, j)``.
It sums model alleles for gene calls in one native batch, retains input order,
rejects unknown calls, and accepts mismatches=0 (exact) or 1 (closed ball).
ALICE and repertoire annotation share this implementation.
