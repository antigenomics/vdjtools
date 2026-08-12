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


Custom germline libraries
-------------------------

:func:`~vdjtools.model.io.from_germline` builds a model scaffold from **any** V(D)J library — your
own FASTA, a population-specific reference, or arda's. ``from_arda`` is a thin wrapper over it.

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


Learning from your own sequences
--------------------------------

:func:`~vdjtools.model.infer.infer_frame` takes a clonotype frame, finds the junction column, and
builds the per-read V/J masks that make EM tractable on a D-bearing locus.

.. code-block:: python

   from vdjtools.model.infer import infer_frame, training_frame

   model, report = infer_frame(template, clones, max_iter=15, gene_prior=1.0)
   training_frame(model)       # run, iter, loglik, n_scoreable, rel_change

Train on **non-functional** reads — out-of-frame *or* stop-codon, since both escaped selection.
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


Likelihood and BIC
------------------

.. code-block:: python

   from vdjtools.model.score import compare_pgen, model_fit, pgen_summary

   model_fit(model, held_out_junctions)            # loglik, k, AIC, BIC
   pgen_summary(compare_pgen(a, b, sequences))     # how differently two models score one set

Two conventions matter. **Likelihoods use nucleotide Pgen**, because ``Σ Pgen_nt`` over all nt CDR3s
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

The bundled ``olga`` model set is a **bit-faithful import** of OLGA's published models: vdjtools
reproduces their marginals exactly, and native Pgen matches ``olga``'s own to machine precision
across all seven loci. Faithful means faithful to the defects too. These are properties of the
source models, verified against OLGA's raw ``model_marginals.txt`` with OLGA's own parser — not
import bugs, and **not** things to "correct" here, because doing so would break the exact-Pgen
invariant that makes the import checkable in the first place.

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
by the collapse, or use the ``learned`` set, which is refit from data on arda germline and does not
inherit the grid. If you need allele resolution *and* an affected gene, be aware Pgen through it is
an underestimate.

Genes with no CDR3-region germline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

OLGA leaves the CDR3-region germline empty for a handful of ORF alleles while still giving them
usage. They can be drawn by the generator but score ``Pgen = 0``, and human ``TRBV23-1`` — 8.6 % of
a real TRB repertoire — is one of them. ``check_model`` reports these as ``unscoreable_gene_mass``.
:func:`~vdjtools.model.io.from_olga` accepts ``derive_orf=True`` to reconstruct the missing germline
from the full germline and the anchor; the ``learned`` models are built that way, while the ``olga``
set keeps it off so it stays an exact Pgen oracle.

V/J usage is protocol-specific
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not a defect, but the most common way to get a wrong answer with these models. OLGA's were fit to
DNA-multiplex data; the ``learned`` set to 5'RACE reads. The two amplify different V genes at very
different rates — TRBV19 is 3.1 % of OLGA's usage and 37 % of these 5'RACE reads — so neither
marginal is right for *your* library. The junction model (trims, insertions, dinucleotides) is the
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

   scored_with = rescale_usage(load_bundled("TRB", "learned"), my_sample)

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

The models shipped in the wheel are fitted to real 5'RACE reads. That whole pipeline —
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

   Real annotated junctions occasionally carry an ambiguous base, and the recombination model is
   defined over A/C/G/T only. Both training entry points substitute ``A`` by default and warn with
   the count (``ambiguous=None`` drops the clonotype instead). Substituting keeps the clonotype —
   one uncertain position in a junction that is otherwise good evidence — and on these reads it
   affects ~0.01% of rows. It is a substitution, not a marginalization, so data with many ambiguous
   positions should use ``ambiguous=None``.


Command line
------------

Every operation above has a ``vdjtools model`` counterpart. A model is named either as a directory
or as ``LOCUS[:source[:organism]]``:

.. code-block:: bash

   vdjtools model check TRB:learned                       # exits 1 on any error
   vdjtools model template --locus TRB -o template/
   vdjtools model learn clones.tsv -t template/ -o fitted/
   vdjtools model log fitted/                              # loglik per iteration
   vdjtools model entropy TRB:olga --table total
   vdjtools model diversity TRB:olga -n 5000
   vdjtools model compare TRB:olga TRB:learned --by gene --dot diff.pdf
   vdjtools model compare-pgen TRB:olga TRB:learned seqs.tsv --summary
   vdjtools model loglik seqs.tsv TRB:learned              # loglik, k, AIC, BIC
   vdjtools model extend fitted/ --locus TRB -o extended/
   vdjtools model rescale TRB:learned my_sample.tsv -o rescaled/
   vdjtools model export TRB:olga --format tsv -o trb_tsv/
   vdjtools model net TRB:olga --format pdf -o bn.pdf
