Command reference
=================

``pip install vdjtools`` installs one ``vdjtools`` command. This page lists every command, what it
is for, and the options that change an answer rather than only a path.

``vdjtools <command> --help`` is the authoritative list of flags and defaults -- it is generated from
the code, so it cannot drift. This page is the map.

Conventions
-----------

These hold for every command, so they are stated once here rather than repeated per command.

**Output.** ``-o/--out`` writes TSV, or Parquet when the path ends in ``.parquet`` or ``.pq``. Omit
it and the result goes to stdout while progress and errors go to stderr, so commands pipe cleanly.

**Input.** Any supported format is detected: native vdjtools, AIRR Rearrangement TSV, Parquet, MiXcr,
MiGec, Adaptive immunoSEQ, IMGT/HighV-QUEST, Vidjil, RTCR, TRUST4, arda. ``--format`` overrides the
sniffer when you need it to.

**Cohorts.** Analysis commands take either a list of sample files or a metadata table with
``-m/--metadata`` plus ``--base-dir``, mirroring the legacy tool's workflow. ``--sample-col`` and
``--file-template`` control how a metadata row becomes a filename.

**Parallelism.** ``-t/--threads`` on analysis commands parallelises over samples; ``0`` means all
cores. ``--cohort`` instead streams one pass over a pre-ingested Parquet store, which is what to use
when the cohort does not fit in memory. On ``signature``, the flag is ``-j/--jobs`` and it means
worker **processes** -- the per-sample work is polars and numpy rather than a GIL-releasing kernel,
so processes are the layer that helps.

**Models.** Wherever a command takes a model, it accepts ``LOCUS``, ``LOCUS:source``,
``LOCUS:source:organism``, or a path to a model directory -- for example ``TRB``, ``TRB:learned``,
``TRB:olga:human``.

Data and pre-processing
-----------------------

``convert`` -- any format to the canonical table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools convert mixcr.txt.gz -o clones.parquet

Read anything supported and write the canonical AIRR-junction table. Converting a cohort to Parquet
once is usually the cheapest thing you can do to a repeated analysis.

``filter`` -- three separate axes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools filter clones.parquet --productive --min-freq 1e-4 -o productive.tsv
   vdjtools filter clones.parquet --nonproductive -o nonproductive.tsv
   vdjtools filter clones.parquet --functional-genes --keep-orf -o f_orf.tsv

The three filtering axes are deliberately separate flags, because they are separate facts:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Flag
     - Asks
   * - ``--productive``
     - Does the **rearrangement** encode a chain? In frame, no stop codon.
   * - ``--functional-genes``
     - Is the **germline gene** real? IMGT F, ORF or P; ``--keep-orf`` admits ORF.
   * - ``--min-len`` / ``--max-len``
     - Is ``junction_aa`` within an inclusive length bound? Defaults sanity-bound 5 to 60.

A productive rearrangement can use a pseudogene V, so filtering one axis says nothing about the
other. ``--recompute-frequencies`` is on by default, renormalising ``frequency`` over the survivors;
``--keep-frequencies`` leaves the file's own values untouched.
``--min-freq``, ``--v``, ``--j`` and ``--remove`` add frequency and segment selection.

``downsample`` -- normalise depth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools downsample clones.parquet 100000 -o ds.tsv
   vdjtools downsample clones.parquet 5000 --clones -o ds.tsv     # unique clonotypes, not reads

Reduces to a common depth by resampling. ``--seed`` makes it reproducible. Prefer
coverage-standardised diversity over downsampling where the question allows it -- see
:doc:`preprocessing`.

``correct-vj`` -- cross-batch usage bias
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools correct-vj s1.tsv s2.tsv s3.tsv s4.tsv -b A,A,B,B \
       --transform sigmoid --usage-out usage.tsv --outdir corrected/

Corrects V/J usage differences between batches and, with ``--outdir``, rewrites the clonotype tables
to match. ``-b/--batches`` is required and parallel to the sample list. ``--transform`` chooses
``location`` (default) or ``sigmoid``; ``--scope`` selects ``v``, ``j`` or ``vj``.

Four more options change the answer rather than a path. ``--winsor-q`` clamps the per-batch mean and
sigma at a quantile --- off by default, matching the published method, with ``0.025`` a robustness
setting for shallow or RNA-seq repertoires. ``--unweighted`` builds usage from distinct clonotypes
instead of reads, which is the right choice when one hyperexpanded clone would otherwise define a
batch's usage. ``--z-cap`` bounds ``|Z|`` under ``--transform sigmoid``. ``--rescale`` rewrites
counts deterministically instead of resampling them, so the output is reproducible without a seed
but no longer integer-sampled. See :doc:`preprocessing`.

``pool`` -- combine samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools pool s1.tsv s2.tsv s3.tsv -o pooled.tsv
   vdjtools pool s1.tsv s2.tsv s3.tsv --join --min-samples 2 -o joint.tsv

Without ``--join``, sums counts into one deeper repertoire. With it, keeps clonotypes seen in at
least ``--min-samples`` samples -- an incidence join, which is the right input to public-clonotype
work. ``--key`` selects nucleotide or amino-acid identity.

Repertoire analytics
--------------------

All four take sample files or ``-m``, and accept ``-t`` and ``--cohort``.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Command
     - What it computes
   * - ``diversity``
     - Observed richness, Chao1, Chao coverage, Efron-Thisted, Shannon, normalised Shannon, inverse
       Simpson, d50 -- one row per sample. ``--on-duplicate sum`` merges rows that repeat a
       clonotype key; the default refuses, because a table with no ``junction_nt`` cannot say
       whether two such rows are one clonotype or two.
   * - ``spectratype``
     - Junction-length distribution. ``--kind aa|nt``, ``--weight reads|clones``.
   * - ``segment-usage``
     - V, D, J, C or VJ usage. ``--segment v`` by default, ``--weight reads|clones``.
   * - ``overlap``
     - Exact-match pairwise overlap (D, F, F2, R) for every sample pair.

.. code-block:: bash

   vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
   vdjtools spectratype   --cohort cohort_parquet/ -o spectra.tsv
   vdjtools segment-usage -m metadata.txt --base-dir samples/ -t 0 -o usage.tsv
   vdjtools overlap       *.tsv -o overlap.tsv

Enrichment and longitudinal
---------------------------

``tcrnet`` and ``alice`` -- neighbourhood enrichment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both ask whether a clonotype has more similar neighbours than expected, and they differ in what
supplies the expectation:

.. code-block:: bash

   vdjtools tcrnet sample.tsv --scope 1,0,0,1 -o tcrnet.tsv   # vs a CONTROL REPERTOIRE
   vdjtools alice  sample.tsv --scope 1,0,0,1 -o alice.tsv    # vs a GENERATION MODEL

``--scope`` is the edit-distance ball as ``substitutions,insertions,deletions,total``. ``alice`` adds
``--q`` (the significance threshold), ``--min-degree`` and ``--min-count``; both take ``--locus``,
``--species`` or ``--source``, and ``--threads``.

``dynamics`` -- what changed between two timepoints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools dynamics day0.tsv day15.tsv -o tracked.tsv

A paired within-donor test classifying each clonotype as emergent, expanded, persistent, contracted
or vanishing. ``--neff`` overrides the effective depth used for the pair, ``--umi`` declares that
counts are molecular barcodes, ``--min-total`` and ``--alpha`` set the testing thresholds.

Signatures
----------

``signature`` -- one row of named features per sample
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools signature --corpus blood samples/*.tsv.gz -o vsig.tsv
   vdjtools signature --corpus blood --components 32 --describe
   vdjtools signature --corpus tissue -m metadata.txt --base-dir samples/ -j 8 -o vsig.tsv

``--corpus`` is **required** and takes a published name or a path; a silently chosen rotation would
make two matrices look comparable when they are not. Key options:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Option
     - Meaning
   * - ``--corpus``
     - ``blood``, ``tissue``, ``deep-tcr``, ``blood-uncapped``, ``tissue-uncapped``,
       ``synthetic-blood``, ``synthetic-tissue``, ``naive``, ``memory``, or a path. A published name
       is downloaded and cached on first use.
   * - ``--components``
     - An integer count (``128``) or a variance fraction (``0.95``). Truncating a wider rotation
       downward is exact; asking for more than was fitted refuses.
   * - ``--winsorize``
     - ``features`` (clamp, then rotate), ``pcs`` (rotate, then clamp) or ``none``. Whatever is
       clamped is reported in ``vsig:qc:-:winsor_frac``.
   * - ``--cstar-target``
     - Coverage level the Hill numbers are read at: a number, ``min`` (default, the cohort's own
       per-locus minimum) or ``own`` (each sample's own, not comparable across samples).
   * - ``--winsor-p``
     - Which stored percentile to clamp at, ``0.01`` or ``0.05``. Both are stored in the artifact,
       so switching needs no refit; the default is whichever it was fitted with.
   * - ``--named``
     - Also emit the raw blocks the rotation is fitted on, under their own names:
       ``div``, ``depth``, ``clon``, ``len``, ``iso``, ``shm``, ``pair``, or ``all``. Off by
       default. The values carry their **declared transform**, not a natural scale --- a ``log10``
       diversity comes back as a ``log10``. See :ref:`sig-named`.
   * - ``--on-duplicate``
     - ``error`` (default) or ``sum``. A frame with no ``junction_nt`` that repeats
       ``(junction_aa, v_call, j_call, c_call)`` cannot say whether those rows are one clonotype or
       two, so the library refuses rather than guessing.
   * - ``--describe``
     - Print the columns **this** invocation emits, and exit. Resolved against ``--corpus``,
       ``--components`` and ``--columns``, so it is the exact header you will get.

``corpus`` -- fetch or build an artifact
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools corpus --fetch all                              # pre-warm the cache
   vdjtools corpus --corpus synthetic-blood -o sb.npz       # build one yourself
   vdjtools corpus --smoke -o /tmp/smoke.npz                # minutes, not hours

Building uses no samples from anybody's cohort: every receptor is drawn from the bundled
recombination models, so the artifact is reproducible by anyone who installs the library, and is
byte-identical across processes and thread counts. ``--samples``, ``--size``, ``--seed``, ``--loci``
and ``--depth-spread`` control the draw, ``--components`` and ``--winsor-p`` control the fit, and all
are recorded in the manifest. See :doc:`signature`.

Recombination models
--------------------

``pgen``, ``generate``, ``models``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   vdjtools models                                 # list the bundled models
   vdjtools generate -m TRB -n 1000 -o gen.tsv     # cf. olga-generate_sequences
   vdjtools pgen seqs.tsv -m TRB -o pgen.tsv       # cf. olga-compute_pgen
   vdjtools pgen seqs.tsv -m TRB --mismatches 1    # and the Hamming-1 ball

``pgen`` takes ``--column`` for the sequence column and ``--v-col`` / ``--j-col`` to condition on the
observed V and J; ``--type`` selects nucleotide or amino acid, ``auto`` by default, and
``--no-header`` writes the values alone for piping. ``generate``
takes ``-n/--number`` for how many sequences to draw, ``--productive``/``--no-productive`` to keep
only productive draws (on by default, at roughly 3.6 times the cost per kept sequence on TRB), and
``--seed`` for reproducibility. Both commands take ``-m/--model`` for a bundled model, or
``--model-path`` for one on disk.

``model`` -- the workshop
^^^^^^^^^^^^^^^^^^^^^^^^^^

Thirteen subcommands for building, fitting, auditing and comparing models. Full walkthrough in
:doc:`model`.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Subcommand
     - What it does
   * - ``list``
     - List the bundled models, as ``vdjtools models`` does.
   * - ``template``
     - Scaffold a model from a germline library -- arda's, or your own FASTA plus anchors.
   * - ``learn``
     - Fit marginals from your own sequences by EM, writing the training log alongside.
       ``--checkpoint`` and ``--resume`` survive a long run.
   * - ``build``
     - Build several chains from the full AIRR read corpus: fetch, arda-map, then EM.
   * - ``check``
     - Audit a model against its manifest, its germline and a reference library. Exits non-zero on
       an error, so it belongs in CI.
   * - ``log``
     - Show the EM training log, log-likelihood per iteration, one block per run.
   * - ``extend``
     - Add alleles from a larger germline library, seeded from what the model already knows.
   * - ``rescale``
     - Replace a model's V/J usage with your own sample's, keeping its junction model.
   * - ``export``
     - Export every probability as tables: a hand-editable directory, or one long frame.
   * - ``entropy``
     - Information content per recombination event: entropy, mutual information, or the total.
   * - ``diversity``
     - Total diversity: scenario entropy, sequence entropy, effective diversity.
   * - ``net``
     - Render the recombination Bayes net, nodes annotated with entropy and edges with mutual
       information.
   * - ``compare`` / ``compare-pgen``
     - Compare two models parameter by parameter, or score one sequence set under both and compare
       the Pgen distributions.
   * - ``loglik``
     - Log-likelihood of a sequence set under a model, with free parameters, AIC and BIC.

.. code-block:: bash

   vdjtools model check TRB:learned
   vdjtools model template --locus TRB -o tmpl/
   vdjtools model learn clones.tsv -t tmpl/ -o fitted/
   vdjtools model compare TRB:olga TRB:learned --by gene --dot diff.pdf
   vdjtools model loglik seqs.tsv TRB:learned

Single cell
-----------

``sc`` -- five subcommands
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Subcommand
     - What it does
   * - ``convert``
     - Read a single-cell contig table into the canonical long frame, or AIRR with ``--airr``.
   * - ``pair``
     - Resolve each cell's chains and emit one row per paired receptor.
       ``--flag-mispairing`` reports rather than drops.
   * - ``qc``
     - Chain-multiplicity quadrants: how many cells carry how many light by heavy chains.
   * - ``pgen``
     - Paired generation probability per cell, Pgen(alpha) times Pgen(beta).
   * - ``export``
     - Export for scirpy, dandelion, scRepertoire or any AIRR consumer.

.. code-block:: bash

   vdjtools sc convert filtered_contig_annotations.csv -o long.tsv
   vdjtools sc pair    filtered_contig_annotations.csv --flag-mispairing -o paired.tsv
   vdjtools sc pgen    filtered_contig_annotations.csv -o paired_pgen.tsv
   vdjtools sc export  filtered_contig_annotations.csv --to scirpy -o adata.h5ad

``--locus-pair`` defaults to ``TRA_TRB``; set it for BCR. Details and the interop matrix are in
:doc:`singlecell`.
