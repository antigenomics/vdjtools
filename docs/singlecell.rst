Single-cell
===========

:mod:`vdjtools.sc` ingests paired-chain single-cell data, cleans and pairs the chains with
doublet / mispairing QC, scores paired generation probability, and hands the result to the
downstream single-cell ecosystem — scirpy, dandelion, scRepertoire — without glue code.

Install the extra::

   pip install "vdjtools[sc]"

That brings ``anndata``, ``awkward`` and ``mudata``. The downstream tools themselves are
**not** required: ``to_scirpy`` / ``to_dandelion`` raise a clear install hint if you ask for
a container their library builds, while :func:`~vdjtools.sc.from_scirpy` and
:func:`~vdjtools.sc.read_h5ddl` read *their* outputs with only ``awkward`` / ``h5py``. In
other words, consuming someone else's scirpy or dandelion result costs you no new install.

Ingestion
---------

Everything lands in one flat, ``cell_id``-keyed frame — one row per productive contig, with
canonical AIRR Rearrangement columns.

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Source
     - File
     - Reader
   * - CellRanger (preferred)
     - ``airr_rearrangement.tsv``
     - :func:`~vdjtools.sc.read_airr_cell`
   * - CellRanger contigs
     - ``filtered_contig_annotations.csv`` / ``all_contig_annotations.csv``
     - :func:`~vdjtools.sc.read_10x`
   * - arda single-cell
     - ``<prefix>.contigs.airr.tsv`` + ``<prefix>.chains.tsv``
     - :func:`~vdjtools.sc.read_arda_cells`
   * - Anything AIRR + a barcode
     - any Rearrangement TSV with ``cell_id``
     - :func:`~vdjtools.sc.read_airr_cell`

.. code-block:: python

   from vdjtools import sc

   cells = sc.read_airr_cell("outs/airr_rearrangement.tsv")   # CellRanger >= 4.0
   cells = sc.read_10x("outs/filtered_contig_annotations.csv")
   cells = sc.read_arda_cells("run")                          # `arda cells -p run`

Prefer ``airr_rearrangement.tsv`` where you have it: it is the same file scirpy, dandelion
and scRepertoire read, so ingestion and interop agree by construction. ``read_10x`` tolerates
CellRanger version drift — the ``fwr*`` / ``cdr1`` / ``cdr2`` columns exist only from CR6,
``exact_subclonotype_id`` from CR4+, and ``sample`` only under ``cellranger multi``.

.. note::

   A single-cell AIRR table will **not** load through :func:`vdjtools.io.read`. That path
   returns *bulk* clonotype frames, so reading a barcoded table through it would collapse
   reads across cells and drop ``cell_id`` with no error. It raises instead, and points here.
   Pass ``fmt="airr"`` explicitly if pooling into one bulk repertoire is genuinely what you
   want.

Pairing and QC
--------------

.. code-block:: python

   sc.chain_multiplicity(cells)                     # (n_light, n_heavy) -> cell_count
   cells  = sc.resolve_chains(cells)                # top heavy + top light (+ a dual alpha)
   paired = sc.pair_chains(cells, locus_pair="TRA_TRB")
   paired = sc.flag_mispairing(paired)              # ambient / non-canonical alpha chains
   paired = sc.paired_pgen(paired)                  # pgen_alpha, pgen_beta, pgen_paired

``locus_pair`` is one of ``TRA_TRB``, ``TRG_TRD``, ``IGH_IGK``, ``IGH_IGL``.

``paired_pgen`` conditions on the V/J calls, and the model is keyed by **allele**
(``TRBV20-1*01``). CellRanger reports **genes** (``TRBV20-1``), so ``paired_pgen`` resolves
each gene to its representative allele (``*01`` where the model has it) before scoring.

.. warning::

   It does **not** silently fall back to marginalising over all V/J. That fallback once
   returned a Pgen 2.38x too high with no error, which is why
   :func:`vdjtools.model.native.pgen_aa` raises on a gene name in the first place. Pass
   ``resolve_genes=False`` to score only exact allele matches (unmatched rows become null),
   or ``condition_vj=False`` to marginalise deliberately. If an entire locus scores null,
   ``paired_pgen`` warns rather than handing back a silent column of nulls.

.. note::

   :func:`~vdjtools.sc.read_arda_cells` reports arda's own per-chain verdict as
   ``arda_status`` and does not act on it. arda's call and ``resolve_chains``' call are
   independent answers to the same question — compare them rather than assuming they agree.

Interop
-------

Every bridge is built on one thing: a flat AIRR Rearrangement table carrying ``sequence_id``
and ``cell_id``. That is what scirpy, dandelion and scRepertoire all read, and
:func:`~vdjtools.sc.to_airr` is the single place the spelling differences are reconciled
(AIRR says ``junction``, vdjtools says ``junction_nt``; scRepertoire reads
``consensus_count`` where the others prefer ``umi_count``, so both are emitted).

.. code-block:: python

   airr = sc.to_airr(cells)                    # the universal contract
   sc.write_airr(cells, "airr_rearrangement.tsv")

scirpy / scverse
~~~~~~~~~~~~~~~~

.. code-block:: python

   adata = sc.to_scirpy(cells)                 # obsm["airr"] awkward layout + chain_indices
   mdata = sc.to_scirpy(cells, gex=gex_adata)  # MuData({"gex": ..., "airr": ...})
   back  = sc.from_scirpy(adata)               # -> canonical frame (needs only awkward)

``to_scirpy`` delegates to ``scirpy.io.read_airr``, so the object carries scirpy's own
``obsm["airr"]`` layout exactly as their version defines it rather than a copy of it that
could drift. ``index_chains=True`` (the default) also runs ``scirpy.pp.index_chains``, which
every scirpy tool needs.

dandelion
~~~~~~~~~

.. code-block:: python

   vdj  = sc.to_dandelion(cells)               # .data contigs, .metadata cells
   back = sc.from_dandelion(vdj)
   back = sc.read_h5ddl("dandelion_data.h5ddl")   # no dandelion install needed

scRepertoire (R)
~~~~~~~~~~~~~~~~

Export only — no R code ships with vdjtools.

.. code-block:: python

   sc.write_screpertoire(cells, "airr_rearrangement.tsv")            # loadContigs(format="AIRR")
   sc.write_screpertoire(cells, "contigs.csv", format="10x")         # loadContigs(format="10X")

.. code-block:: r

   contigs  <- loadContigs("airr_rearrangement.tsv", format = "AIRR")
   combined <- combineTCR(contigs)
   seurat   <- combineExpression(combined, seurat)

.. warning::

   ``combineTCR(samples=, ID=)`` rewrites barcodes to ``sample_ID_barcode``. That is the
   usual cause of a silent barcode-join failure against a Seurat ``meta.data``: pass the same
   ``samples`` / ``ID`` to ``combineExpression``, or omit both.

Augmenting an existing object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To push something vdjtools computed onto an object the rest of an analysis already uses,
:func:`~vdjtools.sc.push_obs` writes into an ``AnnData.obs`` or a ``Dandelion.metadata``:

.. code-block:: python

   scores = sc.paired_pgen(paired).select(["cell_id", "pgen_paired"])
   sc.push_obs(adata, scores)                  # also works on a Dandelion

It needs one row per cell; a frame where a cell contributed two receptor pairs must be
aggregated first (it raises rather than picking one silently).

Two AnnData shapes
~~~~~~~~~~~~~~~~~~

:func:`~vdjtools.sc.to_scirpy` is the scverse-native container — one ``obs`` row per **cell**,
full per-contig AIRR records under ``obsm["airr"]``. :func:`~vdjtools.sc.to_anndata` is a flat
alternative — one ``obs`` row per receptor **pair**, ``alpha_*`` / ``beta_*`` columns — handy
for a quick paired table or for attaching an expression matrix, but scirpy does not read it.
Reach for ``to_scirpy`` unless you specifically want the flat view.

.. note::

   A **bulk cohort** must not go into AnnData: ``obs=clonotype`` over many samples yields an
   enormous, almost-empty sparse matrix. Use :func:`vdjtools.io.scan_cohort` instead. Rule of
   thumb: single-cell (``obs=cell``) to AnnData, bulk cohort to Parquet.

Command line
------------

.. code-block:: bash

   vdjtools sc convert outs/airr_rearrangement.tsv -o cells.tsv
   vdjtools sc convert outs/airr_rearrangement.tsv --airr -o airr.tsv
   vdjtools sc qc      outs/airr_rearrangement.tsv
   vdjtools sc pair    outs/airr_rearrangement.tsv --flag-mispairing -o paired.tsv
   vdjtools sc pgen    outs/airr_rearrangement.tsv -o pgen.tsv
   vdjtools sc export  outs/airr_rearrangement.tsv --to scirpy       -o vdj.h5ad
   vdjtools sc export  outs/airr_rearrangement.tsv --to dandelion    -o vdj.h5ddl
   vdjtools sc export  outs/airr_rearrangement.tsv --to screpertoire -o airr.tsv

Worked example
--------------

``examples/single_cell_interop.py`` is a marimo notebook running the whole path —
CellRanger in, QC and paired Pgen, then out to scirpy, dandelion and scRepertoire, and back::

   pip install "vdjtools[examples,sc]"
   marimo edit examples/single_cell_interop.py

``examples/single_cell.py`` covers the analysis side: dCODE dextramer data, antigen-driven
clustering, and grading it with :func:`~vdjtools.sc.cluster_eval`.
