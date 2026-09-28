Worked examples
===============

Fifteen runnable notebooks under ``examples/``, each a complete piece of work --- its data, its
commands and its conclusions in one file. Several reproduce a published result, so you can see what
the library does on data whose answer is already known.

They are `marimo <https://marimo.io>`_ notebooks, which are **plain Python files** with
``@app.cell`` decorators rather than JSON. That means they diff and review like source code, and
they run three ways:

.. code-block:: bash

   pip install 'vdjtools[examples]'
   marimo edit examples/preprocess.py     # interactive, cells re-run as you edit
   marimo run  examples/preprocess.py     # read-only app, no code shown
   python      examples/preprocess.py     # plain script; prints, no UI

Data is fetched from Hugging Face on first run and cached under ``examples/.data/``, so a fresh
``pip install`` is enough --- no local paths and no pre-staged files. If you already hold a copy,
drop it under ``./data_dump/`` (gitignored) and it is used instead of downloading. Notebooks that
need an extra say so; ``examples/README.md`` gives the exact install line for each.

Start here
----------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Pre-processing
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/preprocess.py

      What to do to real data before measuring anything. The three independent filtering axes ---
      productive, abundance, and sequencing artefact --- on real Britanova samples from three
      sequencing batches, with what each filter removes and what it costs you.

   .. grid-item-card:: Ageing of the TCR repertoire
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/aging.py

      A complete study, read three ways on one cohort: 78 donors aged 0--103. Diversity falls with
      age, repertoires become more clonal, and individual repertoires drift apart --- each measured
      in a way that removes the sequencing-depth confound rather than assuming it away.

Describing a cohort
-------------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Cohorts too large for memory
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/scale_cohort.py

      The out-of-core pattern: every sample streamed once into a partitioned Parquet store, then
      every statistic a single streamed query over it. Slide the sample count toward thousands and
      watch peak memory stay flat. No data to download.

   .. grid-item-card:: Junction physicochemistry
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/cdr_features.py

      Amino-acid properties of the junction --- hydropathy, charge, volume, the Kidera factors ---
      per sample, correlated with donor age across the ageing cohort. Junction hydropathy carries
      the strongest single association.

Comparing repertoires
---------------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Overlap and similarity
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/overlap_similarity.py

      Three different questions that all get called "overlap": identical clonotypes, clonotypes
      within one substitution, and graded sequence similarity. The third finds shared structure
      between repertoires that share **no** identical clonotype, which the first two report as zero.

   .. grid-item-card:: Longitudinal tracking
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/vaccination_tracking.py

      Following clonotypes through a vaccination time course --- yellow fever, influenza, TBE.
      Which clones emerged, expanded, persisted, contracted or vanished, tested within donor, plus
      a recapture model for the clones you did not resample.

Finding sequences associated with a phenotype
---------------------------------------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: CMV and HLA association
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/emerson_biomarker.py

      The core of Emerson et al. (*Nat Genet* 2017) on their published cohort: an incidence-based
      screen for public TCR-beta chains associated with CMV serostatus or HLA-A\*02, checked live
      against VDJdb.

   .. grid-item-card:: The same screen at cohort scale
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/emerson_cmv_hla.py

      The streaming version of that screen over the whole cohort, with peak memory reported, for
      when the candidate set is the whole repertoire rather than a shortlist.

   .. grid-item-card:: Association explorer
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/biomarker_explorer.py

      The interactive superset: condition, test and match scope as controls, with a co-occurrence
      panel for finding chains that travel together --- alpha/beta partners, or one clone's
      relatives.

   .. grid-item-card:: A disease motif, and its confounder
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/ankspond_motif.py

      Reproducing a published TRBV9 junction motif in ankylosing spondylitis, and why it takes an
      HLA-B27-matched comparison to say anything: B27 is 26/27 tied to disease in that cohort, so
      the unmatched contrast cannot separate disease from carriage.

Features for a model
--------------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Signature features
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/signature_features.py

      Three feature choices that are easy to get wrong, each measurable in seconds: why the
      amino-acid block is arcsine-transformed rather than logged, why an ambiguous V call must be
      resolved rather than stripped, and why the number of components kept is a real choice.

   .. grid-item-card:: Signature pipeline (mirpy)
      :link: https://github.com/antigenomics/mirpy/blob/master/examples/signature_pipeline.py

      A folder of AIRR files to one row per sample, joined to your metadata sheet --- both halves
      of the signature, the commands and the join. The end-to-end recipe; lives in the companion
      repository.

Single-cell data
----------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Paired-chain single cell
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/single_cell.py

      10x contigs in, alpha/beta receptors out: chain-multiplicity QC, doublet handling, pairing,
      then unsupervised clustering of beta junctions graded against dextramer antigen labels.

   .. grid-item-card:: scirpy, dandelion, scRepertoire
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/single_cell_interop.py

      vdjtools sitting in the middle of the single-cell ecosystem rather than replacing any of it.
      One flat AIRR table is the seam: hand it to each container, read it back, and attach a score
      to an AnnData you did not build.

Recombination models
--------------------

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Model workshop
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/model_workshop.py

      Build, fit, check, compare and score a V(D)J generation model end to end --- from a custom
      germline library through expectation-maximisation to held-out likelihood. Runs offline in
      seconds on a toy locus.

   .. grid-item-card:: Model explorer
      :link: https://github.com/antigenomics/vdjtools/blob/master/examples/model_explorer.py

      What is actually inside a recombination model: the dependency graph between events,
      per-event entropy, mutual information, and every marginal table. Uses the models shipped in
      the wheel, so there is nothing to download.
