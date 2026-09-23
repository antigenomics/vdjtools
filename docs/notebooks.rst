Notebooks
=========

Worked examples as `marimo <https://marimo.io>`_ notebooks --- **plain Python files** with
``@app.cell`` decorators, so they diff and review like source rather than like JSON, and they run
three ways:

.. code-block:: bash

   pip install 'vdjtools[examples]'
   marimo edit examples/signature_features.py    # interactive
   marimo run  examples/signature_features.py    # read-only app
   python      examples/signature_features.py    # plain script; prints, no UI

Each one bootstraps its own data from Hugging Face
(`isalgo/airr_benchmark <https://huggingface.co/datasets/isalgo/airr_benchmark>`_) and caches it, so
a fresh ``pip install`` is enough --- no local paths and no pre-staged files. Drop a copy under
``./data_dump/`` (gitignored) and it is used instead of downloading.

Signatures
----------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="https://github.com/antigenomics/vdjtools/blob/master/examples/signature_features.py">
       <h3>Signature features</h3>
       <p>The <code>vsig</code> half block by block: what each column measures, how it is
       transformed, and which ones survive a shallow library.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/mirpy/blob/master/examples/signature_pipeline.py">
       <h3>Signature pipeline (mirpy)</h3>
       <p>A folder of AIRR TSVs to one table joined with your metadata, in one command --- both
       halves at once. The end-to-end recipe.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/vdjtools/blob/master/examples/scale_cohort.py">
       <h3>Scaling a cohort</h3>
       <p>Coverage-standardised diversity across samples of wildly different depth, and why the
       standardisation is what makes them comparable.</p>
     </a>
   </div>

Analysis
--------

.. raw:: html

   <div class="proj-card-grid">
     <a class="proj-card" href="https://github.com/antigenomics/vdjtools/blob/master/examples/aging.py">
       <h3>Aging</h3>
       <p>78 donors aged 0-103, three ways: streamed cohort statistics, coverage-standardised
       iNEXT diversity, and repertoire divergence by MDS.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/vdjtools/blob/master/examples/preprocess.py">
       <h3>Pre-processing</h3>
       <p>The three filtering axes, and what each one does to a real cohort.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/vdjtools/blob/master/examples/overlap_similarity.py">
       <h3>Overlap and similarity</h3>
       <p>Pairwise repertoire overlap, the metrics that behave under depth variation, and the ones
       that do not.</p>
     </a>
     <a class="proj-card" href="https://github.com/antigenomics/vdjtools/blob/master/examples/model_workshop.py">
       <h3>Model workshop</h3>
       <p>Build, check and compare V(D)J generation models; score sequences; read the training log.</p>
     </a>
   </div>
