Raw features and sample PCA
===========================

``vdjtools features`` exports the unrotated statistical feature table and its
channels, without a fitted corpus. It includes diversity, isotype and segment/kmer
usage through the existing signature feature contract. Declared log, CLR and logit
transforms remain; raw here means before corpus scaling and PCA.

.. code-block:: bash

   vdjtools features samples/*.tsv --threads 2 --coverage 0.8 -o features.tsv
   vdjtools pca features.tsv --columns feature_columns.txt --fit rotation.npz --components 32 -o pcs.tsv
   vdjtools pca new_features.tsv --columns feature_columns.txt --rotation rotation.npz -o new_pcs.tsv

``--coverage min`` (default) uses the cohort's common attainable coverage; ``own``
uses each sample's observed coverage. Use an explicit shared coverage level when
comparing independently generated tables. ``--columns`` lists the feature columns
to fit, one per line; exclude IDs and pass-through channels. ``--cohort TABLE``
fits the rotation using a compatible reference table. The ``embedding_space``
column must identify the same feature basis in both tables; ``--space ID`` can
supply that identity for an external table lacking it. Incompatible identities
or missing fitted columns raise. IDs and unselected channels pass through unchanged.

The shared implementation is ``vdjtools.signature.corpus.TablePCA``. Its saved
NumPy artifact stores imputation values, scaling, centering, feature order and PCA
components. This optional StandardScaler/PCA workflow is distinct from the existing
robust, winsorized ``Corpus`` rotation used by ``vsig`` and mirpy ``rsig``. Existing
rsig is not RFF. Raw feature rotations and mirpy sample RFF rotations are separate outputs.
Raw outputs remain available independently of either rotation.

Choosing feature columns
------------------------

Select measurements explicitly for statistical or signature tables. Exclude
``sample_id``, ``embedding_space`` and pass-through channels from the fit; selected
columns are replaced by ``pc0``, ``pc1``, and so on. Unselected columns retain their
values and input order. For an RFF table, ``phi:`` columns are selected by default.

.. code-block:: python

   import vdjtools                       # establish runtime limits before Polars/NumPy
   from functools import partial
   from vdjtools.io import read
   from vdjtools.signature import raw_features_cohort, raw_columns
   from vdjtools.signature.features import gene_vocab
   from vdjtools.signature.layout import LOCI, NO_LOCUS
   from vdjtools.signature.corpus import TablePCA

   readers = [("S1", partial(read, "s1.tsv", keep=("v_identity",))),
              ("S2", partial(read, "s2.tsv", keep=("v_identity",)))]
   vocab = {loc: gene_vocab(loc, "human") for loc in LOCI}
   table = raw_features_cohort(readers, vocab=vocab, n_jobs=1, cstar_target=0.8)
   columns = [c for loc in (*LOCI, NO_LOCUS)
              for c in raw_columns("vsig", loc, vocab.get(loc))]
   rotation = TablePCA.fit(table, columns=columns, n_components=2)
   rotation.save("feature_pca.npz")
   pcs = rotation.transform(table)

To use compatible control rows, fit on the control table and transform the target:

.. code-block:: bash

   vdjtools pca features.tsv --cohort control_features.tsv --columns feature_columns.txt \
       --fit control_rotation.npz --components 32 -o control_pcs.tsv

Fit imputes each selected column with its training median; an entirely missing
training column uses zero. Apply always uses those same imputation values.
A saved rotation validates feature names and ``embedding_space`` before projection.
An identity supplied with ``--space`` is your declaration of compatibility, not a
conversion between different vocabularies, coverage targets, weights or corpora.

The same table PCA can project vsig, rsig or mirpy RFF exports. Those remain three
different representations. Existing corpus signatures need explicit feature
selection and a truthful corpus identity when projected again. See :doc:`signature`
for the robust corpus transform and `mirpy's embedding guide
<https://docs.isalgo.dev/mirpy/embeddings.html>`_ for RFF bases.

Resource limits
---------------

Readers run in bounded sample workers with serial inner kernels. Raw feature rows
are preallocated; input frames are released after each row. The final table still
grows with sample count and selected feature width. Table PCA needs a dense fitting
matrix and its numerical workspace; known large allocations check available
container memory before allocation. These checks cannot bound an arbitrarily deep
input frame. See :ref:`cpu-container-limits` for CPU allocation and import order.

Reference rotation bundles
--------------------------

A reference bundle's ``features-pca.npz`` projects the unrotated statistical
measurements, using its saved feature order and the same shared coverage target.
The fitting recipe uses human gene vocabularies and ``log2p1`` weighting; its
default coverage is ``0.8``. The RFF rotation in that bundle is a different file
and requires the bundled mirpy bases. Neither file is a vsig/rsig corpus.

.. code-block:: bash

   vdjtools features samples/*.tsv --coverage 0.8 --threads 2 -o features.tsv
   vdjtools pca features.tsv --rotation reference/features-pca.npz -o pcs.tsv

The `mirpy reference recipe
<https://docs.isalgo.dev/mirpy/embeddings.html#reference-rotation-bundles>`_
freshly generates each reference table, fits with the existing ``TablePCA``,
checks save/load equivalence, and exports rotations and coordinate parameters.
Input membership and fit diagnostics remain private. ``bundle.json`` records
parameters and file checksums; it contains no sample identities or raw input paths.

Extract a reference bundle separately from the library installation and verify
its recorded hashes before use. Choose the blood/tissue RNA-seq, synthetic or
deep TRA/TRB reference appropriate to your analysis; keep the chosen rotation
checksum with your results. The deep reference selects TRA/TRB measurements and
cross-locus channels rather than fitting the other loci. Unselected measurements
and channels pass through the projection.

The first command keeps the unprojected feature table. Applying a rotation does
not recompute those measurements or turn them into a corpus signature. Target
and control tables must use the same vocabulary, weighting and coverage target.
A shared vocabulary alone does not make PCs from different rotations comparable.
