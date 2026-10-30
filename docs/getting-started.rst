Getting started
===============

A first pass through vdjtools in about ten minutes: install it, get a repertoire to work on,
measure it, compare two samples, and reduce a sample to one feature vector. The synthetic
examples run offline; the signature artifact is downloaded on first use.

If you are looking for a specific task rather than an introduction, the :doc:`guides <usage>` are
organised by task, and :doc:`cli` lists every command.

Install
-------

.. code-block:: bash

   pip install vdjtools
   vdjtools --help

Wheels are prebuilt for CPython 3.10-3.13 on Linux, macOS (Apple Silicon) and Windows, and include
the native C++ extension. A source install needs a C++ toolchain (Xcode Command Line Tools on macOS,
``build-essential`` on Linux).

.. _gs-sample:

Step 1 -- a repertoire to work on
---------------------------------

Normally you would load your own data with one call, and vdjtools would detect the format:

.. code-block:: python

   from vdjtools import io as vio

   sample = vio.read("clones.tsv")     # MiXcr, immunoSEQ, AIRR, Parquet, native -- detected

So that this page runs with nothing downloaded, draw a synthetic repertoire from one of the
recombination models that ship in the wheel. This fixed tutorial explicitly uses
the historical OLGA model to preserve the saved example outputs below; ordinary
model calls now default to ``isalgo`` (see :doc:`model` for the candidate release
gates). Generated receptor sequences have a rank-based clone-size distribution:

.. code-block:: python

   import vdjtools                       # establish runtime limits first
   import numpy as np
   import polars as pl
   from vdjtools import io as vio
   from vdjtools.io import schema as S
   from vdjtools.model import load_bundled
   from vdjtools.model.generate import generate

   def demo_sample(locus="TRB", clones=2000, reads=20_000, a=1.0, seed=0):
       """A synthetic repertoire: sequences from a bundled model, clone sizes from a rank law."""
       seqs = generate(load_bundled(locus, "olga-legacy"), clones, seed=seed, productive_only=True)
       rng = np.random.default_rng(seed)
       freq = np.arange(1, seqs.height + 1, dtype=float) ** -a      # f_i proportional to i**-a
       counts = rng.multinomial(reads, freq / freq.sum())
       out = seqs.with_columns(pl.Series(S.COUNT, counts)).filter(pl.col(S.COUNT) > 0)
       return vio.normalize(out, recompute_freq=True)

   sample = demo_sample(seed=1)
   sample.select(["v_call", "j_call", "junction_aa", "duplicate_count", "frequency"]).head(3)

.. code-block:: text

   shape: (3, 5)
   | v_call      | j_call     | junction_aa        | duplicate_count | frequency |
   | TRBV7-2*01  | TRBJ2-2*01 | CASRWRNNGTFTGELFF  |            2420 |    0.1210 |
   | TRBV24-1*01 | TRBJ1-5*01 | CATSDSAGTGVGNQPQHF |            1172 |    0.0586 |
   | TRBV5-1*01  | TRBJ2-7*01 | CASSQGGPAYEQYF     |             820 |    0.0410 |

This is the **canonical clonotype frame**, and it is the only table vdjtools passes around: a polars
``DataFrame`` on AIRR column names, one row per clonotype. Whatever format you read, you get these
columns, and repertoire analyses consume this schema and return Polars result tables --
which lets results feed into the next step without a format converter.

.. note::

   vdjtools uses ``junction_nt`` / ``junction_aa``; AIRR spells the nucleotide field
   ``junction``. A **junction** includes the conserved Cys104 and Phe/Trp118 anchor residues;
   an IMGT **CDR3** excludes them, so
   ``junction_aa`` is two residues longer than the CDR3 for the same receptor. This matters before
   any Pgen, clustering or database-matching work -- see :doc:`glossary`.

Step 2 -- measure one repertoire
--------------------------------

.. code-block:: python

   from vdjtools import stats

   stats.diversity_stats(sample)

.. code-block:: text

   reads   observed_diversity   chao1     efron_thisted   shannon_wiener   normalized_shannon_wiener   inverse_simpson   d50
   20000                 1790   1978.04          1977.0          268.463                      0.7467           41.5159   0.9810

Reading that row: 20,000 reads resolved into 1,790 distinct clonotypes; Chao1 estimates 1,978
including what was missed, so this sample is close to saturated. ``inverse_simpson`` of 41.5 is the
effective number of equally-abundant clones -- far below the 1,790 observed, because the
distribution is skewed. ``d50`` of 0.981 is the fraction of clonotypes *not* needed to reach half
the reads, so 1.9 percent of clones (34 of them) carry half the sample.

Richness depends on how deeply you sequenced, so two samples at different depths are not directly
comparable. Compare them at a matched **coverage** instead:

.. code-block:: python

   stats.inext(sample, q=(0, 1, 2))     # Hill numbers with bootstrap intervals

Gene usage and the junction-length spectrum:

.. code-block:: python

   stats.segment_usage(sample, "v")     # or "j", "d", "c"
   stats.vj_usage(sample)                # joint V/J usage
   stats.spectratype(sample)            # junction length distribution

Step 3 -- compare two repertoires
----------------------------------

.. code-block:: python

   from vdjtools import overlap

   other = demo_sample(seed=2)
   overlap.overlap_metrics(sample, other)

.. code-block:: text

   {'D': 0.0, 'F': 0.0, 'F2': 0.0, 'R': None, 'd1': 1790, 'd2': 1763, 'd12': 0}

Zero shared clonotypes, which is the correct answer and a useful calibration: two repertoires drawn
independently from the same recombination model overlap not at all on exact nucleotide match. So a
non-zero overlap in real data is evidence of something -- shared exposure, relatedness between
donors, or cross-sample contamination -- rather than a baseline you have to subtract. ``R``, the
abundance correlation over shared clonotypes, is ``None`` because there are none to correlate.

For a whole cohort, load a metadata sheet once and let vdjtools do the pairing:

.. code-block:: python

   meta = vio.read_metadata("metadata.txt")               # sample_name + phenotype columns
   cohort = vio.read_samples(meta, base_dir="samples/")
   overlap.overlap_matrix(cohort)

Step 4 -- the same thing without Python
----------------------------------------

Every analysis is also a command. Commands take sample files or a metadata table, and write TSV or
Parquet:

.. code-block:: bash

   vdjtools diversity     sampleA.tsv sampleB.tsv -o diversity.tsv
   vdjtools segment-usage *.tsv --segment v -o usage.tsv
   vdjtools overlap       *.tsv -o overlap.tsv

   # a cohort from a metadata sheet, in parallel over samples
   vdjtools diversity -m metadata.txt --base-dir samples/ --threads 8 -o diversity.tsv

Omit ``-o`` and the result goes to stdout, so commands pipe. :doc:`cli` documents every command.

Step 5 -- generation probability
---------------------------------

How likely was a receptor to be produced by V(D)J recombination in the first place? A low
probability means a sequence you are unlikely to see twice by chance, which is what makes
convergent recombination interesting.

.. code-block:: python

   from vdjtools.model import load_bundled, native

   model = load_bundled("TRB")                    # isalgo working default
   native.pgen_aa(model, "CASSLAPGATNEKLFF")       # amino-acid generation probability
   native.pgen_aa(model, "CASSLAPGATNEKLFF", mismatches=1)   # and its whole Hamming-1 ball

Bundled model loading requires no download or OLGA installation. The replacement
isalgo catalog targets seven human and five mouse loci; consult :doc:`model` for
asset availability and the release gate. Native Pgen is checked against OLGA on identical probability tables and raw
germlines. The :ref:`model guide <olga-caveats>` explains germline repairs and oracle
workspace differences. vdjtools also supports tandem-D recombination. You can fit a model to your own
non-productive reads -- see :doc:`model`.

Step 6 -- one vector per repertoire
------------------------------------

For unrotated diversity, isotype and usage features with optional PCA, use
:doc:`features`. For a fixed, published reference system, use a signature.

When the question is a *model* rather than a statistic -- classify these donors, cluster these
samples -- you need each repertoire as a fixed-width row of numbers, on a scale a classifier can
consume. That is a :doc:`signature`:

.. code-block:: bash

   vdjtools signature --corpus blood samples/*.tsv.gz -o vsig.tsv

The ``--corpus`` argument is required and is the whole point: features are standardised against a
published reference corpus of real repertoires, so your matrix and a collaborator's are in the same
coordinate system without either of you fitting a scaler. Nine corpora are published; ``blood`` and
``tissue`` are fitted on 11,117 and 21,131 real samples respectively. The artifact is downloaded
and cached on first use.

Add the geometry half from `mirpy <https://github.com/antigenomics/mirpy>`_ and join on
``sample_id`` for the full vector:

.. code-block:: bash

   pip install mirpy-lib
   mir signature --corpus blood samples/*.tsv.gz -o rsig.tsv

Where to go next
----------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Next
     - Page
   * - Clean and normalise real data before analysing it
     - :doc:`preprocessing`
   * - Every analysis module, with runnable examples
     - :doc:`usage`
   * - Build, fit, check and compare recombination models
     - :doc:`model`
   * - Signatures: which corpus, and what the columns mean
     - :doc:`signature` | :doc:`channels`
   * - 10x and single-cell data
     - :doc:`singlecell`
   * - A term you did not recognise here
     - :doc:`glossary`
