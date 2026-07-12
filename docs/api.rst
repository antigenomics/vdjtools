API reference
=============

.. note::

   The API is filling in phase by phase. Module pages are added here as each subpackage
   lands (``vdjtools.model``, ``vdjtools.stats``, ``vdjtools.features``, …).

vdjtools
--------

.. automodule:: vdjtools
   :members:
   :undoc-members:
   :show-inheritance:

Model engine (``vdjtools.model``)
---------------------------------

The native V(D)J recombination engine: a model is a directory of tidy ``polars`` marginal
tables plus a ``manifest.json`` declaring the recombination Bayes net. It supersedes OLGA
(generation probability, sampling) and IGoR (EM inference), adds tandem-D (D-D) support, and
exposes information-theoretic diagnostics.

Model container, schema and events
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.model
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.schema
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.events
   :members:
   :undoc-members:
   :show-inheritance:

Import, germline reference and stitching
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.io
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.bundled
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.reference
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.stitch
   :members:
   :undoc-members:
   :show-inheritance:

Generation probability, sampling and inference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.pgen
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.native
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.generate
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: vdjtools.model.infer
   :members:
   :undoc-members:
   :show-inheritance:

Tandem-D (D-D) extension
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.dd
   :members:
   :undoc-members:
   :show-inheritance:

Model diagnostics
~~~~~~~~~~~~~~~~~~

.. automodule:: vdjtools.model.analyze
   :members:
   :undoc-members:
   :show-inheritance:
