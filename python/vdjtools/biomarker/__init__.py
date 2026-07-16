"""vdjtools.biomarker — incidence association, co-occurrence, and metaclonotype grouping.

Clonotype-association testing across a cohort of repertoires (Emerson 2017, Howie 2015,
De Witt 2018, Vlasova 2026). :func:`association` tests each clonotype feature's subject
incidence against a condition (binary / category / stratified) with a choice of test;
:func:`cooccurrence` tests feature-vs-feature co-occurrence (α-β pairing, same-chain);
:mod:`condition` builds the phenotype design; :func:`fisher_association` is the Emerson
Fisher special case; :func:`metaclonotypes` groups 1-mismatch variants.
"""
from . import condition, stats
from .association import association, select_candidates
from .cooccurrence import cooccurrence
from .fisher import fisher_association
from .metaclonotype import metaclonotypes

__all__ = [
    "association",          # general incidence association (any test / condition / scope)
    "cooccurrence",         # feature-vs-feature co-occurrence (α-β + same-chain)
    "select_candidates",    # public features over an incidence count / fraction threshold
    "fisher_association",   # Emerson-2017 Fisher special case (back-compat)
    "metaclonotypes",       # edit-scope metaclonotype grouping (1mm matching)
    "condition",            # phenotype-design builders (binary / categorical / HLA / CMH)
    "stats",                # vectorised 2×2 test kernels
]
