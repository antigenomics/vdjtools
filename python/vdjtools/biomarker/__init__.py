"""vdjtools.biomarker — Incidence-based association (Fisher) and metaclonotype grouping."""
from .fisher import fisher_association
from .metaclonotype import metaclonotypes

__all__ = [
    "fisher_association",   # Fisher-exact clonotype-incidence vs phenotype association
    "metaclonotypes",       # edit-scope metaclonotype grouping (1mm matching)
]
