"""vdjtools.model — native V(D)J recombination model: Pgen, generation, EM inference.

A model is a directory of tidy polars marginal tables + a ``manifest.json`` declaring the
recombination Bayes net (:mod:`~vdjtools.model.events`). Bootstrap models are imported from
OLGA's format with :func:`from_olga`; native models round-trip through :func:`save_model` /
:func:`load_model` with no OLGA dependency.
"""
from . import analyze, check, score
from .bundled import list_bundled, load_bundled
from .check import check_model
from .collapse import collapse_alleles
from .events import Event, EventKind
from .io import (
    from_arda,
    from_germline,
    from_olga,
    load_model,
    marginals_frame,
    save_model,
    set_marginals,
)
from .model import Model
from .reference import (
    cut_segment,
    load_germline,
    read_germline_fasta,
    reconcile_olga,
    reverse_complement,
    translate,
    validate_germline,
)
from .rescale import rescale_usage
from .schema import Manifest
from .stitch import stitch_contig, stitch_frame
from .viterbi import Scenario, best_scenario

__all__ = [
    "analyze",
    "check",
    "score",
    "check_model",
    "marginals_frame",
    "set_marginals",
    "Event",
    "EventKind",
    "Manifest",
    "Model",
    "from_olga",
    "from_arda",
    "from_germline",
    "read_germline_fasta",
    "validate_germline",
    "load_model",
    "save_model",
    "load_bundled",
    "list_bundled",
    "collapse_alleles",
    "rescale_usage",
    "load_germline",
    "cut_segment",
    "reconcile_olga",
    "reverse_complement",
    "translate",
    "stitch_contig",
    "stitch_frame",
    "Scenario",
    "best_scenario",
]
