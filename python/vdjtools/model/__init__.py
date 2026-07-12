"""vdjtools.model — native V(D)J recombination model: Pgen, generation, EM inference.

A model is a directory of tidy polars marginal tables + a ``manifest.json`` declaring the
recombination Bayes net (:mod:`~vdjtools.model.events`). Bootstrap models are imported from
OLGA's format with :func:`from_olga`; native models round-trip through :func:`save_model` /
:func:`load_model` with no OLGA dependency.
"""
from . import analyze
from .bundled import list_bundled, load_bundled
from .events import Event, EventKind
from .io import from_olga, load_model, save_model
from .model import Model
from .reference import cut_segment, load_germline, reconcile_olga, reverse_complement, translate
from .schema import Manifest
from .stitch import stitch_contig, stitch_frame

__all__ = [
    "analyze",
    "Event",
    "EventKind",
    "Manifest",
    "Model",
    "from_olga",
    "load_model",
    "save_model",
    "load_bundled",
    "list_bundled",
    "load_germline",
    "cut_segment",
    "reconcile_olga",
    "reverse_complement",
    "translate",
    "stitch_contig",
    "stitch_frame",
]
