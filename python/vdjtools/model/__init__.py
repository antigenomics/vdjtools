"""vdjtools.model — native V(D)J recombination model: Pgen, generation, EM inference.

A model is a directory of tidy polars marginal tables + a ``manifest.json`` declaring the
recombination Bayes net (:mod:`~vdjtools.model.events`). Bootstrap models are imported from
OLGA's format with :func:`from_olga`; native models round-trip through :func:`save_model` /
:func:`load_model` with no OLGA dependency.
"""
from .events import Event, EventKind
from .io import from_olga, load_model, save_model
from .model import Model
from .schema import Manifest

__all__ = ["Event", "EventKind", "Manifest", "Model", "from_olga", "load_model", "save_model"]
