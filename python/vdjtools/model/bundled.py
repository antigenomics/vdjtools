"""Precomputed recombination models shipped with vdjtools.

Two model sets, seven human loci each, live under ``vdjtools/model/_bundled/``:

- ``olga`` — imported from OLGA's default models (the exact-Pgen bootstrap; single-D).
- ``learned`` — EM-inferred from real out-of-frame reads (HuggingFace), tandem-D on the D-bearing
  loci (IGH/TRD/TRB). These carry a learned ``P(n_D=2)`` and broader trim/insertion distributions
  than the synthetic OLGA models.

Each model is a directory of parquet marginal tables + ``manifest.json`` (see :mod:`vdjtools.model`).
Provenance and the build command are recorded in ``SOURCES.md``.
"""
from __future__ import annotations

from importlib.resources import as_file, files

from .io import load_model
from .model import Model

#: The two shipped model sets.
SOURCES = ("olga", "learned")
#: The seven human loci with bundled models.
LOCI = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")


def load_bundled(locus: str, source: str = "olga") -> Model:
    """Load a precomputed model shipped with the package (no OLGA/HuggingFace at runtime).

    Args:
        locus: One of ``TRA TRB TRG TRD IGH IGK IGL`` (case-insensitive).
        source: ``"olga"`` (OLGA bootstrap, exact Pgen) or ``"learned"`` (EM-inferred from real
            out-of-frame reads; tandem-D on IGH/TRD/TRB).

    Returns:
        The :class:`~vdjtools.model.model.Model`.

    Raises:
        ValueError: If ``source`` is not one of :data:`SOURCES`.
        FileNotFoundError: If no bundled model exists for ``(source, locus)``.
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {SOURCES}, got {source!r}")
    root = files("vdjtools.model") / "_bundled" / source / locus.upper()
    if not root.is_dir():
        raise FileNotFoundError(f"no bundled {source!r} model for locus {locus!r}")
    with as_file(root) as path:
        return load_model(path)


def list_bundled() -> dict[str, list[str]]:
    """Return the available bundled models as ``{source: [locus, ...]}``."""
    out: dict[str, list[str]] = {}
    for src in SOURCES:
        base = files("vdjtools.model") / "_bundled" / src
        out[src] = sorted(d.name for d in base.iterdir() if d.is_dir()) if base.is_dir() else []
    return out
