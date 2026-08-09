"""Precomputed recombination models shipped with vdjtools.

Three model sets live under ``vdjtools/model/_bundled/``:

- ``olga`` — imported from OLGA's default models (the exact-Pgen bootstrap; single-D).
  Seven human loci, OLGA's germline namespace.
- ``learned`` — EM-inferred from real non-functional reads (out-of-frame + stop-codon, HuggingFace), tandem-D on the D-bearing
  loci (IGH/TRD/TRB). These carry a learned ``P(n_D=2)`` and broader trim/insertion distributions
  than the synthetic OLGA models. Seven human loci, OLGA's germline namespace.
- ``arda`` — the same EM inference on the **arda** IMGT allele namespace (:func:`from_arda`
  scaffold refit by :func:`~vdjtools.model.infer.infer_native`). Nine models: the seven human
  loci **plus mouse TRA/TRB — the only bundled set with a non-human organism**. Use this one when
  the rest of your pipeline is arda-annotated, so generated sequences share one allele namespace
  with your query data instead of needing a name fallback.

``arda`` models are keyed by ``{organism}_{LOCUS}`` on disk; the other two by ``{LOCUS}`` alone.
:func:`load_bundled` hides that — pass ``organism=`` and it picks the right key.

Each model is a directory of parquet marginal tables + ``manifest.json`` (see :mod:`vdjtools.model`).
Provenance and the build command are recorded in ``SOURCES.md``.
"""
from __future__ import annotations

from importlib.resources import as_file, files

from .io import load_model
from .collapse import collapse_alleles
from .model import Model

#: The three shipped model sets.
SOURCES = ("olga", "learned", "arda")
#: The seven human loci with bundled models (``arda`` adds mouse TRA/TRB).
LOCI = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")
#: Sets whose directories are keyed ``{organism}_{LOCUS}`` rather than ``{LOCUS}``.
_ORGANISM_KEYED = ("arda",)


def _bundled_key(source: str, locus: str, organism: str) -> str:
    """Directory name for one bundled model, per that set's naming scheme."""
    if source in _ORGANISM_KEYED:
        return f"{organism.lower()}_{locus.upper()}"
    if organism.lower() != "human":
        raise ValueError(
            f"the {source!r} model set is human-only (got organism={organism!r}); "
            f"use source='arda' for mouse"
        )
    return locus.upper()


def load_bundled(locus: str, source: str = "olga", *, organism: str = "human",
                 collapse: bool = True) -> Model:
    """Load a precomputed model shipped with the package (no OLGA/HuggingFace at runtime).

    Args:
        locus: One of ``TRA TRB TRG TRD IGH IGK IGL`` (case-insensitive).
        source: ``"olga"`` (OLGA bootstrap, exact Pgen), ``"learned"`` (EM-inferred from real
            non-functional reads; tandem-D on IGH/TRD/TRB), or ``"arda"`` (the same EM inference on
            the **arda IMGT allele namespace**, and the only set covering mouse).
        organism: ``"human"`` (default) or ``"mouse"``. Only the ``arda`` set ships a non-human
            organism; asking the others for one is an error rather than a silent human model.
        collapse: If ``True`` (default), collapse each gene to a single ``*01`` allele via
            :func:`~vdjtools.model.collapse.collapse_alleles` — the working gene-level resolution,
            in which Pgen also collapses a clonotype's allele to ``*01`` (short-read aligners cannot
            resolve alleles reliably, so the suffix is noise). Pass ``collapse=False`` for full
            allele resolution, e.g. the exact-OLGA-Pgen fidelity check.

    Returns:
        The :class:`~vdjtools.model.model.Model`.

    Raises:
        ValueError: If ``source`` is not one of :data:`SOURCES`, or a non-human ``organism`` is
            asked of a human-only set.
        FileNotFoundError: If no bundled model exists for ``(source, organism, locus)``.

    Example:
        >>> m = load_bundled("TRB")                                  # olga, human
        >>> m_arda = load_bundled("TRB", "arda")                     # arda namespace, human
        >>> m_mouse = load_bundled("TRB", "arda", organism="mouse")  # the only mouse models
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {SOURCES}, got {source!r}")
    key = _bundled_key(source, locus, organism)
    root = files("vdjtools.model") / "_bundled" / source / key
    if not root.is_dir():
        raise FileNotFoundError(
            f"no bundled {source!r} model for locus {locus!r} (organism {organism!r}); "
            f"available: {list_bundled().get(source, [])}"
        )
    with as_file(root) as path:
        m = load_model(path)
    return collapse_alleles(m) if collapse else m


def list_bundled() -> dict[str, list[str]]:
    """Return the available bundled models as ``{source: [key, ...]}``.

    Keys are bare loci for ``olga``/``learned`` and ``{organism}_{LOCUS}`` for ``arda``, matching
    what :func:`load_bundled` accepts for each set.
    """
    out: dict[str, list[str]] = {}
    for src in SOURCES:
        base = files("vdjtools.model") / "_bundled" / src
        out[src] = sorted(d.name for d in base.iterdir() if d.is_dir()) if base.is_dir() else []
    return out
