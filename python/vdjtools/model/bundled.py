"""Precomputed models loaded without network access or an OLGA runtime dependency.

``learned`` is the working default: precomputed parameters for twelve contexts,
seven human loci and mouse TRA/TRB/IGH/IGK/IGL. Each manifest describes its
parameter source and reference. Human keys are bare loci;
mouse keys are ``mouse_{LOCUS}``. ``arda`` retains its legacy organism-prefixed
paths; its corrected fitted assets remain available explicitly. ``olga`` retains
the seven human bootstrap models for exact oracle and historical reproduction.
The installed asset directories define availability; missing fits are never synthesized.
"""
from __future__ import annotations

from importlib.resources import as_file, files

from .io import load_model
from .collapse import collapse_alleles
from .model import Model

#: The three shipped model sets.
SOURCES = ("olga", "learned", "arda")
#: The seven human receptor loci; learned mouse bundles cover five loci.
LOCI = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")
#: Sets whose directories are keyed ``{organism}_{LOCUS}`` rather than ``{LOCUS}``.
_ORGANISM_KEYED = ("arda",)


def _bundled_key(source: str, locus: str, organism: str) -> str:
    """Directory name for one bundled model, per that set's naming scheme."""
    if source in _ORGANISM_KEYED:
        return f"{organism.lower()}_{locus.upper()}"
    if source == "learned" and organism.lower() == "mouse":
        return f"mouse_{locus.upper()}"
    if organism.lower() != "human":
        raise ValueError(
            f"vdjtools bundles the {source!r} model set for human only (got "
            f"organism={organism!r}); the learned set carries mouse TRA/TRB/IGH/IGK/IGL. "
            f"⚠ This is a gap in what is VENDORED, not in the upstream set: OLGA ships "
            f"mouse_T_alpha, mouse_T_beta, mouse_B_heavy, mouse_B_kappa and mouse_B_lambda, and "
            f"all five import cleanly through `from_olga`."
        )
    return locus.upper()


def load_bundled(locus: str, source: str = "learned", *, organism: str = "human",
                 collapse: bool = True) -> Model:
    """Load a precomputed model shipped with the package (no OLGA/HuggingFace at runtime).

    Args:
        locus: One of ``TRA TRB TRG TRD IGH IGK IGL`` (case-insensitive).
        source: ``"learned"`` (default; EM fits on ARDA references), ``"olga"``
            (human bootstrap for exact Pgen reproduction), or legacy ``"arda"`` paths.
        organism: ``"human"`` (default) or ``"mouse"``. The learned mouse catalog covers TRA/TRB/IGH/IGK/IGL.
            No mouse TRG/TRD model is bundled. Their data acquisition is a separate
            future project. Missing combinations fail rather than using a human model.
        collapse: If ``True`` (default), collapse multi-allele models to one ``*01`` entry per
            gene via :func:`~vdjtools.model.collapse.collapse_alleles`. Already gene-level learned
            fits retain their selected reference labels and tables. Pass ``collapse=False`` for
            full allele resolution, e.g. the exact-OLGA-Pgen fidelity check.

    Returns:
        The :class:`~vdjtools.model.model.Model`.

    Raises:
        ValueError: If ``source`` is not one of :data:`SOURCES`, or a non-human ``organism`` is
            asked of a set vendored here for human only.
        FileNotFoundError: If no bundled model exists for ``(source, organism, locus)``.

    Example:
        >>> m = load_bundled("TRB")                                  # learned, human
        >>> m_arda = load_bundled("TRB", "arda")                     # legacy path, human
        >>> m_mouse = load_bundled("TRB", "arda", organism="mouse")  # legacy mouse path
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
    if source == "learned" and all(g["gene"].n_unique() == g.height for g in m.genomic.values()):
        return m  # Already fitted at gene resolution; keep the selected reference allele labels.
    return collapse_alleles(m) if collapse else m


def list_bundled() -> dict[str, list[str]]:
    """Return the available bundled models as ``{source: [key, ...]}``.

    Keys are bare human loci for ``olga``/``learned``, ``mouse_{LOCUS}`` for learned mouse,
    and ``{organism}_{LOCUS}`` for ``arda``, matching
    what :func:`load_bundled` accepts for each set.
    """
    out: dict[str, list[str]] = {}
    for src in SOURCES:
        base = files("vdjtools.model") / "_bundled" / src
        out[src] = sorted(d.name for d in base.iterdir() if d.is_dir()) if base.is_dir() else []
    return out
