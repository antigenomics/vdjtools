"""Offline bundled models: fresh ``isalgo`` fits and frozen ``olga-legacy`` assets.

Compatibility ``learned`` and ``arda`` select fresh isalgo fits; ``olga`` selects
frozen OLGA assets. Old parameter families require a pinned historical release.
Missing new fits fail explicitly; they are never synthesized or substituted.
"""
from __future__ import annotations

from importlib.resources import as_file, files

from .io import load_model
from .collapse import collapse_alleles
from .model import Model

#: Canonical names and explicitly supported compatibility selectors.
SOURCES = ("isalgo", "olga-legacy", "olga", "learned", "arda")
_CANONICAL_SOURCES = ("isalgo", "olga-legacy")
#: The seven human receptor loci; isalgo mouse bundles cover five loci.
LOCI = ("TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL")
def _source_family(source: str) -> str:
    """Resolve supported historical spellings without selecting historical learned weights."""
    return {"learned": "isalgo", "arda": "isalgo", "olga": "olga-legacy"}.get(source, source)


def _bundled_key(source: str, locus: str, organism: str) -> str:
    """Directory name for one bundled model, per that set's naming scheme."""
    source = _source_family(source)
    if source == "isalgo" and organism.lower() == "mouse":
        return f"mouse_{locus.upper()}"
    if organism.lower() != "human":
        raise ValueError(
            f"vdjtools bundles the {source!r} model set for human only (got "
            f"organism={organism!r}); the isalgo set carries mouse TRA/TRB/IGH/IGK/IGL. "
            f"⚠ This is a gap in what is VENDORED, not in the upstream set: OLGA ships "
            f"mouse_T_alpha, mouse_T_beta, mouse_B_heavy, mouse_B_kappa and mouse_B_lambda, and "
            f"all five import cleanly through `from_olga`."
        )
    return locus.upper()


def load_bundled(locus: str, source: str = "isalgo", *, organism: str = "human",
                 collapse: bool = True) -> Model:
    """Load a precomputed model shipped with the package (no OLGA/HuggingFace at runtime).

    Args:
        locus: One of ``TRA TRB TRG TRD IGH IGK IGL`` (case-insensitive).
        source: ``"isalgo"`` (default; fresh ARDA-reference EM fits), ``"olga-legacy"``
            (frozen human bootstrap), or compatibility ``"olga"``, ``"learned"``, ``"arda"``.
        organism: ``"human"`` (default) or ``"mouse"``. The isalgo mouse catalog covers TRA/TRB/IGH/IGK/IGL.
            No mouse TRG/TRD model is bundled. Their data acquisition is a separate
            future project. Missing combinations fail rather than using a human model.
        collapse: If ``True`` (default), collapse multi-allele models to one ``*01`` entry per
            gene via :func:`~vdjtools.model.collapse.collapse_alleles`. Already gene-level isalgo
            fits retain their selected reference labels and tables. Pass ``collapse=False`` for
            full allele resolution, e.g. the exact-OLGA-Pgen fidelity check.

    Returns:
        The :class:`~vdjtools.model.model.Model`.

    Raises:
        ValueError: If ``source`` is not one of :data:`SOURCES`, or a non-human ``organism`` is
            asked of a set vendored here for human only.
        FileNotFoundError: If no bundled model exists for ``(source, organism, locus)``.

    Example:
        >>> m = load_bundled("TRB")                                  # isalgo, human
        >>> m_arda = load_bundled("TRB", "arda")                     # alias of isalgo, human
        >>> m_mouse = load_bundled("TRB", "arda", organism="mouse")  # alias of isalgo, mouse
    """
    if source not in SOURCES:
        raise ValueError(f"source must be one of {SOURCES}, got {source!r}")
    source = _source_family(source)
    key = _bundled_key(source, locus, organism)
    directory = "olga" if source == "olga-legacy" else source
    root = files("vdjtools.model") / "_bundled" / directory / key
    if not root.is_dir():
        raise FileNotFoundError(
            f"no bundled {source!r} model for locus {locus!r} (organism {organism!r}); "
            f"available: {list_bundled().get(source, [])}"
        )
    with as_file(root) as path:
        m = load_model(path)
    if collapse and source == "isalgo" and all(g["gene"].n_unique() == g.height for g in m.genomic.values()):
        return m  # Already fitted at gene resolution; keep the selected reference allele labels.
    return collapse_alleles(m) if collapse else m


def list_bundled() -> dict[str, list[str]]:
    """Return canonical bundled families as ``{source: [key, ...]}``.

    Compatibility names resolve to these families; old weights require a pinned prior release.

    Keys are bare human loci, plus ``mouse_{LOCUS}`` for isalgo mouse fits.
    Empty canonical families remain visible so missing release assets cannot be hidden.
    """
    out: dict[str, list[str]] = {}
    for src in _CANONICAL_SOURCES:
        directory = "olga" if src == "olga-legacy" else src
        base = files("vdjtools.model") / "_bundled" / directory
        out[src] = sorted(d.name for d in base.iterdir() if d.is_dir()) if base.is_dir() else []
    return out
