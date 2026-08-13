"""The ``vsig:kmer`` block: a V+k-mer profile projected onto a frozen basis.

The feature space itself lives in :mod:`vdjtools.features.kmer_space`; this is the thin layer
that gives it column names, a tier, and a place in the layout.

It is **not** in :data:`vdjtools.signature.layout._BLOCKS`, and that is deliberate. The block's
width is the SVD rank that was actually fitted, and its columns mean nothing without the
vocabulary and IDF they were fitted with -- so unlike every other block it cannot be declared
before its artifact exists. Declaring it early would put N permanently-nan columns into the
contract, which anyone downstream would read as "my samples were too shallow", not as "this was
never computed". :func:`register_kmer` is how it arrives.
"""
from __future__ import annotations

import numpy as np

from ..features.kmer_space import KmerSpace
from . import layout as L


def kmer_block(df, locus: str, space: KmerSpace | None, *, weight: str = "freq") -> dict:
    """Project one repertoire onto its locus's frozen k-mer basis.

    Args:
        df: The locus's work frame -- already filtered, and carrying whatever ``frequency`` the
            caller's clone-weight ladder wrote.
        locus: Locus name, for the error message only.
        space: The fitted space for this locus, or ``None`` if none was fitted.
        weight: Clone weight ladder.

    Returns:
        ``{"PC01": ..., ...}``. All-``nan`` when no space is available or the frame is empty --
        a hole, which the mask reports, rather than a zero, which is a coordinate.
    """
    if space is None or space.n_components == 0:
        n = 0 if space is None else space.n_columns
        return {f"PC{i + 1:02d}": np.nan for i in range(n)}
    names = [f"PC{i + 1:02d}" for i in range(space.n_components)]
    if df is None or df.height == 0:
        return dict.fromkeys(names, np.nan)
    return dict(zip(names, (float(v) for v in space.transform(df, weight=weight))))


def kmer_spec(spaces: dict[str, KmerSpace], *, tier: str = "full") -> L.Block:
    """The :class:`~vdjtools.signature.layout.Block` describing a set of fitted spaces.

    Every locus must have been fitted to the same rank: a block whose width varies by locus is
    not one block, and the column contract is a fixed-width index subset by construction.

    ``attributable=True``: unlike a Hill number, "which clonotypes drive this column" is a
    well-posed question here -- a k-mer column has a literal clonotype pre-image, the clones
    whose junction contains it under that V gene.
    """
    if not spaces:
        raise ValueError("no fitted k-mer spaces")
    ranks = {loc: sp.n_components for loc, sp in spaces.items()}
    if len(set(ranks.values())) != 1:
        raise ValueError(f"every locus must be fitted to the same rank; got {ranks}")
    d = next(iter(ranks.values()))
    if d < 1:
        raise ValueError("fitted spaces carry no components; fit with n_components > 0")
    return L.Block("vsig", "kmer",
                   L.feats(tier, "none", *(f"PC{i + 1:02d}" for i in range(d))),
                   loci=tuple(sorted(spaces)), attributable=True)


def register_kmer(spaces: dict[str, KmerSpace], *, tier: str = "full") -> L.Block:
    """Build the block from fitted spaces and add it to the layout registry.

    The transform is ``"none"``: the projection is already onto an orthonormal basis of
    L2-normalised TF-IDF rows, so the coordinates are on a common scale by construction and a
    variance-stabilising transform on top would only distort it. The frozen per-column reference
    still rescales them like anything else.
    """
    block = kmer_spec(spaces, tier=tier)
    L.register(block)
    return block
