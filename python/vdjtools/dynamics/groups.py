"""Metaclonotype-grouped dynamics — cluster near-variant clonotypes, then test the group.

A vaccine- or antigen-driven response is often *convergent*: many T-cell clones with slightly
different CDR3s recognise the same epitope, so the expansion signal is spread thin across a
family of near-variants and a per-exact-clonotype test (:func:`vdjtools.dynamics.test_pair`) sees
each one only weakly. Collapsing a CDR3-neighbourhood into **one metaclonotype** first, then
running the same paired test on the group, concentrates that signal.

The grouping is delegated to :func:`vdjtools.biomarker.metaclonotypes` (native ``seqtree`` /
``vdjmatch`` fuzzy search, single-linkage components) — ``scope="1,0,0,1"`` is a 1-substitution
(Hamming) ball, ``"1,1,1,1"`` a 1-edit (Levenshtein) ball. This is the VDJtrack 1-mismatch
annotation-expansion (``example.Rmd`` / ``vaccination.Rmd``) generalised into the test itself.
"""
from __future__ import annotations

import polars as pl

from ..io.schema import J_CALL, JUNCTION_AA, V_CALL
from .paired import test_pair


def test_metaclonotypes(a: pl.DataFrame, b: pl.DataFrame, *, scope: str = "1,0,0,1",
                        match_v: bool = True, match_j: bool = True, threads: int = 0,
                        **test_kw) -> pl.DataFrame:
    """Group CDR3-neighbour clonotypes into metaclonotypes, then :func:`test_pair` the groups.

    The union of both samples' clonotypes is clustered once into metaclonotypes; each sample's
    clonotypes are relabelled by ``meta_id``; and the paired within-donor test runs on the
    per-``meta_id`` summed counts — so a convergent expansion counts as a single feature.

    Args:
        a: The earlier sample (canonical clonotype frame).
        b: The later sample.
        scope: vdjmatch edit scope ``"subs,ins,dels,total"`` — ``"1,0,0,1"`` (default) is a
            1-substitution (Hamming) ball, ``"1,1,1,1"`` a 1-edit (Levenshtein) ball,
            ``"0,0,0,0"`` reduces to exact grouping.
        match_v: Require the same ``v_call`` to group two clonotypes.
        match_j: Require the same ``j_call`` to group.
        threads: Worker threads for the native search (``0`` = all cores).
        **test_kw: Forwarded to :func:`test_pair` (e.g. ``neff``, ``min_total``, ``alpha``).

    Returns:
        :func:`test_pair`'s output keyed by ``meta_id``, plus a representative ``junction_aa``
        (and ``v_call``/``j_call`` when matched) and ``n_variants`` per group, sorted by
        ``q_value``.

    Raises:
        ImportError: If vdjmatch is not importable (it is a base dependency).
    """
    from ..biomarker.metaclonotype import metaclonotypes

    grp_cols = [JUNCTION_AA] + ([V_CALL] if match_v else []) + ([J_CALL] if match_j else [])
    pooled = pl.concat(
        [a.select(grp_cols), b.select(grp_cols)], how="vertical_relaxed").unique()
    meta = metaclonotypes(pooled, scope=scope, match_v=match_v, match_j=match_j,
                          threads=threads)
    join_on = [c for c in grp_cols if c in meta.columns]

    a2 = a.join(meta, on=join_on, how="inner")
    b2 = b.join(meta, on=join_on, how="inner")
    res = test_pair(a2, b2, key=("meta_id",), **test_kw)

    reps = ([pl.col(JUNCTION_AA).first().alias("junction_aa")]
            + ([pl.col(V_CALL).first().alias("v_call")] if V_CALL in meta.columns else [])
            + ([pl.col(J_CALL).first().alias("j_call")] if J_CALL in meta.columns else []))
    rep = meta.group_by("meta_id").agg(*reps, pl.len().alias("n_variants"))
    return res.join(rep, on="meta_id", how="left").sort("q_value")


def _demo() -> None:
    """Self-check: 1-mismatch CDR3 variants collapse to one metaclonotype and test as a group."""
    import numpy as np

    from ..io.schema import COUNT

    try:
        from ..biomarker.metaclonotype import _require_vdjmatch
        _require_vdjmatch()
    except ImportError:
        print("groups._demo skipped: vdjmatch not importable")
        return

    rng = np.random.default_rng(0)
    AA = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    # A metaclonotype family: CASSLxPGATNEKLFF, one substitution apart at position x.
    fam = [f"CASSL{c}PGATNEKLFF" for c in "AGSTNDEQ"]
    # Random valid-amino-acid background clonotypes (>1 edit from the family and each other).
    bg = ["CASR" + "".join(rng.choice(AA, 6)) + "QYF" for _ in range(300)]
    v = "TRBV7-2"
    j = "TRBJ1-1"

    def frame(cdrs, counts):
        return pl.DataFrame({JUNCTION_AA: cdrs, V_CALL: [v] * len(cdrs),
                             J_CALL: [j] * len(cdrs), COUNT: counts})

    a = frame(fam + bg, [3] * len(fam) + list(rng.integers(50, 200, len(bg))))
    b = frame(fam + bg, [300] * len(fam) + list(rng.integers(50, 200, len(bg))))
    res = test_metaclonotypes(a, b, neff=None)          # counts are already the draw
    fam_row = res.filter(pl.col("n_variants") >= len(fam))
    assert fam_row.height == 1, res
    assert fam_row["dynamics"][0] == "expanded", fam_row
    print("groups._demo OK:\n", fam_row.select("meta_id", "n_variants", "count_a",
                                               "count_b", "dynamics", "q_value"))


if __name__ == "__main__":
    _demo()
