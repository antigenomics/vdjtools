"""Tests for vdjtools.overlap.tcrnet — neighbourhood-enrichment via vdjmatch e-values.

A synthetic planted clique (convergent cluster) must light up as enriched against a
sea of random CDR3s; a real antigen-specific TRB set (VDJdb GILGFVFTL) must show more
enrichment than a size-matched random OLGA-generated set. vdjmatch (the ``overlap``
extra) is guarded; the real-data files skip cleanly when absent.
"""
import gzip
import os
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from scipy.stats import poisson

from vdjtools.io import schema as S
from vdjtools import overlap as O

pytest.importorskip("vdjmatch")

# Optional real-data assets (gitignored); the test below skips cleanly when absent. Drop them
# (or symlink) under ./data_dump/mirpy-assets/, or point $MIRPY_ASSETS at a checkout.
_MIRPY_ASSETS = Path(os.environ.get("MIRPY_ASSETS", "data_dump/mirpy-assets"))
_AA = "ACDEFGHIKLMNPQRSTVWY"


def _sample_from_cdr3(cdr3, v="TRBV7-9", j="TRBJ2-1"):
    n = len(cdr3)
    df = pl.DataFrame({S.V_CALL: [v] * n, S.J_CALL: [j] * n,
                       S.JUNCTION_AA: cdr3, S.COUNT: [1] * n})
    return S.add_locus(S.normalize(df, recompute_freq=True))


def test_tcrnet_planted_cluster_enriched():
    """A clique of 8 mutually-1-substitution CDR3s (single varying position) buried
    among random CDR3s must be flagged as neighbourhood-enriched; the random
    background must not be."""
    rng = np.random.default_rng(1)
    tmpl = list("CASSXPGELFFYT")
    pos = 4
    clique = ["".join(tmpl[:pos] + [r] + tmpl[pos + 1:]) for r in "ACDEFGHI"]
    noise = ["".join(rng.choice(list(_AA), 13)) for _ in range(250)]
    sample = _sample_from_cdr3(clique + noise)

    res = O.tcrnet(sample, scope="1,0,0,1")
    assert set(res.columns) >= {S.JUNCTION_AA, "n_neighbors", "n_control", "E",
                                "p_enrichment"}
    clu = res.filter(pl.col(S.JUNCTION_AA).is_in(clique))
    noi = res.filter(~pl.col(S.JUNCTION_AA).is_in(clique))

    # Every clique member sees the other 7; the random background sees essentially none.
    assert clu["n_neighbors"].min() >= 7
    assert (clu["p_enrichment"] < 0.001).all()
    assert clu["p_enrichment"].max() < noi["p_enrichment"].min()


def _neighbors(res):
    return dict(zip(res[S.JUNCTION_AA].to_list(), res["n_neighbors"].to_list()))


def test_tcrnet_explicit_locus_and_exclude_exact_toggle():
    """Locus can be passed explicitly; exclude_exact toggles self-hit removal, which
    is pinned by the neighbour counts (not just the row count)."""
    # CASSLAPGELFF <-1-> CASSLAPGELFY (1 sub); CASSQQTGELFF is 3 subs from both.
    sample = _sample_from_cdr3(["CASSLAPGELFF", "CASSLAPGELFY", "CASSQQTGELFF"])
    incl = O.tcrnet(sample, locus="TRB", scope="1,0,0,1", exclude_exact=True)
    excl = O.tcrnet(sample, locus="TRB", scope="1,0,0,1", exclude_exact=False)
    # exclude_exact drops the distance-0 self-hit; without it every clonotype gains +1.
    assert _neighbors(incl) == {"CASSLAPGELFF": 1, "CASSLAPGELFY": 1, "CASSQQTGELFF": 0}
    assert _neighbors(excl) == {"CASSLAPGELFF": 2, "CASSLAPGELFY": 2, "CASSQQTGELFF": 1}


def test_tcrnet_counts_public_clones_in_the_background():
    """The control-side self-hit fix, pinned.

    The query is a member of its own sample's index but NOT of the control, so the two sides
    need opposite treatment. vdjmatch's ``exclude_exact=True`` punctures distance-0 hits on
    *both* (its own docstring: "use when queries may be members of the target/control"), which
    silently discarded genuine public-clone neighbours from the background, understating E and
    inflating significance — anti-conservative, and worst exactly on the public clonotypes a
    TCRnet run is most likely to report.

    Here the control CONTAINS the query verbatim. That exact match is a real background
    neighbour: a background repertoire that carries this clonotype is evidence the clonotype is
    unremarkable, not evidence it is enriched.
    """
    import seqtree

    q = "CASSLAPGELFF"
    sample = _sample_from_cdr3([q, "CASSLAPGELFY", "CASSQQTGELFF"])
    # A background carrying the query itself plus one 1-mismatch variant of it.
    control = seqtree.Index.build([q, "CASSLAPGELFW", "CWWWWWWWWWWF"], alphabet="aa")
    res = O.tcrnet(sample, control=control, locus="TRB", scope="1,0,0,1")
    row = res.filter(pl.col(S.JUNCTION_AA) == q).row(0, named=True)

    # 2 background neighbours: the exact copy AND the 1-sub variant. Dropping the exact copy
    # (the bug) would leave 1 and halve E.
    assert row["n_control"] == 2
    assert row["E"] == pytest.approx(2 * len(sample) / 3)
    # ...and the p-value is computed from the self-corrected degree against that E.
    assert row["n_neighbors"] == 1                      # CASSLAPGELFY, self excluded
    assert row["p_enrichment"] == pytest.approx(poisson.sf(0, row["E"]))


def test_tcrnet_returns_a_q_value():
    """Raw Poisson tails over ~1e5 clonotypes are not a threshold anyone can apply; the
    multiple-testing burden spans every clonotype scored in one call."""
    sample = _sample_from_cdr3(["CASSLAPGELFF", "CASSLAPGELFY", "CASSQQTGELFF"])
    res = O.tcrnet(sample, locus="TRB", scope="1,0,0,1")
    assert "q_value" in res.columns
    assert (res["q_value"].to_numpy() >= res["p_enrichment"].to_numpy() - 1e-12).all()
    assert ((res["q_value"] >= 0) & (res["q_value"] <= 1)).all()


def test_tcrnet_mixed_locus_scored_per_locus():
    """A sample mixing TRA and TRB clonotypes must be scored PER LOCUS — each locus
    against its own background — not all against a single inferred control."""
    trb = _sample_from_cdr3(
        ["".join(list("CASSXPGELFF")[:4] + [r] + list("CASSXPGELFF")[5:]) for r in "ACDEF"],
        v="TRBV7-9", j="TRBJ2-1")
    tra = _sample_from_cdr3(
        ["".join(list("CAVRXDDKIIF")[:4] + [r] + list("CAVRXDDKIIF")[5:]) for r in "ACDEF"],
        v="TRAV1-2", j="TRAJ33")
    res = O.tcrnet(pl.concat([trb, tra]), scope="1,0,0,1")

    assert set(res[S.LOCUS].unique().to_list()) == {"TRA", "TRB"}
    tra_rows = res.filter(pl.col(S.LOCUS) == "TRA")
    trb_rows = res.filter(pl.col(S.LOCUS) == "TRB")
    # Both planted cliques enriched; and TRA is scored against the TRA background
    # (real n_control neighbours) — under the old single-TRB-background bug a TRA
    # sequence would have ~0 control neighbours and be spuriously (mis)calibrated.
    assert (tra_rows["p_enrichment"] < 0.05).all()
    assert (trb_rows["p_enrichment"] < 0.05).all()
    assert tra_rows["n_control"].max() > 0


def _read_cdr3_gz(path: Path, col: int = 0):
    with gzip.open(path, "rt") as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    seqs = [r[col] if len(r) > col else r[0] for r in rows]
    return [s for s in dict.fromkeys(seqs) if s and s.isalpha()]


def _frac_significant(cdr3):
    sample = _sample_from_cdr3(cdr3, v="TRBV", j="TRBJ")
    res = O.tcrnet(sample, locus="TRB", species="human", scope="1,0,0,1")
    return float((res["p_enrichment"] < 0.05).mean())


def test_tcrnet_gilgfvftl_vs_random_real_data():
    """VDJdb GILGFVFTL-specific TRB CDR3s (antigen-driven, convergent) must show more
    neighbourhood enrichment than a size-matched random OLGA-generated TRB set."""
    gil_path = _MIRPY_ASSETS / "gilgfvftl_trb_junctions.txt.gz"
    olga_path = _MIRPY_ASSETS / "olga_humanTRB_1000.txt.gz"
    if not (gil_path.exists() and olga_path.exists()):
        pytest.skip("mirpy antigen-specific / OLGA assets not available")

    gil = _read_cdr3_gz(gil_path, col=0)          # one junction per line
    olga = _read_cdr3_gz(olga_path, col=1)         # cols: nt, aa, v, j
    frac_gil = _frac_significant(gil)
    frac_olga = _frac_significant(olga)

    # Antigen-specific set is markedly more enriched than random (empirically ~0.29
    # vs ~0.10); assert a clear, version-robust margin.
    assert frac_gil > 0.20
    assert frac_gil > 1.5 * frac_olga
