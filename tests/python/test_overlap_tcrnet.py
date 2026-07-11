"""Tests for vdjtools.overlap.tcrnet — neighbourhood-enrichment via vdjmatch e-values.

A synthetic planted clique (convergent cluster) must light up as enriched against a
sea of random CDR3s; a real antigen-specific TRB set (VDJdb GILGFVFTL) must show more
enrichment than a size-matched random OLGA-generated set. vdjmatch (the ``overlap``
extra) is guarded; the real-data files skip cleanly when absent.
"""
import gzip
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from vdjtools.io import schema as S
from vdjtools import overlap as O

pytest.importorskip("vdjmatch")

_MIRPY_ASSETS = Path("/Users/mikesh/vcs/code/mirpy/tests/assets")
_AA = "ACDEFGHIKLMNPQRSTVWY"


def _sample_from_cdr3(cdr3, v="TRBV7-9", j="TRBJ2-1"):
    n = len(cdr3)
    df = pl.DataFrame({S.V_CALL: [v] * n, S.J_CALL: [j] * n,
                       S.CDR3_AA: cdr3, S.COUNT: [1] * n})
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
    assert set(res.columns) >= {S.CDR3_AA, "n_neighbors", "n_control", "E",
                                "p_enrichment"}
    clu = res.filter(pl.col(S.CDR3_AA).is_in(clique))
    noi = res.filter(~pl.col(S.CDR3_AA).is_in(clique))

    # Every clique member sees the other 7; the random background sees essentially none.
    assert clu["n_neighbors"].min() >= 7
    assert (clu["p_enrichment"] < 0.001).all()
    assert clu["p_enrichment"].max() < noi["p_enrichment"].min()


def test_tcrnet_explicit_locus_and_no_exact():
    """Locus can be passed explicitly; exclude_exact toggles self-hit removal."""
    sample = _sample_from_cdr3(["CASSLAPGELFF", "CASSLAPGELFY", "CASSQQTGELFF"])
    res = O.tcrnet(sample, locus="TRB", scope="1,0,0,1", exclude_exact=True)
    assert res.height == 3


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
