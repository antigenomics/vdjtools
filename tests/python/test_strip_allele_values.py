"""``strip_allele_values`` must equal the expression it is built from, at every input shape.

The expression is the only statement of the *rule*; this resolves it once per distinct call
instead of once per row, because a repertoire has 10^5-10^6 rows and 10^1-10^2 distinct segment
calls. So the test that matters is equality with the expression -- if the two ever disagree the
faster one is simply wrong, and the ambiguity handling is exactly where that would show up first:
``strip_allele`` keeps a genuine cross-gene tie (``IGHV3-23,IGHV3-23D``) and collapses an
allele-level tie within one gene (``IGHV1-2*02,IGHV1-2*04`` -> ``IGHV1-2``).
"""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjtools.io.schema import strip_allele, strip_allele_values


def _expr(s: pl.Series) -> list:
    return pl.DataFrame({"g": s}).select(strip_allele(pl.col("g")).alias("g"))["g"].to_list()


#: Every shape that has ever mattered here: plain, ambiguous both ways, whitespace, no allele,
#: a bare suffix, a lone separator, empty string and null.
ADVERSARIAL = [
    None, "", "TRBV12-3*01", "TRBV12-3", "A*01,A*02", "A*01,B*01", "B*01,A*01",
    "IGHV3-23*01,IGHV3-23D*01", "IGHV1-2*02,IGHV1-2*04", " X*01 , Y*02 ", "NOSTAR",
    "*01", ",", "IGHV1-69*01,IGHV1-69D*01,IGHV1-69*04",
]


def test_it_equals_the_expression_on_every_adversarial_shape():
    s = pl.Series("v_call", ADVERSARIAL, dtype=pl.Utf8)
    assert strip_allele_values(s).to_list() == _expr(s)


@pytest.mark.parametrize("n", [1, 2, 1000, 20_000])
def test_it_equals_the_expression_at_scale(n):
    """Scale matters because the dedup path only engages when values repeat."""
    rng = np.random.default_rng(n)
    genes = [f"TRBV{i}-{j}" for i in range(1, 12) for j in range(1, 4)]
    s = pl.Series("v_call", [f"{g}{a}" for g, a in
                             zip(rng.choice(genes, n), rng.choice(["*01", "*02", ""], n))])
    assert strip_allele_values(s).to_list() == _expr(s)


def test_the_degenerate_columns_do_not_raise():
    """An empty or all-null column has nothing distinct to resolve, which is not an error."""
    assert strip_allele_values(pl.Series("v", [], dtype=pl.Utf8)).to_list() == []
    allnull = pl.Series("v", [None] * 7, dtype=pl.Utf8)
    assert strip_allele_values(allnull).to_list() == [None] * 7 == _expr(allnull)


def test_a_null_survives_beside_real_calls():
    """The dedup keys on non-null values, so a null has to be carried by the default, not dropped."""
    s = pl.Series("v", ["TRBV6-5*01", None, "TRBV6-5*02", None], dtype=pl.Utf8)
    assert strip_allele_values(s).to_list() == ["TRBV6-5", None, "TRBV6-5", None]


def test_it_keeps_the_name_and_dtype_of_its_input():
    """Callers alias it straight into a frame, so a renamed column would land in the wrong place."""
    out = strip_allele_values(pl.Series("j_call", ["TRBJ2-7*01"], dtype=pl.Utf8))
    assert out.name == "j_call" and out.dtype == pl.Utf8
