"""Clonotype recapture model — the VDJtrack size-bucket capture test.

Where :func:`vdjtools.dynamics.test_pair` calls *each* clonotype expanded/contracted, the
capture model asks a **group** question: does an annotated set of clonotypes (e.g. antigen-
specific, or "emerging") **persist / recapture** across a time course more than the rest,
*after* controlling for what drives recapture on its own — clonotype **size** and repertoire
**diversity**?

Under Poisson sampling a clonotype of population frequency ``f`` is recaptured in a sample of
depth ``R`` with probability ``P = 1 − exp(−f·R)`` (:func:`poisson_capture`), so recapture rate
rises monotonically with pre-sample size — **singleton < doubleton < tripleton < large**. The
model bins clonotypes into those four size classes, measures the recapture fraction per
(donor, group, size class), puts a ``Beta(captured, missing)`` posterior on it for uncertainty,
and tests the **group** effect with a log-linear model ``log(recapture) ~ size + group +
log(div_ratio)`` (and a per-bucket paired test across donors).

Port of the group's VDJtrack R pipeline (``github.com/antigenomics/vdjtrack``); the method is
Pavlova, Zvyagin & Shugay, *Front Immunol* 2024, DOI 10.3389/fimmu.2024.1321603. The two-step
sampling correction and the per-clonotype test live in :mod:`vdjtools.dynamics.paired`.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import beta as _beta
from scipy.stats import t as _t
from scipy.stats import ttest_rel

from ..io.schema import COUNT
from .paired import DEFAULT_KEY

#: Pre-sample size classes, ordered (singleton=1, doubleton=2, tripleton=3, large=4+).
SIZE_CLASSES = ("singleton", "doubleton", "tripleton", "large")


def size_class(count: "str | pl.Expr" = COUNT) -> pl.Expr:
    """Polars expression classifying a clonotype count into a :data:`SIZE_CLASSES` bucket."""
    c = pl.col(count) if isinstance(count, str) else count
    return (pl.when(c == 1).then(pl.lit("singleton"))
              .when(c == 2).then(pl.lit("doubleton"))
              .when(c == 3).then(pl.lit("tripleton"))
              .otherwise(pl.lit("large")))


def poisson_capture(freq, depth) -> np.ndarray:
    """Poisson recapture probability ``P = 1 − exp(−f·R)``.

    The probability that a clonotype of population frequency ``freq`` yields at least one
    molecule in a sample of ``depth`` (total reads/UMIs) — the generative model the recapture
    curve is read against (Pavlova 2024; aging_capture_model.Rmd).

    Args:
        freq: Clonotype population frequency (scalar or array).
        depth: Sampling depth ``R`` (scalar or array).

    Returns:
        The capture probability, same shape as the broadcast of ``freq`` and ``depth``.
    """
    return 1.0 - np.exp(-np.asarray(freq, dtype=float) * np.asarray(depth, dtype=float))


def _presence(df: pl.DataFrame, key: list[str]) -> pl.DataFrame:
    """Unique clonotype keys of a sample (its recapture target set)."""
    return df.select(key).unique()


def capture_rates(pre: pl.DataFrame, post: pl.DataFrame, *, key=DEFAULT_KEY,
                  group_col: str | None = None, donor: str | None = None) -> pl.DataFrame:
    """Recapture fraction per (size class[, group]) between one ``pre`` and ``post`` sample.

    Each ``pre`` clonotype is binned by its pre-sample size (:data:`SIZE_CLASSES`) and marked
    **captured** if its ``key`` reappears in ``post``. Rows are then counted per size class
    (and ``group_col`` level, if given): ``n_captured`` = α, ``n_total`` = α+β, and the recapture
    rate gets a ``Beta(α+1, β+1)`` posterior (Laplace) with a 95% credible interval.

    Args:
        pre: The earlier sample (canonical clonotype frame).
        post: The later sample.
        key: Clonotype match key (default CDR3 aa + V + J).
        group_col: A column of ``pre`` carrying a group label (e.g. antigen specificity). If
            ``None``, all clonotypes are one group ``"all"``.
        donor: If given, added as a constant ``donor`` column (for :func:`capture_test` over a
            cohort assembled by concatenating per-donor calls).

    Returns:
        One row per (``donor``,) ``group``, ``size_class`` with ``n_captured``, ``n_total``,
        ``capture_rate`` (posterior mean), ``ci_lo``, ``ci_hi``.
    """
    key = list(key)
    keep = [*key, *( [group_col] if group_col else [] )]
    agg = pre.group_by(keep, maintain_order=True).agg(pl.col(COUNT).sum().alias("_c"))
    agg = agg.with_columns(size_class("_c").alias("size_class"))
    captured = agg.join(_presence(post, key).with_columns(pl.lit(True).alias("_cap")),
                        on=key, how="left").with_columns(pl.col("_cap").fill_null(False))
    grp = ["size_class"] + ([group_col] if group_col else [])
    out = (captured.group_by(grp, maintain_order=True)
           .agg(pl.col("_cap").sum().cast(pl.Int64).alias("n_captured"),
                pl.len().alias("n_total")))
    alpha = out["n_captured"].to_numpy() + 1.0        # Beta(α+1, β+1), Laplace prior
    beta = (out["n_total"].to_numpy() - out["n_captured"].to_numpy()) + 1.0
    out = out.with_columns(
        pl.Series("capture_rate", alpha / (alpha + beta)),
        pl.Series("ci_lo", _beta.ppf(0.025, alpha, beta)),
        pl.Series("ci_hi", _beta.ppf(0.975, alpha, beta)),
    )
    if group_col is None:
        out = out.with_columns(pl.lit("all").alias("group"))
    else:
        out = out.rename({group_col: "group"})
    if donor is not None:
        out = out.with_columns(pl.lit(donor).alias("donor"))
    cols = (["donor"] if donor is not None else []) + ["group", "size_class",
            "n_captured", "n_total", "capture_rate", "ci_lo", "ci_hi"]
    return out.select(cols)


def _ols(y: np.ndarray, X: np.ndarray, names: list[str]) -> pl.DataFrame:
    """Ordinary least squares with per-coefficient t-tests (no statsmodels)."""
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(y) - X.shape[1]
    resid = y - X @ coef
    if dof <= 0:                                       # saturated / underdetermined
        se = np.full(X.shape[1], np.nan)
        p = np.full(X.shape[1], np.nan)
    else:
        sigma2 = float(resid @ resid) / dof
        cov = np.linalg.pinv(X.T @ X) * sigma2
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
        with np.errstate(divide="ignore", invalid="ignore"):
            tstat = np.where(se > 0, coef / se, np.nan)
        p = 2.0 * _t.sf(np.abs(tstat), dof)
    return pl.DataFrame({"term": names, "estimate": coef, "std_error": se, "p_value": p})


def capture_test(rates: pl.DataFrame, *, group_col: str = "group",
                 div_ratio: dict | None = None) -> pl.DataFrame:
    """Test the **group** effect on recapture with a log-linear model.

    Fits ``log(capture_rate) ~ size_class + group [+ log(div_ratio)]`` by OLS (the VDJtrack
    ``lm``; example.Rmd), so the ``group`` coefficient's p-value asks whether the annotated
    group recaptures differently than baseline *after* accounting for clonotype size and (if
    supplied) the per-donor diversity ratio. Size classes and group are dummy-coded (first level
    = reference); a numeric ``log(div_ratio)`` is added per donor when ``div_ratio`` is given.

    Args:
        rates: Stacked :func:`capture_rates` output (one row per donor × group × size class).
        group_col: The group column (default ``"group"``); coded against its first level.
        div_ratio: Optional ``{donor: div_after/div_before}`` — adds ``log(div_ratio)`` as a
            covariate (needs a ``donor`` column in ``rates``).

    Returns:
        A coefficient table (``term``, ``estimate``, ``std_error``, ``p_value``); the ``group``
        row is the effect of interest.
    """
    df = rates.filter(pl.col("capture_rate") > 0)
    y = np.log(df["capture_rate"].to_numpy())
    cols = [np.ones(len(y))]
    names = ["intercept"]
    # size-class dummies (reference = first present in canonical order)
    present = [s for s in SIZE_CLASSES if s in df["size_class"].unique().to_list()]
    for s in present[1:]:
        cols.append((df["size_class"] == s).to_numpy().astype(float))
        names.append(f"size:{s}")
    # group dummies (reference = first level)
    glevels = df[group_col].unique(maintain_order=True).to_list()
    for g in glevels[1:]:
        cols.append((df[group_col] == g).to_numpy().astype(float))
        names.append(f"group:{g}")
    if div_ratio is not None:
        if "donor" not in df.columns:
            raise ValueError("div_ratio needs a 'donor' column in rates")
        lr = np.log(np.array([div_ratio[d] for d in df["donor"].to_list()], dtype=float))
        cols.append(lr)
        names.append("log_div_ratio")
    return _ols(y, np.column_stack(cols), names)


def capture_paired_test(rates: pl.DataFrame, *, group_col: str = "group",
                        donor_col: str = "donor") -> pl.DataFrame:
    """Per-size-class paired t-test of log recapture, group vs baseline across donors.

    The VDJtrack per-bucket test (vaccination.Rmd): within each size class, pair donors and
    compare the two group levels' ``log(capture_rate)`` with a paired t-test — robust to the
    donor-to-donor baseline that the log-linear model absorbs into an intercept. Requires
    exactly two group levels.

    Args:
        rates: Stacked :func:`capture_rates` output.
        group_col: The two-level group column.
        donor_col: The donor column to pair on.

    Returns:
        One row per size class: ``size_class``, ``n_pairs``, ``t``, ``p_value`` (null where
        fewer than two complete donor pairs exist).
    """
    levels = rates[group_col].unique(maintain_order=True).to_list()
    if len(levels) != 2:
        raise ValueError(f"paired test needs exactly two {group_col} levels, got {levels}")
    a, b = levels
    rows = []
    for sc in [s for s in SIZE_CLASSES if s in rates["size_class"].unique().to_list()]:
        sub = rates.filter(pl.col("size_class") == sc)
        wide = (sub.pivot(values="capture_rate", index=donor_col, on=group_col)
                .drop_nulls([a, b])) if sub.height else sub
        n = wide.height if sub.height else 0
        if n >= 2:
            res = ttest_rel(np.log(wide[a].to_numpy()), np.log(wide[b].to_numpy()))
            rows.append((sc, n, float(res.statistic), float(res.pvalue)))
        else:
            rows.append((sc, n, None, None))
    return pl.DataFrame(rows, schema=["size_class", "n_pairs", "t", "p_value"], orient="row")


def _demo() -> None:
    """Self-check: a planted specific-group persistence signal is recovered."""
    rng = np.random.default_rng(0)
    from ..io.schema import J_CALL, JUNCTION_AA, V_CALL

    def sample(cdrs, counts):
        return pl.DataFrame({JUNCTION_AA: cdrs, V_CALL: ["TRBV9"] * len(cdrs),
                             J_CALL: ["TRBJ2-3"] * len(cdrs), COUNT: counts})

    all_rates = []
    div = {}
    for d in range(6):
        n = 400
        cdrs = [f"CASS{i:04d}F" for i in range(n)]
        # specific = first 100; they persist (captured) far more often than background
        pre_counts = rng.integers(1, 6, n)
        pre = sample(cdrs, pre_counts).with_columns(
            pl.Series("specific", ["yes"] * 100 + ["no"] * 300))
        p_cap = np.where(np.arange(n) < 100, 0.9, 0.3)
        keep = rng.random(n) < p_cap
        post = sample([c for c, k in zip(cdrs, keep) if k],
                      [int(x) for x, k in zip(pre_counts, keep) if k])
        all_rates.append(capture_rates(pre, post, group_col="specific", donor=f"D{d}"))
        div[f"D{d}"] = 1.0
    rates = pl.concat(all_rates)
    coef = capture_test(rates, group_col="group")
    g = coef.filter(pl.col("term").str.starts_with("group:"))
    # sign depends on which level became the dummy reference; magnitude + p is the signal
    assert abs(g["estimate"][0]) > 0.5 and g["p_value"][0] < 0.05, coef
    paired = capture_paired_test(rates)
    assert (paired.drop_nulls("p_value")["p_value"] < 0.05).any(), paired
    assert abs(poisson_capture(1e-5, 1e6) - (1 - np.exp(-10))) < 1e-12
    print("capture._demo OK:\n", coef, "\n", paired)


if __name__ == "__main__":
    _demo()
