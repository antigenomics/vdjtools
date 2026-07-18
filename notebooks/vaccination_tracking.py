# vdjtools — longitudinal clonotype tracking across a vaccination time course.
# Reactive marimo app over three public vaccination cohorts (yellow fever `isalgo/airr_yfv19`,
# influenza `isalgo/airr_flu_vac`, TBE virus `isalgo/airr_tbev_vac`): pick a vaccine + donor,
# track every clonotype across the time course, classify each pre->post pair with the
# effective-sample-size paired test (`vdjtools.dynamics.test_pair`), and read off the VDJtrack
# recapture model (size buckets x Poisson capture). Run with:
#     marimo edit notebooks/vaccination_tracking.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _(mo):
    mo.md(
        r"""
        # Clonotype tracking in vaccination — a vdjtools v2 explorer

        A vaccine drives a **clonal expansion**: antigen-specific T-cell clones divide, so their
        clonotypes rise in frequency around the response peak (day ~15 for yellow fever) and a
        subset **persist** into memory. `vdjtools.dynamics` makes that measurable within one
        donor:

        - **`test_pair`** — the Ayestaran (2024) effective-sample-size paired test: it reads each
          pair's own noise off its mean–variance scaling (`N_eff`), so it does not manufacture
          hits the way a naive Fisher-on-read-counts does. Every clonotype is called *emergent /
          expanded / persistent / contracted / vanishing*.
        - **the recapture model** (VDJtrack; Pavlova, Zvyagin & Shugay 2024) — under Poisson
          sampling a clone of frequency `f` is recaptured with probability `1 − exp(−f·R)`, so
          recapture rises with clone **size** (singleton → large). Binning by size and comparing
          an annotated group (here: the clones that expanded at the peak) to that baseline shows
          whether the *response* persists more than background.

        Data auto-loads from HuggingFace (`isalgo/airr_yfv19` · `airr_flu_vac` · `airr_tbev_vac`),
        preferring a local `~/hf/` or `./` copy if present.
        """
    )
    return


@app.cell
def _():
    # --- imports & configuration (single cell so every name is defined once) ---
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from vdjtools.io.batch import read
    from vdjtools.overlap.track import track_clonotypes
    from vdjtools.dynamics import SIZE_CLASSES, capture_rates, capture_test, test_pair
    from vdjtools.dynamics.paired import CLASSES

    # Each vaccine: HF repo, metadata sheet, file-name folder prefix, the column giving the
    # per-person subject id, and a (column, value) filter selecting one library type per sample.
    VACCINES = {
        "Yellow fever (YFV)": dict(repo="isalgo/airr_yfv19", meta="metadata.txt", folder="",
                                   subject="donor", day="day", filt=("replica", "F1")),
        "Influenza": dict(repo="isalgo/airr_flu_vac", meta="metadata.tsv", folder="samples/",
                          subject="donor", day="day", filt=("cell_subset", "bulk")),
        "TBE virus": dict(repo="isalgo/airr_tbev_vac", meta="metadata.tsv", folder="samples/",
                          subject="donor", day="day", filt=("cell_subset", "bulk")),
    }
    # Colour-blind-safe (Okabe–Ito) per dynamics class.
    CLASS_COLOR = {"emergent": "#009E73", "expanded": "#D55E00", "persistent": "#8C8C8C",
                   "contracted": "#0072B2", "vanishing": "#CC79A7", "untested": "#E5E5E5"}

    def fetch(repo, filename):
        """Local-first data fetch: ./<file>, ~/hf/<repo-basename>/<file>, else HuggingFace."""
        for root in (Path.cwd(), Path.home() / "hf" / repo.split("/")[-1]):
            p = root / filename
            if p.exists():
                return str(p)
        import huggingface_hub as hub
        return hub.hf_hub_download(repo, filename, repo_type="dataset")

    return (CLASSES, CLASS_COLOR, Path, SIZE_CLASSES, VACCINES, capture_rates, capture_test,
            fetch, mo, np, pl, plt, read, test_pair, track_clonotypes)


@app.cell
def _(VACCINES, mo):
    vaccine = mo.ui.dropdown(list(VACCINES), value="Yellow fever (YFV)", label="Vaccine")
    vaccine
    return (vaccine,)


@app.cell
def _(VACCINES, fetch, mo, pl, vaccine):
    # Read the chosen vaccine's metadata; offer its subjects (donors) with a full time course.
    cfg = VACCINES[vaccine.value]
    _mp = fetch(cfg["repo"], cfg["meta"])
    meta = pl.read_csv(_mp, separator="\t", infer_schema_length=0)
    fcol, fval = cfg["filt"]
    if fcol in meta.columns:
        meta = meta.filter(pl.col(fcol).str.contains(fval))
    meta = meta.with_columns(pl.col(cfg["day"]).cast(pl.Int64, strict=False).alias("_day"))
    subjects = (meta.group_by(cfg["subject"]).agg(pl.col("_day").n_unique().alias("n"))
                .filter(pl.col("n") >= 3).sort("n", descending=True))
    donor = mo.ui.dropdown(subjects[cfg["subject"]].to_list(),
                           value=subjects[cfg["subject"]][0], label="Donor")
    mo.vstack([mo.md(f"**{vaccine.value}** — {subjects.height} donors with ≥3 timepoints."), donor])
    return cfg, donor, meta


@app.cell
def _(cfg, donor, fetch, meta, mo, pl, read):
    # Load the donor's time course: {day: clonotype frame}, ordered by day.
    _sub = meta.filter(pl.col(cfg["subject"]) == donor.value).sort("_day")
    frames, order = {}, []
    for _row in _sub.iter_rows(named=True):
        _dy = int(_row["_day"])
        frames[_dy] = read(fetch(cfg["repo"], cfg["folder"] + _row["file_name"]))
        order.append(_dy)
    sizes = pl.DataFrame({"day": order,
                          "clonotypes": [frames[d].height for d in order],
                          "reads": [int(frames[d]["duplicate_count"].sum()) for d in order]})
    mo.vstack([mo.md(f"Loaded **{len(order)} timepoints** for donor **{donor.value}** "
                     f"(days {order})."), sizes])
    return frames, order, sizes


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 1 · Track every clonotype, classify each pre→post transition

        `track_clonotypes` pivots per-sample frequency into one column per timepoint;
        `test_pair` then classifies each clonotype between the **baseline** (first timepoint) and
        each later timepoint. The stacked area shows the read-fraction in each dynamics class over
        the course — the vaccine peak is the emergent/expanded swell.
        """
    )
    return


@app.cell
def _(CLASSES, frames, order, pl, test_pair, track_clonotypes):
    tracked = track_clonotypes(frames, order=[str(d) for d in order])
    # Classify baseline (first day) vs every later timepoint; tabulate read-fraction per class.
    base = order[0]
    comp = {}
    for d in order[1:]:
        res = test_pair(frames[base], frames[d])
        w = res.with_columns(pl.col("count_b").alias("_w"))
        frac = (w.group_by("dynamics").agg(pl.col("_w").sum())
                .with_columns((pl.col("_w") / pl.col("_w").sum()).alias("frac")))
        comp[d] = {r["dynamics"]: r["frac"] for r in frac.iter_rows(named=True)}
    class_frac = pl.DataFrame(
        [{"day": d, **{c: comp[d].get(c, 0.0) for c in CLASSES}} for d in order[1:]])
    return base, class_frac, comp, tracked


@app.cell
def _(CLASS_COLOR, CLASSES, base, class_frac, np, order, plt):
    _days = class_frac["day"].to_list()
    _classes = [c for c in CLASSES if c != "untested"]
    _stack = np.array([class_frac[c].to_numpy() for c in _classes])
    fig1, ax1 = plt.subplots(figsize=(7.2, 4.2))
    ax1.stackplot(_days, _stack, labels=_classes,
                  colors=[CLASS_COLOR[c] for c in _classes], alpha=0.9)
    ax1.set(xlabel=f"Day (vs baseline day {base})", ylabel="Read fraction in class",
            title="Repertoire dynamics across the time course")
    ax1.legend(loc="upper left", fontsize=8, frameon=False, ncol=2)
    ax1.set_xticks(order[1:])
    ax1.spines[["top", "right"]].set_visible(False)
    fig1.tight_layout()
    fig1
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 2 · Trajectories of the top responding clonotypes

        The clonotypes that **expanded** at the peak, followed across every timepoint — the
        classic "rise at day 15, contract into memory" vaccine trajectory.
        """
    )
    return


@app.cell
def _(base, frames, order, pl, plt, test_pair, tracked):
    _peak = order[len(order) // 2] if len(order) > 2 else order[-1]        # ~day 15
    _res = test_pair(frames[base], frames[_peak])
    _top = (_res.filter(pl.col("dynamics") == "expanded")
            .sort("q_value").head(12)["junction_aa"].to_list())
    _fcols = [f"freq_{d}" for d in order]
    _traj = tracked.filter(pl.col("junction_aa").is_in(_top))
    fig2, ax2 = plt.subplots(figsize=(7.2, 4.2))
    for _row in _traj.iter_rows(named=True):
        ax2.plot(order, [max(_row[c], 1e-7) for c in _fcols], "-o", ms=3, lw=1, alpha=0.8)
    ax2.set_yscale("log")
    ax2.set(xlabel="Day", ylabel="Clonotype frequency (log)",
            title=f"Top {len(_top)} clones expanded at day {_peak}")
    ax2.set_xticks(order)
    ax2.spines[["top", "right"]].set_visible(False)
    fig2.tight_layout()
    fig2
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ## 3 · The recapture model — do the responders persist?

        Take the clones that **expanded at the peak** as the "response" group and everything else
        as background, then measure the recapture rate from the peak into the **latest** timepoint
        per size class (singleton / doubleton / tripleton / large). Recapture rises with size for
        *both* groups (Poisson sampling), but the response group sits far above the size baseline —
        it persists. `capture_test` fits `log(recapture) ~ size + group` and reports the group
        effect.
        """
    )
    return


@app.cell
def _(capture_rates, capture_test, donor, frames, order, pl, test_pair):
    peak = order[len(order) // 2] if len(order) > 2 else order[0]
    late = order[-1]
    _src = frames[peak] if peak != late else frames[order[0]]
    _res = test_pair(_src, frames[late])
    _expanded = set(_res.filter(pl.col("dynamics").is_in(["expanded", "emergent"]))
                    ["junction_aa"].to_list())
    _pre = frames[peak].with_columns(
        pl.when(pl.col("junction_aa").is_in(list(_expanded))).then(pl.lit("response"))
        .otherwise(pl.lit("background")).alias("group"))
    rates = capture_rates(_pre, frames[late], group_col="group", donor=str(donor.value))
    coef = capture_test(rates, group_col="group")
    return coef, late, peak, rates


@app.cell
def _(SIZE_CLASSES, coef, late, mo, peak, plt, rates):
    _order = [s for s in SIZE_CLASSES if s in rates["size_class"].to_list()]
    fig3, ax3 = plt.subplots(figsize=(7.0, 4.2))
    for _grp, _color, _marker in [("background", "#8C8C8C", "o"), ("response", "#D55E00", "s")]:
        _sub = rates.filter(rates["group"] == _grp)
        _x = [_order.index(s) for s in _sub["size_class"].to_list()]
        ax3.errorbar(_x, _sub["capture_rate"].to_numpy(),
                     yerr=[(_sub["capture_rate"] - _sub["ci_lo"]).to_numpy(),
                           (_sub["ci_hi"] - _sub["capture_rate"]).to_numpy()],
                     fmt=_marker, color=_color, capsize=3, label=_grp, ms=7)
    ax3.set_xticks(range(len(_order)))
    ax3.set_xticklabels(_order)
    ax3.set(xlabel="Pre-sample size class", ylabel=f"Recapture rate (day {peak}→{late})",
            title="Response clones persist above the size baseline")
    ax3.set_yscale("log")
    ax3.legend(frameon=False)
    ax3.spines[["top", "right"]].set_visible(False)
    fig3.tight_layout()
    _g = coef.filter(coef["term"].str.starts_with("group:"))
    mo.vstack([mo.md(f"**Group effect** (log-linear, size-adjusted): "
                     f"β = {_g['estimate'][0]:.2f}, p = {_g['p_value'][0]:.2g}"), fig3])
    return


@app.cell
def _(mo):
    mo.md(
        r"""
        ---
        **Takeaway.** One `test_pair` call classifies every clonotype's within-donor change with a
        pair-specific noise model; the recapture model then asks the survival question the raw
        counts cannot. Swap the **vaccine** and **donor** above to compare yellow fever, influenza
        and TBE responses. The per-clonotype test is `vdjtools.dynamics.test_pair` (Ayestaran
        2024); the capture model and the 1-mismatch metaclonotype grouping
        (`vdjtools.dynamics.capture_rates`, `test_metaclonotypes`) port the VDJtrack pipeline
        (Pavlova, Zvyagin & Shugay, *Front Immunol* 2024).
        """
    )
    return


if __name__ == "__main__":
    app.run()
