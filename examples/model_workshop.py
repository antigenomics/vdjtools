# vdjtools — the recombination-model workshop, end to end.
# Build a model on a custom germline library, fit it to sequences, check it, compare it against a
# reference model, score a held-out set (log-likelihood / BIC), estimate the diversity it describes,
# extend its allele library, and re-weight its V/J usage. Everything runs offline in seconds on a
# toy locus; toggle "use real reads" to pull a small pre-annotated TRB/TRA subset from HuggingFace.
# Run with:  marimo edit examples/model_workshop.py
import marimo

__generated_with = "0.23.14"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import subprocess

    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    from vdjtools.model import check_model, load_bundled, load_germline, rescale_usage
    from vdjtools.model.analyze import (
        compare_models,
        compare_net_dot,
        compare_usage,
        total_entropy,
    )
    from vdjtools.model.generate import generate
    from vdjtools.model.infer import extend_alleles, infer_frame, training_frame
    from vdjtools.model.io import from_germline, marginals_frame
    from vdjtools.model.score import compare_pgen, diversity, model_fit, pgen_spectrum, pgen_summary
    return (
        check_model,
        compare_models,
        compare_net_dot,
        compare_pgen,
        compare_usage,
        diversity,
        extend_alleles,
        from_germline,
        generate,
        infer_frame,
        load_bundled,
        load_germline,
        marginals_frame,
        model_fit,
        np,
        pgen_spectrum,
        pgen_summary,
        pl,
        plt,
        rescale_usage,
        subprocess,
        total_entropy,
        training_frame,
    )


@app.cell
def _(mo):
    mo.md(
        """
        # The recombination-model workshop

        A model here is a **manifest** (the recombination Bayes net), a set of **tidy polars
        marginal tables**, and the **germline** those tables are keyed against. Because the
        probabilities are ordinary DataFrames, every step below is table-in / table-out.

        We build one from scratch on a three-gene toy locus, fit it, check it, compare it, score
        sequences under it, and ask how much diversity it actually describes.
        """
    )
    return


@app.cell
def _(mo):
    mo.md("## 1 · A custom germline library")
    return


@app.cell
def _(pl):
    # A hand-written library. `sequence` is the CDR3-REGION germline: V from the conserved Cys104
    # codon to the 3' end, J from the 5' end through the [FW]118 codon. Real work would use
    # `read_germline_fasta("V.fasta", "J.fasta")`, or `load_germline("TRB", "human")` for arda's.
    toy_germline = pl.DataFrame([
        {"allele": "TOYV1*01", "segment": "V", "sequence": "TGTGCCAGCAGCTTA"},
        {"allele": "TOYV1*02", "segment": "V", "sequence": "TGTGCCAGCAGCTTG"},
        {"allele": "TOYV2*01", "segment": "V", "sequence": "TGTGCTTCCTCACTG"},
        {"allele": "TOYJ1*01", "segment": "J", "sequence": "AACACTGAAGCTTTCTTT"},
        {"allele": "TOYJ2*01", "segment": "J", "sequence": "AACGAGCAGTACTTT"},
    ])
    toy_germline
    return (toy_germline,)


@app.cell
def _(mo, toy_germline):
    from vdjtools.model import validate_germline

    issues = validate_germline(toy_germline)
    mo.vstack([
        mo.md("`validate_germline` audits the library **before** it becomes a model. The anchor-frame "
              "checks matter most: an anchor one codon off shifts every deletion profile by a "
              "constant and nothing downstream complains."),
        issues if issues.height else mo.md("*No issues — the library is clean.*").callout("success"),
    ])
    return


@app.cell
def _(from_germline, toy_germline):
    # A D allele in the frame would make this VDJ; without one it is VJ. The marginals start as
    # placeholders, and their support ranges bound what EM can later learn -- hence ins_max.
    template = from_germline(toy_germline, locus="TOY", ins_max=6)
    template
    return (template,)


@app.cell
def _(mo):
    mo.md(
        """
        ## 2 · Training data

        Real models are fitted to **non-functional** reads — out-of-frame *or* stop-codon, since
        both escaped selection. Keeping only the out-of-frame half would condition the training set
        on junction length modulo 3, which the insertion-length model would then happily learn.

        Here we simulate from a "truth" model so the fit can be checked against a known answer.
        """
    )
    return


@app.cell
def _(mo):
    use_real = mo.ui.checkbox(
        value=False,
        label="Use the shipped real TRB reads instead (arda-mapped, in the source tree, offline)",
    )
    use_real
    return (use_real,)


@app.cell
def _(generate, mo, pl, template, use_real):
    # Truth = the template with a deliberately non-uniform V usage and a peaked insertion length,
    # so "did EM recover it" is a question with an answer.
    from vdjtools.model.model import Model

    truth = Model(
        manifest=template.manifest,
        tables={
            **template.tables,
            "v_choice": pl.DataFrame({"v_allele": ["TOYV1*01", "TOYV1*02", "TOYV2*01"],
                                      "p": [0.6, 0.1, 0.3]}),
            "vj_ins": pl.DataFrame({"length": pl.Series([0, 1, 2, 3, 4, 5, 6], dtype=pl.Int16),
                                    "p": [0.05, 0.10, 0.30, 0.30, 0.15, 0.07, 0.03]}),
        },
        genomic=template.genomic,
    )
    clones = generate(truth, 4000, seed=0).rename({"junction_nt": "junction"}).select(
        ["junction", "v_call", "j_call"])
    note = mo.md("")
    if use_real.value:
        try:
            from vdjtools.model.data import load_prepared

            clones = load_prepared("human", "TRB", "nonfunctional")
            note = mo.md(
                f"Using **{clones.height:,}** real out-of-frame TRB clonotypes — 5'RACE reads "
                f"mapped with arda, shipped in `tests/python/fixtures/model_reads/`. There is no "
                f"ground-truth model for these, so the recovery check below stops applying."
            ).callout("success")
        except FileNotFoundError as e:   # running from an installed wheel, not a checkout
            note = mo.md(f"Shipped reads unavailable ({e}); staying on the simulated set.").callout("warn")
    mo.vstack([note, clones.head(5)])
    return clones, truth


@app.cell
def _(mo):
    mo.md("## 3 · Fit, and read the training log")
    return


@app.cell
def _(clones, infer_frame, template, training_frame):
    fitted, report = infer_frame(template, clones, max_iter=12, tol=1e-5)
    log = training_frame(fitted)
    log
    return fitted, log, report


@app.cell
def _(log, mo, plt, report):
    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(log["iter"], log["loglik"], marker="o")
    ax.set_xlabel("EM iteration")
    ax.set_ylabel("mean log Pgen per sequence")
    ax.set_title(f"converged={report.converged} after {report.n_iter} iterations")
    fig.tight_layout()
    mo.vstack([
        mo.md("The log is stored **on the model** (`model.training`) and saved beside it as "
              "`training.json`, so a model can always say what it was fitted on. A warm-start "
              "refit appends a second run rather than overwriting the first."),
        ax,
    ])
    return


@app.cell
def _(fitted, mo, pl, truth):
    recovered = (
        truth.tables["v_choice"].rename({"p": "truth"})
        .join(fitted.tables["v_choice"].rename({"p": "fitted"}), on="v_allele")
        .with_columns(pl.col("fitted").round(3), pl.col("truth").round(3))
    )
    mo.vstack([mo.md("### Did EM recover the V usage it was given?"), recovered])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 4 · Check the model

        `check_model` returns a tidy issue frame instead of raising, so every problem is visible at
        once. `severity == "error"` means the model will crash or score wrongly; `warn` means
        suspicious but usable. Each check exists because it once produced a silently wrong answer.
        """
    )
    return


@app.cell
def _(check_model, fitted, mo, pl):
    fit_issues = check_model(fitted, germline="none")
    errors = fit_issues.filter(pl.col("severity") == "error")
    mo.vstack([
        mo.md(f"**{errors.height}** error(s), "
              f"**{fit_issues.filter(pl.col('severity') == 'warn').height}** warning(s)"),
        fit_issues if fit_issues.height else mo.md("*Clean.*").callout("success"),
    ])
    return


@app.cell
def _(check_model, load_bundled, mo, pl):
    # The same check on a shipped model. The findings here are real and worth knowing: OLGA leaves a
    # few ORF alleles with no CDR3-region germline while still giving them usage, so their Pgen is
    # always 0; and a shared deletion-bin grid strands some mass on the shortest alleles.
    real_issues = check_model(load_bundled("TRB", "learned"))
    mo.vstack([
        mo.md("### The same audit on the bundled human TRB model"),
        real_issues.group_by(["severity", "check"]).len().sort(["severity", "len"],
                                                               descending=[False, True]),
        real_issues.filter(pl.col("severity") != "info").head(5),
    ])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 5 · Information content and total diversity

        Two different questions. **Scenario entropy** is the information in one recombination
        *event*, summed over the Bayes net. **Sequence entropy** is the entropy of the junction
        distribution itself — always smaller, because different scenarios can produce the same
        junction — estimated by Monte Carlo (sequences drawn from the model *are* distributed as
        Pgen, so `E[-log2 Pgen]` is unbiased and the standard error comes free).
        """
    )
    return


@app.cell
def _(load_bundled, mo, total_entropy):
    trb = load_bundled("TRB", "olga")
    per_event = total_entropy(trb)
    mo.vstack([
        mo.md(f"### Where human TRB's **{per_event['contribution_bits'].sum():.1f} bits** "
              f"per rearrangement come from"),
        per_event.sort("contribution_bits", descending=True),
    ])
    return (trb,)


@app.cell
def _(diversity, mo, trb):
    div = diversity(trb, n=3000, seed=0)
    d = div.to_dicts()[0]
    mo.vstack([
        div,
        mo.md(
            f"- scenario entropy **{d['scenario_entropy_bits']:.1f} bits** "
            f"→ {d['scenario_diversity']:.3g} distinct rearrangements\n"
            f"- sequence entropy **{d['sequence_entropy_bits']:.2f} ± "
            f"{d['sequence_entropy_se_bits']:.2f} bits**\n"
            f"- Hill *q*=1 (`2^H`) — the usual headline figure — **{d['diversity_shannon']:.3g}** "
            f"distinct junctions\n"
            f"- Hill *q*=2 (`1/E[Pgen]`) — draws before two coincide — "
            f"**{d['diversity_simpson']:.3g}**"
        ),
    ])
    return


@app.cell
def _(mo, pgen_spectrum, plt, trb):
    spec = pgen_spectrum(trb, n=3000, seed=0, bins=30)
    fig2, ax2 = plt.subplots(figsize=(6, 3))
    ax2.bar(spec["bin_mid"], spec["frac"], width=(spec["bin_right"] - spec["bin_left"])[0])
    ax2.set_xlabel("log10 Pgen")
    ax2.set_ylabel("fraction of generated sequences")
    ax2.set_title("TRB Pgen spectrum")
    fig2.tight_layout()
    mo.vstack([mo.md("### The Pgen distribution the model implies"), ax2])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 6 · Compare two models

        OLGA's TRB model was fit to DNA-multiplex data; the `learned` one to 5'RACE reads. They
        describe the same recombination machinery but were amplified differently, so **V usage
        should disagree and neither is wrong** — while the junction model (trims, insertions,
        dinucleotides) is the genuinely comparable part.
        """
    )
    return


@app.cell
def _(compare_models, load_bundled, mo, trb):
    learned_trb = load_bundled("TRB", "learned")
    diff = compare_models(trb, learned_trb, labels=("olga", "learned"), by="gene")
    mo.vstack([
        mo.md("Jensen-Shannon is the headline: symmetric, bounded by one bit, and **finite when the "
              "supports differ** — which is exactly the case here. `tv_max` reports the worst single "
              "group, which finds the one broken gene an average hides."),
        diff.sort("jsd_bits", descending=True, nulls_last=True),
    ])
    return (learned_trb,)


@app.cell
def _(compare_net_dot, learned_trb, mo, subprocess, trb):
    # Render the comparison graph inline via the graphviz `dot` CLI (no python-graphviz dep).
    cdot = compare_net_dot(trb, learned_trb, labels=("olga", "learned"))
    try:
        out = subprocess.run(["dot", "-Tsvg"], input=cdot, capture_output=True, text=True,
                             check=True).stdout
        graph = mo.Html(out[out.index("<svg"):])
    except (FileNotFoundError, subprocess.CalledProcessError):
        graph = mo.md("`dot` (graphviz) not found — install it to see the comparison graph.").callout("warn")
    mo.vstack([mo.md("### Comparison graph &nbsp; *(node shade = JSD; dashed/dotted = one side only)*"),
               graph])
    return


@app.cell
def _(compare_usage, learned_trb, mo, plt, trb):
    usage = compare_usage(trb, learned_trb, "v").head(20)
    fig3, ax3 = plt.subplots(figsize=(7, 4))
    y = range(usage.height)
    ax3.barh([i + 0.2 for i in y], usage["p_a"], height=0.4, label="olga")
    ax3.barh([i - 0.2 for i in y], usage["p_b"], height=0.4, label="learned")
    ax3.set_yticks(list(y)); ax3.set_yticklabels(usage["name"].to_list())
    ax3.invert_yaxis(); ax3.set_xlabel("P(V gene)"); ax3.legend()
    fig3.tight_layout()
    mo.vstack([mo.md("### V usage is protocol-dependent, not a defect"), ax3])
    return


@app.cell
def _(mo):
    mo.md("## 7 · Score sequences: Pgen distributions, log-likelihood and BIC")
    return


@app.cell
def _(compare_pgen, generate, learned_trb, mo, pgen_summary, trb):
    # A held-out set, drawn from neither model's training data.
    held_out = generate(trb, 600, seed=99)["junction_nt"].to_list()
    cmp = compare_pgen(trb, learned_trb, held_out, labels=("olga", "learned"), use_calls=False)
    summary = pgen_summary(cmp, labels=("olga", "learned"))
    mo.vstack([
        mo.md("`only_a_scoreable` / `only_b_scoreable` are the numbers that usually matter: one "
              "model assigning Pgen 0 to sequences the other scores fine is the finding, and a "
              "mean-delta-only report would hide it entirely."),
        summary,
    ])
    return cmp, held_out


@app.cell
def _(cmp, mo, np, plt):
    a = cmp["log10_olga"].to_numpy()
    b = cmp["log10_learned"].to_numpy()
    ok = ~(np.isnan(a) | np.isnan(b))
    fig4, ax4 = plt.subplots(figsize=(4.5, 4.5))
    ax4.scatter(a[ok], b[ok], s=6, alpha=0.4)
    lim = [min(a[ok].min(), b[ok].min()), max(a[ok].max(), b[ok].max())]
    ax4.plot(lim, lim, "k--", lw=1)
    ax4.set_xlabel("log10 Pgen (olga)"); ax4.set_ylabel("log10 Pgen (learned)")
    fig4.tight_layout()
    mo.vstack([mo.md("### The same sequences under both models"), ax4])
    return


@app.cell
def _(held_out, learned_trb, mo, model_fit, pl, trb):
    fits = pl.concat([
        model_fit(trb, held_out, use_calls=False).with_columns(model=pl.lit("olga")),
        model_fit(learned_trb, held_out, use_calls=False).with_columns(model=pl.lit("learned")),
    ]).select(["model", "n", "n_scoreable", "loglik_sum", "k", "aic", "bic"])
    mo.vstack([
        mo.md(
            "Likelihoods use **nucleotide** Pgen, because `Σ Pgen_nt = 1` over all nt CDR3s — so "
            "`log Pgen_nt` is a proper log-likelihood and BIC means something. Amino-acid Pgen sums "
            "only the in-frame, stop-free fiber of a translation, so its log-likelihood is "
            "unnormalized and the missing constant *differs between models*.\n\n"
            "The free-parameter count `k` is **support-based**: occupied cells minus one per "
            "normalization group, dropping undefined and unreachable ones. Counting rows instead "
            "would put TRB's `v_3_del` at ~3,600 parameters when ~700 are real."
        ),
        fits,
        mo.md("NOTE: These sequences were drawn from `olga`, so it *should* win here. A log-likelihood "
              "on a model's own training data is that model's EM objective, which EM increases by "
              "construction — it validates nothing. Always score a held-out set.").callout("warn"),
    ])
    return


@app.cell
def _(mo):
    mo.md("## 8 · Extend the allele library, and re-weight usage")
    return


@app.cell
def _(extend_alleles, learned_trb, load_germline, mo, pl):
    bigger = extend_alleles(learned_trb, load_germline("TRB", "human"))

    def gene_usage(m):
        from vdjtools.model.analyze import gene_marginal

        out = {}
        for allele, p in gene_marginal(m, "v").items():
            out[allele.split("*")[0]] = out.get(allele.split("*")[0], 0.0) + p
        return out

    before, after = gene_usage(learned_trb), gene_usage(bigger)
    drift = pl.DataFrame([{"gene": g, "before": before[g], "after": after.get(g, 0.0)}
                          for g in before]).with_columns(
        delta=(pl.col("after") - pl.col("before")).abs()).sort("delta", descending=True)
    mo.vstack([
        mo.md(f"V alleles **{learned_trb.genomic['genes_v'].height} → "
              f"{bigger.genomic['genes_v'].height}**. Each pre-existing *gene* keeps its total "
              f"usage: alleles of one gene are alternative versions of the same gene — a diploid "
              f"carries at most two — so a richer library splits a gene's mass more finely rather "
              f"than multiplying it."),
        drift.head(8),
        mo.md("This *seeds*; it does not estimate. Follow it with "
              "`infer_frame(bigger, clones, init='template')`."),
    ])
    return


@app.cell
def _(generate, learned_trb, mo, plt, rescale_usage, trb):
    # V/J usage is protocol-dependent, the junction model is not. Learn the junction model once,
    # then set P(V) from the repertoire you are actually about to score.
    sample = generate(trb, 3000, seed=7)
    rescaled = rescale_usage(learned_trb, sample)
    from vdjtools.model.analyze import compare_usage as _cu

    shift = _cu(learned_trb, rescaled, "v").head(15)
    fig5, ax5 = plt.subplots(figsize=(7, 4))
    yy = range(shift.height)
    ax5.barh([i + 0.2 for i in yy], shift["p_a"], height=0.4, label="learned (5'RACE)")
    ax5.barh([i - 0.2 for i in yy], shift["p_b"], height=0.4, label="rescaled to the sample")
    ax5.set_yticks(list(yy)); ax5.set_yticklabels(shift["name"].to_list())
    ax5.invert_yaxis(); ax5.set_xlabel("P(V gene)"); ax5.legend()
    fig5.tight_layout()
    mo.vstack([mo.md("### `rescale_usage` — same junction model, your protocol's V usage"), ax5])
    return


@app.cell
def _(mo):
    mo.md("## 9 · Every probability as one table")
    return


@app.cell
def _(fitted, marginals_frame, mo):
    flat = marginals_frame(fitted)
    mo.vstack([
        mo.md("`marginals_frame` flattens every marginal into one self-describing long table, and "
              "`set_marginals` takes it back — so a hand-edited TSV is a first-class model input. "
              "`save_model(m, path, fmt='tsv')` writes a whole model directory you can read in any "
              "tool; `load_model` detects the format and restores the dtypes."),
        flat.head(12),
    ])
    return


@app.cell
def _(mo):
    mo.md(
        """
        ---
        **Building from the real read corpus.** The bundled `learned` models are fitted to raw
        5'RACE FASTQ from the private `isalgo/airr_model_read` dataset: fetch, map with arda,
        collapse to unique clonotypes, then EM. That whole pipeline is one command, parallel across
        chains (it needs HuggingFace access and arda's mmseqs2, and takes minutes per chain):

        ```bash
        vdjtools model build --chains TRB,TRA,IGH --workers 4 -o models/
        ```
        """
    )
    return


if __name__ == "__main__":
    app.run()
