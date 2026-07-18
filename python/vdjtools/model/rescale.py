"""Rescale a model's V/J usage to a sample's own protocol, keeping the junction model.

V and J usage are **protocol-dependent**, the junction model is not. 5'RACE and DNA-multiplex
amplify different V genes at very different rates -- in the bundled `learned` TRB model's own
5'RACE training reads TRBV19 is 37% of functional reads, against OLGA's DNA-multiplex 3.1% --
so neither usage marginal is "right" in general, and shipping one as universal silently imposes
one protocol's amplification bias on everybody else's data.

The recombination machinery underneath *is* shared: trims, insertion lengths, insertion
dinucleotides, D usage and P(J) come out close between the two (measured TV 0.09-0.19 per
parent). That is the reusable part, and it is what these reads are good for -- they cover every
V/J combination, so the junction model and J|V are learnable from them even where the usage is
skewed.

So: learn the junction model once, then set P(V) (and P(J|V) on a VJ locus) from **your own**
sample before scoring it.

    from vdjtools.model import load_bundled, rescale_usage
    m = rescale_usage(load_bundled("TRB", "learned"), my_sample)
"""
from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable

import polars as pl

from ..io.schema import J_CALL, V_CALL
from .model import Model

Frame = "pl.DataFrame | pl.LazyFrame"


def _gene(col: str) -> pl.Expr:
    return pl.col(col).str.split("*").list.first()


def _vote_counts(sample, call_col: str) -> dict[str, float]:
    """Raw per-gene fractional votes from one frame (``DataFrame`` **or** ``LazyFrame``).

    Gene, not allele: allele calls on short reads are mismapping-prone (arda splits TRBV19 into
    *03/*01 on 151bp reads), so allele-resolution usage is noise. One vote per clonotype (row) rather
    than per read: usage is a property of the recombination, and read counts additionally carry
    clonal expansion and PCR.

    Ambiguous comma-separated calls (``IGHV3-23*01,IGHV3-23D*01``) are split FRACTIONALLY, each named
    gene getting ``1/k`` of the clonotype's vote. Dropping them is not missing-at-random: duplicated
    loci (IGHV3-23/IGHV3-23D, IGHV1-69/IGHV1-69D) are near-identical, so the aligner ties on them far
    more than on unique genes -- measured, 91% of real IGHV3-23 calls are ties, so dropping them
    understates P(IGHV3-23), the most common human IGHV gene, by 10x, and that error multiplies
    straight into Pgen (P(V) is a root factor). Fractional allocation is the honest estimate under
    "the truth is one of these, equally likely a priori".

    Only the small per-gene grouped result is materialized — the frame's other columns are never read
    — so this scales to a whole dataset streamed one sample at a time (see :func:`_dataset_usage`).
    """
    lf = sample.lazy() if isinstance(sample, pl.DataFrame) else sample
    genes = (pl.col(call_col).str.split(",")
             .list.eval(pl.element().str.strip_chars().str.split("*").list.first()))
    g = (lf.filter(pl.col(call_col).is_not_null())
           .select(genes.alias("g"))
           .with_columns(w=pl.lit(1.0) / pl.col("g").list.len())
           .explode("g").group_by("g").agg(pl.col("w").sum())
           .collect())
    return {r["g"]: r["w"] for r in g.iter_rows(named=True)}


def _empirical(sample, call_col: str) -> dict[str, float]:
    """P(gene) from one clonotype frame: normalized :func:`_vote_counts`."""
    counts = _vote_counts(sample, call_col)
    tot = sum(counts.values())
    if tot == 0:
        raise ValueError(f"no usable {call_col} values in the sample — cannot rescale usage")
    return {g: w / tot for g, w in counts.items()}


def _dataset_usage(sample, call_cols: list[str], aggregate: str) -> dict[str, dict[str, float]]:
    """Target usage ``{call_col: {gene: P}}`` from one frame, or a whole dataset streamed sample by sample.

    A single frame (``DataFrame`` or ``LazyFrame``) delegates to :func:`_empirical` per column. An
    **iterable of frames** is a dataset/cohort processed in ONE streaming pass — each sample is reduced
    to small per-gene vote counts and discarded, so the whole dataset is never held in memory. Pass a
    generator, or per-sample ``pl.scan_parquet(path)`` LazyFrames, to stream straight from disk. Both
    requested columns (V and J) are accumulated in the same pass, so a one-shot generator is safe.

    * ``"pool"`` — count every clonotype across the dataset once (usage weighted by each sample's
      clonotype richness); ``"mean"`` — average the per-sample usage (each sample equal, depth-robust).
    """
    if isinstance(sample, (pl.DataFrame, pl.LazyFrame)):
        return {c: _empirical(sample, c) for c in call_cols}
    if aggregate not in ("pool", "mean"):
        raise ValueError(f"aggregate must be 'pool' or 'mean', got {aggregate!r}")
    pooled: dict[str, dict[str, float]] = {c: defaultdict(float) for c in call_cols}
    summed: dict[str, dict[str, float]] = {c: defaultdict(float) for c in call_cols}
    n = 0
    for frame in sample:
        counts = {c: _vote_counts(frame, c) for c in call_cols}
        if not any(counts.values()):
            continue
        n += 1
        for c in call_cols:
            tot = sum(counts[c].values())
            for gene, w in counts[c].items():
                pooled[c][gene] += w
                if tot:
                    summed[c][gene] += w / tot
    if n == 0:
        raise ValueError(f"no usable {'/'.join(call_cols)} values in the dataset — cannot rescale usage")
    out: dict[str, dict[str, float]] = {}
    for c in call_cols:
        if aggregate == "pool":
            tot = sum(pooled[c].values())
            out[c] = {gene: w / tot for gene, w in pooled[c].items()} if tot else {}
        else:
            out[c] = {gene: w / n for gene, w in summed[c].items()}
    return out


def _reweight(table: pl.DataFrame, allele_col: str, target: dict[str, float],
              group_keys: list[str]) -> pl.DataFrame:
    """Set each gene's total mass to ``target[gene]``, preserving the within-gene allele split.

    The model is keyed by allele and the target is per gene, so a gene's new mass is divided
    among its alleles **in the model's existing proportions** — the sample cannot resolve alleles
    (see :func:`_empirical`), so the model's own split is the best available and is preserved
    rather than replaced by a uniform guess. A gene the model gives zero mass to is spread
    uniformly over its alleles, since there is no existing split to preserve.
    """
    t = table.with_columns(_gene(allele_col).alias("_g"))
    gene_tot = pl.col("p").sum().over([*group_keys, "_g"])
    n_alleles = pl.len().over([*group_keys, "_g"])
    share = pl.when(gene_tot > 0).then(pl.col("p") / gene_tot).otherwise(1.0 / n_alleles)
    t = t.with_columns(
        p=share * pl.col("_g").replace_strict(target, default=0.0, return_dtype=pl.Float64)
    )
    # Renormalize within each conditioning group: genes absent from the sample get 0, so the
    # remaining mass must be rescaled back to 1 (and a group the sample never uses stays 0).
    tot = pl.col("p").sum().over(group_keys) if group_keys else pl.col("p").sum()
    return t.with_columns(p=pl.when(tot > 0).then(pl.col("p") / tot).otherwise(0.0)).drop("_g")


def rescale_usage(model: Model, sample: pl.DataFrame | list[pl.DataFrame], *,
                  v: bool = True, j: bool = True, aggregate: str = "pool") -> Model:
    """Return a copy of ``model`` with V/J usage taken from ``sample``, junction model untouched.

    Args:
        model: Source model (e.g. ``load_bundled("TRB", "learned")``).
        sample: Either one clonotype frame (canonical schema), or a **list of frames** (a dataset /
            cohort) to rescale to the whole set's usage. ``v_call``/``j_call`` define the target usage;
            one vote per row, so pass unique clonotypes rather than reads if the frame is read-level.
            Pass the repertoire you are actually going to score: an **out-of-frame** sample is not
            interchangeable with a functional one, because a pseudogene's rearrangements are never
            productive and are therefore enriched out of frame (measured on this TRB data: TRBV23-1
            is 0.8% of functional reads but 31% of out-of-frame clonotypes).
        v: Rescale ``v_choice``.
        j: Rescale ``j_choice`` — ``P(J)`` on a VDJ locus, ``P(J|V)`` on a VJ locus, where the
            per-V conditional is set to the sample's overall J usage (the sample cannot support a
            reliable per-V J estimate at typical depths).
        aggregate: How to combine a **list** of frames (ignored for a single frame). ``"pool"``
            (default) counts every clonotype across the dataset once (usage weighted by each sample's
            richness); ``"mean"`` averages the per-sample usage (each sample equal, robust to depth).

    Returns:
        A new :class:`Model`; the input is not modified.

    Raises:
        ValueError: If the sample/dataset has no usable V/J calls, or ``aggregate`` is unknown.

    Example:
        >>> m = rescale_usage(load_bundled("TRB", "learned"), my_sample)
        >>> cohort = rescale_usage(load_bundled("TRB", "learned"), [s1, s2, s3], aggregate="mean")
        >>> native.pgen_aa(m, "CASSIRSSYEQYF", "TRBV19*01", "TRBJ2-7*01")
    """
    call_cols = ([V_CALL] if v else []) + ([J_CALL] if j else [])
    usage = _dataset_usage(sample, call_cols, aggregate)   # one streaming pass over the dataset
    tables = dict(model.tables)
    if v:
        tables["v_choice"] = _reweight(tables["v_choice"], "v_allele", usage[V_CALL], [])
    if j:
        # VJ loci key j_choice on (v_allele, j_allele) = P(J|V); VDJ on (j_allele) = P(J).
        keys = ["v_allele"] if "v_allele" in tables["j_choice"].columns else []
        tables["j_choice"] = _reweight(tables["j_choice"], "j_allele", usage[J_CALL], keys)
    return Model(manifest=model.manifest, tables=tables, genomic=model.genomic)
