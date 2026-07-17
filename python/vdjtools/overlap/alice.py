"""ALICE neighbourhood enrichment — the *generative* null, against TCRnet's *empirical* one.

ALICE (Pogorelyy et al., PLoS Biol 2019) asks the same question as :mod:`vdjtools.overlap.tcrnet`
— does this clonotype have more close CDR3 neighbours than chance? — but answers it from a V(D)J
**generation model** rather than from a control repertoire. For each amino-acid CDR3 ``σ`` within
one V–J class::

    λ_σ  =  n · Σ_{σ' : Hamming(σ,σ') ≤ 1}  Q · Pgen(σ')
    p    =  P(Poisson(λ_σ) ≥ d(σ))

where ``n`` is the number of unique **nucleotide** clonotypes in that V–J class, ``d(σ)`` counts
the nucleotide clonotypes whose CDR3aa lies within one mismatch of ``σ`` (σ included — different
nucleotide variants of one amino-acid sequence are genuine neighbours), and ``Q`` rescales for
thymic selection, which removes a fraction ``1 − 1/Q`` of generated sequences.

The two nulls are complements, not rivals, and that is why both are here: ALICE's generative null
controls for the intrinsic biases of V(D)J recombination but knows nothing about selection or
about which clonotypes are *already* common in people; TCRnet's control repertoire absorbs thymic
selection and endemic-pathogen expansions but needs a large, HLA-matched cohort to do it.

Both share a known blind spot worth stating plainly: neighbourhood enrichment **cannot see a
monoclonal expansion**. A single hyperexpanded clone has no near neighbours by definition, so it
scores as unremarkable no matter how dominant it is. That is precisely why this module and
:mod:`vdjtools.dynamics` are complementary — enrichment measures *breadth*, the paired test
measures *magnitude*.

The Hamming-1 ball sum is exact and already native: ``pgen_aa_batch(..., mismatches=1)`` computes
``Σ_k Pgen(a_{k→*}) − (L−1)Pgen(a)``, which is the **closed** ball with the centre counted exactly
once — verified here against brute-force enumeration of all 19L neighbours to 2e-16.
"""
from __future__ import annotations

import numpy as np
import polars as pl
from scipy.stats import poisson

from ..biomarker.stats import fdr_bh
from ..io.schema import COUNT, J_CALL, JUNCTION_AA, JUNCTION_NT, LOCUS, V_CALL, add_locus

#: ALICE's default thymic-selection factor — the paper's average over V–J combinations.
DEFAULT_Q = 9.41

_COLS = [JUNCTION_AA, V_CALL, J_CALL, COUNT, "n_neighbors", "n_group", "pgen_ball",
         "E", "p_enrichment", "q_value", LOCUS]


def _alleles_for(idx_of: dict[str, int], name: str) -> list[str]:
    """The model alleles a V/J call denotes: itself if it is one, else the gene's alleles.

    Real repertoires carry gene-level calls (``TRBV19``); the model is keyed by allele. Passing a
    gene straight to the native Pgen raises (it used to silently marginalise over *every* allele —
    measured 183x too high). Summing a gene's alleles is not an approximation: ``Σ_a Pgen(σ, V=a)``
    over a gene's alleles *is* the gene marginal, verified to ratio 1.0000000000 against the
    all-allele sum.
    """
    if name in idx_of:
        return [name]
    alleles = sorted(a for a in idx_of if a.split("*")[0] == name)
    if not alleles:
        raise KeyError(f"V/J call {name!r} is not in the model and is not a gene of any of its "
                       f"{len(idx_of)} alleles")
    return alleles


def _ball_pgen(model, seqs: list[str], v: str, j: str, threads: int) -> np.ndarray:
    """Σ Pgen over the closed Hamming-1 ball of each ``seq``, at gene *or* allele resolution."""
    from ..model.native import pack, pgen_aa_batch

    _, vi, ji = pack(model)
    vs, js = _alleles_for(vi, v), _alleles_for(ji, j)
    total = np.zeros(len(seqs), dtype=float)
    for va in vs:                       # a gene-level call sums over its alleles (exact)
        for ja in js:
            total += np.asarray(pgen_aa_batch(model, seqs, [va] * len(seqs), [ja] * len(seqs),
                                              mismatches=1, threads=threads))
    return total


def _degrees(unique_aa: list[str], expanded_aa: list[str], scope: str, threads: int) -> np.ndarray:
    """Nucleotide-clonotype degree of each unique CDR3aa, self included.

    ``expanded_aa`` holds one entry per *nucleotide* clonotype, so a hit count over it is already
    in ALICE's units. ``vdjmatch.cluster.overlap`` returns self-pairs, so the ball is closed —
    matching ``λ``, which sums Pgen over the closed ball.
    """
    from vdjmatch.cluster import overlap

    pairs = overlap(unique_aa, expanded_aa, scope=scope, threads=threads)
    counts = (pairs.group_by("a_idx").len()
                   .rename({"len": "d"}).sort("a_idx"))
    out = np.zeros(len(unique_aa), dtype=np.int64)
    out[counts["a_idx"].to_numpy()] = counts["d"].to_numpy()
    return out


def alice(sample: pl.DataFrame, model=None, *, locus: str | None = None,
          source: str = "olga", scope: str = "1,0,0,1", selection_q: float = DEFAULT_Q,
          min_degree: int = 3, min_count: int = 2, threads: int = 0) -> pl.DataFrame:
    """Per-clonotype neighbourhood enrichment against a V(D)J generation model.

    Args:
        sample: Clonotype frame (canonical schema). ``junction_nt`` is used when present so
            degrees are in nucleotide-clonotype units, as ALICE defines them; without it each
            amino-acid sequence counts once.
        model: A :class:`~vdjtools.model.Model`. ``None`` loads the bundled model for ``locus``.
        locus: Locus to score (e.g. ``"TRB"``). ``None`` infers it from ``v_call``; a sample
            spanning several loci is scored per locus.
        source: Bundled model source. **Pinned to ``"olga"`` by default and you should leave it
            there**: the ``"learned"`` models are EM-fit on ~2k clonotypes with no gene-usage
            pseudocount, so 68 of 89 bundled TRB V alleles have ``P(V) = 0`` (vs 8 for OLGA) —
            those clonotypes get ``λ = 0`` and come back infinitely significant. Their ball-Pgen
            scale is also ~16x below OLGA's, which ``Q`` is calibrated against.
        scope: vdjmatch edit scope defining the neighbourhood (default one substitution).
        selection_q: Thymic-selection factor ``Q``. The default is the paper's average over V–J
            combinations; it is calibrated on OLGA's **TRB** Pgen scale and there is no evidence
            it transfers to the other loci, so the benchmark titrates it per locus.
        min_degree: Only clonotypes with at least this many neighbours (self included) are
            tested; the rest are dropped. ALICE tests ``d(σ) > 2``.
        min_count: Clonotypes below this count do not participate at all. Low-count variants of
            an abundant clonotype are usually sequencing error, and counting them as neighbours
            inflates every degree in their neighbourhood.
        threads: Worker threads (``0`` = all cores).

    Returns:
        One row per tested clonotype: the key columns, ``duplicate_count``, ``n_neighbors``
        (``d(σ)``, self included), ``n_group`` (``n``, the V–J class's nucleotide-clonotype
        count), ``pgen_ball`` (``Σ Pgen`` over the closed ball), ``E`` (``λ_σ``),
        ``p_enrichment`` (Poisson tail), ``q_value`` (BH FDR over the clonotypes tested in this
        call) and ``locus`` — sorted by ascending ``p_enrichment``. **No threshold is applied**:
        ALICE's own papers use BH < 0.001 while the TCRnet framework paper used 0.05 for the same
        family of tests, and that was never reconciled — so the caller chooses.

    Raises:
        ValueError: If no locus can be resolved, or the sample has no usable clonotypes.
        KeyError: If a V/J call is neither a model allele nor a gene of one.
    """
    from ..model.bundled import load_bundled

    df = sample.filter(pl.col(COUNT) >= min_count)
    if df.is_empty():
        raise ValueError(f"no clonotypes with duplicate_count >= {min_count}")
    df = add_locus(df).filter(pl.col(LOCUS).is_not_null())
    if locus is not None:
        df = df.filter(pl.col(LOCUS) == locus)
    if df.is_empty():
        raise ValueError("cannot infer locus from v_call; pass locus=")

    parts = []
    for loc in sorted(df[LOCUS].unique().to_list()):
        mdl = model if model is not None else load_bundled(loc, source)
        parts.append(_score_locus(df.filter(pl.col(LOCUS) == loc), mdl, loc, scope,
                                  selection_q, min_degree, threads))
    out = pl.concat([p for p in parts if not p.is_empty()]) if any(
        not p.is_empty() for p in parts) else parts[0]
    if out.is_empty():
        return out.select(_COLS)
    q = fdr_bh(out["p_enrichment"].to_numpy())
    return out.with_columns(pl.Series("q_value", q)).select(_COLS).sort("p_enrichment")


def _score_locus(df: pl.DataFrame, model, loc: str, scope: str, selection_q: float,
                 min_degree: int, threads: int) -> pl.DataFrame:
    """Score one locus, V–J class by V–J class (ALICE's null is defined per V–J)."""
    has_nt = JUNCTION_NT in df.columns and df[JUNCTION_NT].null_count() < df.height
    nt_key = [JUNCTION_NT, JUNCTION_AA, V_CALL, J_CALL] if has_nt else [JUNCTION_AA, V_CALL, J_CALL]
    nt = df.group_by(nt_key).agg(pl.col(COUNT).sum())      # unique nucleotide clonotypes

    rows = []
    for (v, j), grp in nt.group_by([V_CALL, J_CALL]):
        if v is None or j is None:
            continue
        expanded = grp[JUNCTION_AA].to_list()              # one entry per nt clonotype
        agg = (grp.group_by(JUNCTION_AA).agg(pl.col(COUNT).sum())
                  .sort(JUNCTION_AA))
        uniq = agg[JUNCTION_AA].to_list()
        deg = _degrees(uniq, expanded, scope, threads)
        keep = deg >= min_degree                           # ALICE tests d(sigma) > 2
        if not keep.any():
            continue
        seqs = [s for s, k in zip(uniq, keep) if k]
        ball = _ball_pgen(model, seqs, v, j, threads)
        lam = len(expanded) * selection_q * ball
        rows.append(agg.filter(pl.Series(keep)).with_columns(
            pl.Series("n_neighbors", deg[keep]),
            pl.lit(len(expanded), dtype=pl.Int64).alias("n_group"),
            pl.Series("pgen_ball", ball),
            pl.Series("E", lam),
            pl.Series("p_enrichment", poisson.sf(deg[keep] - 1, lam)),
            pl.lit(v).alias(V_CALL), pl.lit(j).alias(J_CALL),
        ))
    if not rows:
        return pl.DataFrame(schema={c: pl.Utf8 for c in _COLS})
    return pl.concat(rows).with_columns(pl.lit(loc, dtype=pl.Utf8).alias(LOCUS))
