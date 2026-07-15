# Suggested vdjtools edits

Repo-local notes on gaps found while using vdjtools "by hand" for analyses. Each entry: what's
missing, why, reference implementation, and a proposed API. Promote to GitHub issues as needed.

## 1. V-J usage batch correction — incomplete vs Vlasova et al. 2026 — ✅ IMPLEMENTED

**Status (2026-07-16): done.** `preprocess/batch.py` now provides both pieces:
- `correct_vj_usage(..., transform="sigmoid", z_cap=6.0)` — the σ-standardised z-score with the
  **grand-mean-preserving sigmoid** `P_final = 2·P_avg/(1+exp(−Z))` (owner-confirmed formula;
  legacy mirpy's own `compute_batch_corrected_gene_usage` uses `p·exp(Z)` instead — this follows
  the paper's Methods, not the legacy code). `transform="location"` is unchanged (default).
- `apply_vj_correction(sample_df, corrected_usage, *, scope, weighted, resample=True, seed)` —
  rescales each clonotype by `P_final(G)/P(G)` and roulette-wheel resamples to a new integer-count
  table (multinomial at the original read depth), or `resample=False` for deterministic expected
  counts. Port of legacy mirpy v2 `resample_to_gene_usage`.

Tests: `tests/python/test_preprocess_batch.py` (divergence removal, grand-mean preservation,
value-pin, total-read preservation, determinism). Original notes retained below for provenance.



**Context.** Reproducing Vlasova, Nekrasova, Komkov, … Britanova, Shugay, *Inference of SARS-CoV-2
exposure biomarkers using large-scale T-cell repertoire profiling*, **Genome Medicine 2026;18:20**
(DOI 10.1186/s13073-025-01589-4). The paper's batch-effect correction operates on **clonotype
tables**, not just gene-usage profiles.

**Current state.** `preprocess/batch.py::correct_vj_usage` is a **location-only ComBat** adjustment:
per-`(locus, gene, batch)` winsorized mean of `log p` → subtract batch mean, add grand mean →
`p_corrected` (renormalised gene usage). It **stops at corrected gene usage** and deliberately omits
the scale (σ) term.

**Paper's full method (Methods, "Batch-effect correction and data normalization"):**
1. Per gene per batch, model usage `P(gene,sample) ~ LogNormal(μ, σ | batch)`; compute
   **Z-score** `Z = (log P − μ)/σ` — **uses the σ/scale term** (validated normality via Shapiro–Wilk).
2. Map Z back to `[0,1]` with a **sigmoid that preserves the grand-mean usage** `P_avg(gene)`:
   `P_final(gene,sample) = 2·P_avg(gene) / (1 + exp(−Z(gene,sample)))`.
3. **Rescale clonotype frequencies** `f_i' = f_i · P_final(G,S)/P(G,S)` (G = clonotype's V or J gene),
   then **resample** the clonotype composition (roulette-wheel selection from `U[0,1]` scaled to read
   count) → a new integer-count clonotype table with corrected V/J usage.

**Gap.** vdjtools has (1) location only — no σ term, (2) no sigmoid-preserving-`P_avg` map, and (3) **no
step that applies the corrected usage back to clonotype frequencies / resamples the clonotype table.**
The full method is in **legacy mirpy v2**: `mir.basic.gene_usage` (`compute_batch_corrected_gene_usage`)
+ `mir.common.sampling.resample_to_gene_usage`.

**Proposed API (preprocess):**
- extend `correct_vj_usage(..., scale=False, transform="location"|"sigmoid")` — add the σ Z-score +
  `P_avg`-preserving sigmoid path to match the paper;
- add `apply_vj_correction(sample_df, corrected_usage, *, resample=True, seed=0) -> pl.DataFrame` —
  rescale `duplicate_count`/`freq` by `p_corrected/p` per clonotype and (optionally) roulette-wheel
  resample to integer counts. Port `resample_to_gene_usage` from legacy mirpy v2.

Until this lands, the reproduction implements the rescale+resample step by hand in the analysis script.
