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

## 2. Adaptive/immunoSEQ → IMGT gene names are wrong for 100 of 161 tokens — ✅ IMPLEMENTED

**Status (2026-08-09): done, shipped in 3.2.0.** `io/convert.py::_adaptive_to_imgt` looks the raw
Adaptive token up in the shipped table and falls back to the legacy rewrite on a miss. Tests:
`tests/python/test_convert.py::test_adaptive_to_imgt` (the design note's self-check plus the
off-table fallback) and `::test_read_immunoseq_uses_the_validated_gene_map` (V/D/J end-to-end).
The optional `strict=True` keyword was **not** implemented — nothing needs it yet.
Original notes below for provenance. Design note: **`appendix/adaptive_imgt_map.md`**.
Reference table: **`python/vdjtools/resources/adaptive_imgt_map.tsv`** (provenance in `SOURCES.md`).
Builder: **`appendix/build_adaptive_imgt_map.py`**.

**What's missing.** `io/convert.py::_adaptive_to_imgt` normalises Adaptive names with a global
`re.sub(r"0([1-9])", r"\1", …)`, which always re-emits the trailing `-01` as an IMGT *subgroup*. On
the 161 distinct Adaptive tokens of the IMMREP25 release + the pairSEQ mock cohort (44,000 gene
calls) it disagrees with the CDR-validated map on **100 tokens / 22,058 calls**, and **every one of
those 100 outputs is a gene name absent from the IMGT human reference** — `TCRAJ39-01 → "TRAJ39-1"`,
`TCRBV09-01 → "TRBV9-1"`, `TCRAV22-01 → "TRAV22-1"`, `TCRBD01-01 → "TRBD1-1"`,
`TCRBV03-01/03-02 → "TRBV3-1/3-2"`, `TCRBV20-X → "TRBV20-X"`, `TCRAV38-02 → "TRAV38-2"`
(IMGT: `TRAV38-2/DV8`).

**Why it matters.** Any consumer that resolves gene names against a germline reference silently
loses the rows. tcrdist3 dropped **100 %** of both affected cohorts
(`input rows=1000; kept=0 (dropped 1000 unknown-gene)`); with the table, 1000/1000 and 993/1000.
The legacy Groovy `CommonUtil.extractVDJImmunoSeq` has the same defect, so this is a v1 bug inherited
by v2, not a porting regression.

**Why a regex can't fix it.** Whether the trailing group is a subgroup or an allele is a per-family
fact: `TRAV1` has subgroups so `TCRAV01-01 = TRAV1-1`; `TRAV22` has none so `TCRAV22-01 = TRAV22`;
no human `TRAJ` has one (61 genes, 68 alleles, zero `-N`). The two readings are textually identical.

**Proposed API (io.convert).** Keep the signature; look the token up in the shipped table first and
fall back to the current regex on a miss (so unknown tokens behave exactly as today). Optional
`strict=True` returns `None` instead of falling back. Full patch in the design note; blast radius is
`_adaptive_call` → `read_immunoseq` (V, D **and** J calls) → `io.batch` `fmt="immunoseq"`; no other
reader touches it, and the existing `immunoseq`/`immunoseqv2` test oracles (`TRBV29-1`, `TRBJ2-6`)
are unchanged by the table.
