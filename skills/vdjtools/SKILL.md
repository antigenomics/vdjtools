# SKILL: vdjtools (v2)

Public API surface of **vdjtools v2** — a Python + C++ (pybind11 `_core`) rewrite for TCR/BCR
repertoire analysis on the **AIRR schema + polars**. Check here for an existing function before
writing new code. Keep this file current when a subpackage's public API changes.

## Canonical data model

Every reader emits and every analysis consumes one flat **clonotype frame** (`vdjtools.io.schema`):
`v_call, d_call, j_call, c_call, junction_aa, junction_nt, duplicate_count, frequency` (+ derived
`locus`). Columns are **AIRR junction** (conserved Cys104…Phe/Trp118 anchors *included*) — names are
`junction_nt` / `junction_aa` (v2.2.0 rename from `cdr3_nt`/`cdr3_aa`; readers still accept the old
names, strict-AIRR `junction`, IMGT `cdr3`, and native `cdr3nt`/`cdr3aa` as input aliases).
Minimal OO — free functions returning `pl.DataFrame`.

## Build / test / run

```bash
conda env create -f environment.yml && conda activate vdjtools   # or reuse .venv
pip install -e ".[dev,test]"                                     # builds the _core C++ ext
pytest tests/python -q -m "not slow"                             # fast suite
sphinx-build -W --keep-going -b html docs docs/_build/html       # docs (zero-warning gate)
```
Iterating on C++: `cmake --build build/<wheel_tag>` then copy `_core.*.so` into the venv's
`site-packages/vdjtools/`, or re-run `pip install -e .`.

## API by subpackage

### `vdjtools.io` — IO, schema, format converters
- **Readers**: `read(path, fmt="auto")` (sniffs + dispatches), `read_airr`, `read_vdjtools`,
  `read_parquet`; **converters** `read_mixcr` (v1/2+v3/4, incl. C-gene/BCR isotype), `read_migec`,
  `read_immunoseq` (Adaptive v1/v2), `read_imgt` (IMGT/HighV-QUEST), `read_vidjil` (JSON),
  `read_rtcr`, `read_trust4` (`*_report.tsv`), `read_arda` (arda AIRR output, delegates to
  `read_airr`) (`vdjtools.io.convert`; ported from the legacy Groovy parsers + tool docs, incl.
  Adaptive→IMGT gene conversion).
- **Cohorts**: `read_metadata`, `read_samples`, `iter_samples` (streaming), `sniff_format`;
  `ingest_cohort` / `scan_cohort` (hive-partitioned Parquet, lazy).
- **Schema**: `SCHEMA, COLUMNS`, constants `V_CALL D_CALL J_CALL C_CALL JUNCTION_AA JUNCTION_NT COUNT
  FREQ LOCUS`; helpers `normalize, add_locus, locus_of, recompute_frequency`.

### `vdjtools.model` — native V(D)J recombination engine (supersedes OLGA + IGoR)
- **Load**: `load_bundled(locus, source="olga"|"learned")`, `list_bundled`; `from_olga`, `load_model`,
  `save_model`. `Model`, `Manifest`, `Event`, `EventKind`.
- **Pgen (native)** `vdjtools.model.native`: `pgen_nt`, `pgen_aa(m, aa, v=None, j=None, mismatches=0)`
  (0=exact, 1=Hamming-1 ball; v/j=None marginalises), **`pgen_aa_batch(m, seqs, v=, j=, mismatches=,
  threads=)`** (thread-parallel across sequences, bitwise-identical to serial, ~11× on 16 cores).
  Pure-Python reference impls in `vdjtools.model.pgen`.
- **Generate**: `vdjtools.model.generate.generate(model, n, seed=, productive_only=)` → `pl.DataFrame`.
- **Infer (EM)**: `vdjtools.model.infer.infer` / `infer_native(template, seqs, masks=, dd_allowed=,
  nd_prior=, single_d=)`.
- **Germline (arda = single source of truth)** `vdjtools.model.reference`: `load_germline(locus,
  organism)` (CDR3-region + anchor), **`load_full_vj_germline(organism)`** and
  **`arda_full_germline(locus, organism)`** (full-length V/J germline + stitch anchor, from arda
  scaffolds), `reconcile_olga`, `cut_segment`, `translate`, `reverse_complement`.
- **Stitch**: `stitch_contig(model, v, j, cdr3_nt)`, `stitch_frame`.
- **Diagnostics** `vdjtools.model.analyze`: Bayes-net→graphviz DOT, entropy / mutual-information tables.
- Tandem-D (D-D) supported throughout (`vdjtools.model.dd`).

### `vdjtools.stats` — diversity, spectratype, usage
- `diversity_stats` (all indices), individual: `observed_richness, chao1, chao_e, efron_thisted,
  shannon_wiener, normalized_shannon_wiener, inverse_simpson, d50`.
- Rarefaction: `rarefaction`, `inext` (Hill q=0/1/2, size+coverage), `inext_batch`,
  `rarefaction_batch`, `inext_coverage`, `asymptotic_diversity`, `coverage`, `sample_coverage`,
  `estimate_d`.
- `segment_usage`, `vj_usage`, `spectratype`, `vj_spectratype`.

### `vdjtools.features` — CDR features
`physchem_profile` (region × property), `kmer_profile`, `v_kmer_c_profile`, `load_property_table`,
`DEFAULT_PROPERTIES`.

### `vdjtools.overlap` — overlap + TCRnet (delegates to vdjmatch/seqtree)
`overlap_metrics`, `overlap_pair`, `DEFAULT_KEY`; `fuzzy_overlap`, `fuzzy_overlap_metrics`;
`similarity_overlap`, `similarity_matrix`, `SimilarityMatrices` (TINA / Leinster-Cobbold);
`tcrnet`; `pairwise_distances`, `cluster_samples`; `track_clonotypes`.

### `vdjtools.preprocess`
`downsample`, `select_top`; `filter_functional`, `filter_frequency`, `filter_segment`,
`filter_by_sample`; `correct` (freq error-correction), `decontaminate`; `pool_samples`,
`join_samples`, `resolve_key`; `correct_vj_usage` (VJ batch-effect: `transform="location"` ComBat
default or `"sigmoid"` = Vlasova 2026 z-score + grand-mean-preserving sigmoid), `apply_vj_correction`
(rescale + roulette-wheel resample the clonotype table to the corrected usage).

### `vdjtools.biomarker`
Incidence contingency testing across a cohort (Emerson 2017 / Howie 2015 / De Witt 2018 / Vlasova 2026).
- `association(cohort, design, *, test=, level_col=, stratum_col=, key=, match=, min_incidence[_frac]=,
  candidates=, alternative=)` — feature-vs-condition; `test` ∈ {`fisher`,`chi2`,`bayes_logodds`,
  `bayes_bf`,`permutation`} (str or list → long output w/ `test` col); category via `level_col` (one-vs-rest),
  paired via `stratum_col` (Cochran–Mantel–Haenszel). Match scope = `key` (`(junction_aa,)`/`+v`/`+v+j`) × `match` (`exact`/`1mm`).
- `cooccurrence(cohort, *, chain_a=, chain_b=, test=, min_incidence[_frac]=, min_cooccurrence=, evalue=, depth_strata=10)` —
  **depth-conditioned by default** (CMH over repertoire-depth strata): a deep repertoire carries more of
  everything, so a pooled test is badly miscalibrated (measured FPR 0.46 on independent pairs at the
  incidence regime `max_features` selects). `depth_strata=0` restores the pooled test. Adds `or_mh`/`chi2`.
  feature-vs-feature θ=n·n_AB/(n_A·n_B) + Fisher/χ² + FDR; α-β pairing (chain_a≠chain_b) or same-chain (chain_b=None).
- `condition` builders: `binary`, `categorical`, `hla_alleles`, `zygosity`, `stratified` → design frame (`_pos`/`_level`/`_stratum`).
- `select_candidates` (public features over incidence count/fraction), `stats` (vectorised 2×2 kernels),
  `fisher_association` (Emerson Fisher shortcut, legacy schema), `metaclonotypes` (1mm grouping).

### `vdjtools.sc` — single-cell (AIRR Cell / 10x)
`read_10x`, `read_airr_cell`, `write_airr_cell`; `resolve_chains`, `pair_chains`,
`chain_multiplicity`, `flag_mispairing`; **`paired_pgen(paired, source=, condition_vj=)`**
(`Pgen(α)·Pgen(β)` via the native model); `cluster_eval` (+ `purity`, `homogeneity`, `parsimony`,
`q_measure`, …); `to_anndata`.

### `vdjtools.cli`
The `vdjtools` typer app: `models`, `generate`, `pgen`, `diversity`, `overlap`, `segment-usage`,
`spectratype`. Inputs auto-detected; TSV to `-o`/stdout.

## Conventions
- **arda germline is the single source of germline truth** — resolve V/D/J germline + CDR3 anchors
  by allele name via `model.reference`. Never mix germline sources within a model.
- Delegate: overlap/TCRnet → **vdjmatch**; annotation/markup/germline → **arda**; search/e-value →
  **seqtree**. Only `arda`/`vdjmatch`/`seqtree` are dependencies; everything else under `~/vcs/code/`
  is read-only oracle/reference (never modify).
- Native C++ only for the hot loops (Pgen DP, generation sampler, EM E-step) via the single `_core`
  ext; everything else is polars.
