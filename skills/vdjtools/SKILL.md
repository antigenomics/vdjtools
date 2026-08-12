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
uv venv && source .venv/bin/activate && uv pip install -e ".[dev,test]"   # builds _core (default)
# or `bash setup.sh` (uv-first, portable bash/zsh); conda env.yml only for mmseqs2 + slow arda tests
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
  `read_airr`) (`vdjtools.io.convert`; ported from the legacy Groovy parsers + tool docs).
  Adaptive→IMGT gene names come from the shipped CDR-validated `resources/adaptive_imgt_map.tsv`
  (subgroup-vs-allele is a per-family fact: `TCRAV01-01`→`TRAV1-1` but `TCRAV22-01`→`TRAV22`;
  family calls, slash ties and `/DVn` co-locus names resolve too), legacy zero-strip off-table.
- **Cohorts**: `read_metadata`, `read_samples`, `iter_samples` (streaming), **`map_samples(fn,
  items, *, workers=)`** (thread-parallel per-sample reduce, `O(workers)`-sample RAM, input-order),
  `sniff_format`; `ingest_cohort` / `scan_cohort` (hive-partitioned Parquet, lazy).
- **Schema**: `SCHEMA, COLUMNS`, constants `V_CALL D_CALL J_CALL C_CALL JUNCTION_AA JUNCTION_NT COUNT
  FREQ LOCUS`; helpers `normalize, add_locus, locus_of, recompute_frequency`.

### `vdjtools.model` — native V(D)J recombination engine (supersedes OLGA + IGoR)
- **Load**: `load_bundled(locus, source="olga"|"learned"|"arda")`, `list_bundled`; `from_olga`,
  `load_model(path, validate=False)`, `save_model(m, path, fmt="parquet"|"tsv"|"csv")`.
  `Model` (`.manifest .tables .genomic .training`), `Manifest`, `Event`, `EventKind`.
- **Build from any reference**: **`from_germline(germline_df, locus=, organism=, ins_max=, strict=)`**
  — the custom-library entry point; `from_arda(locus, organism)` is a thin wrapper over it.
  `reference.read_germline_fasta(v, j, d=None, anchors=)` builds the frame from your own FASTA;
  `reference.validate_germline(df)` audits it (anchor frame, IUPAC, duplicate/gene-level names).
  Required germline columns `allele segment sequence`; optional `gene functional cdr3_anchor
  full_germline`. A D allele present ⇒ VDJ, absent ⇒ VJ.
- **Tables in/out**: **`marginals_frame(m)`** → one long frame of every marginal;
  **`set_marginals(m, frame)`** → back to a `Model` (a hand-edited TSV is a first-class input).
- **Pgen (native)** `vdjtools.model.native`: `pgen_nt`, `pgen_aa(m, aa, v=None, j=None, mismatches=0)`
  (0=exact, 1=Hamming-1 ball; v/j=None marginalises), **`pgen_aa_batch(m, seqs, v=, j=, mismatches=,
  threads=)`** (thread-parallel across sequences, bitwise-identical to serial, ~11× on 16 cores).
  Pure-Python reference impls in `vdjtools.model.pgen`.
- **Generate**: `vdjtools.model.generate.generate(model, n, seed=, productive_only=)` → `pl.DataFrame`.
  NOTE: `seed=` is reproducible **across processes** only from **3.3.0**: `collapse_alleles` used
  unordered polars `group_by`, so the collapsed table's row order varied per process and the same
  seed drew a different allele. Expectations recorded from `generate()` before 3.3.0 are stale.
- **Infer (EM)**: `vdjtools.model.infer.infer` / `infer_native(template, seqs, masks=, dd_allowed=,
  nd_prior=, single_d=, init="align"|"uniform"|"template")`; **`infer_frame(template_or_locus,
  clones_df)`** takes a clonotype frame and builds the V/J masks for you. `init="template"` is the
  warm start = fine-tuning. **Long fits**: `progress=print_progress()` reports per-iteration loglik
  + rel_change; `checkpoint=DIR`/`checkpoint_every=N` saves each iteration and `resume(DIR, seqs)`
  continues (exact — a resumed run matches an uninterrupted one, tables and all). CLI:
  `model learn -v --checkpoint DIR --resume DIR`. **Training log**: every fit appends to
  `model.training["runs"]`, saved
  as a `training.json` sidecar; read it with **`training_frame(model)`** (`run iter loglik
  n_scoreable rel_change`). Bundled models have `training is None`.
- **Extend the allele library**: **`extend_alleles(model, germline_df, weight=1.0)`** — seeds new
  alleles from a gene-mate (or the germline-nearest allele for a brand-new gene, at a floor mass),
  clipping copied deletion rows to the new germline. **Preserves each pre-existing gene's total
  usage** — alleles of a gene are alternatives, not extra genes. Seeds only; follow with
  `infer_native(..., init="template")`. `augment_from_oracle(learned, oracle)` fills gaps from
  another *model* instead.
- **Collapse to gene level** (`load_bundled(..., collapse=True)`, the default): representative
  germline = **longest CDR3-region germline first, usage second** (a truncated IMGT allele must
  never define the gene's trim range — see IGKV3-20), then the deletion conditionals are projected
  onto that germline's reachable support and renormalized. `collapse=False` for exact-OLGA fidelity.
- **Check**: **`check_model(m, germline="auto"|"none"|df, raise_on=None)`** → tidy issue frame
  `severity check event segment allele detail value`. Catches unnormalized or out-of-range
  probabilities, alleles missing from the germline (or vice versa), functional genes stuck at P=0,
  mass on unreachable deletions (as a per-allele **fraction**), incomplete dinucleotide tables, and
  VDJ/VJ event-set mismatches. `severity == "error"` means the model is broken.
- **Score** `vdjtools.model.score`: `pgen_frame(m, seqs_or_frame, kind="auto", use_calls=)`,
  **`model_fit(m, seqs, weights=)`** (loglik/AIC/**BIC**; uses **nt** Pgen so the likelihood is
  normalized — aa is a relative score only; Pgen 0 is counted, never `-inf`),
  `free_params(m, by_event=)` (support-based, drops undefined and unreachable conditionals),
  `compare_pgen(a, b, seqs)` + `pgen_summary(cmp)` (KS, Spearman, and the headline
  `only_a_scoreable`/`only_b_scoreable`), `pgen_spectrum(m, n=, bins=)`,
  **`diversity(m, n=, seed=)`** (scenario entropy, Monte-Carlo sequence entropy, Hill q=1 `2^H`
  and q=2 `1/E[Pgen]`; TRB ≈ 52 bits scenario / 45 bits sequence / ~3e13).
- **Germline (arda = single source of truth)** `vdjtools.model.reference`: `load_germline(locus,
  organism)` (CDR3-region + anchor), **`load_full_vj_germline(organism)`** and
  **`arda_full_germline(locus, organism)`** (full-length V/J germline + stitch anchor, from arda
  scaffolds), `reconcile_olga`, `cut_segment`, `translate`, `reverse_complement`.
- **Scenario (argmax)** `vdjtools.model.viterbi`: **`best_scenario(model, cdr3_nt, v=, j=)`** →
  `Scenario(cdr3_nt, v_call, j_call, v_end, j_start, d_call, d_start, d_end, scenario_p, ...)` — the
  single most likely recombination for a KNOWN nt CDR3, i.e. the V/D/J boundary markup. Max-product
  over the same loops `pgen_nt` sums, so the D obeys `p_d_given_j` (TRBD2×TRBJ1 = 0). Coordinates
  are 0-based half-open in CDR3-nt space.
- **aa → nt** **`infer_nt(model, cdr3_aa, v=, j=, n_best=8)`** → the same `Scenario`, with `pgen`
  the exact `pgen_nt` of the returned sequence — the VDJdb case, where a record has `(V, J, CDR3aa)`
  and no nucleotides. Two stages: the argmax over every scenario `pgen_aa` sums (germline pinned,
  free N-region positions scored by the VD/DJ/VJ dinucleotide model), then a `pgen_nt` re-score of
  the survivors — stage 1 maximises the *joint*, the re-score answers about the *marginal*.
  Reproduces `infer_nt_bruteforce` (the exact but exponential ORACLE, tests only) 25/25 TRG and
  19/19 TRA; the one-scenario "best codon per residue" shortcut manages 9/25 and 4/19, because a
  trim chosen before the codons pins a codon the true optimum would have trimmed away.
  **`native.best_aa_scenarios(model, aa, v=, j=, k=8)`** is stage 1 alone (top-k scenarios as
  `(w, v, len_v, j, len_j, d, idx5, idx3, pos)`): the same Pi_L*Pi_R transfer matrix as `pgen_aa`
  with `max` for the sums. **2.5 ms/TRB, 0.5 ms/TRA — all of VDJdb in ~3 min.**
  `v=`/`j=` take one allele, a list or comma-separated string (ambiguous `v_call`), or `None`
  (marginalize — barely slower, 0.26 vs 0.23 ms). NOTE: pass a `Model`, not a `prepare()`-d one:
  a prepared model selects the pure-Python reference search, ~600x slower on TRB (tests cross-check
  the two). WARNING: Do not pin germline flanks to shrink the search — it drops trimmed-germline
  sequences and the true max can be one. Tandem-D is not enumerated in stage 1 (it cannot add
  candidates, only reorder them) but is fully counted by the stage-2 `pgen_nt`.
- **Stitch**: `stitch_contig(model, v, j, cdr3_nt)`, `stitch_frame`.
- **Usage re-weighting**: `rescale_usage(model, sample_or_list, v=, j=, aggregate="pool"|"mean")` —
  V/J usage is protocol-dependent (5'RACE vs DNA-multiplex), the junction model is not.
- **Diagnostics** `vdjtools.model.analyze`: `entropy_table`, `mutual_information`,
  **`total_entropy`** (per-event contribution to the scenario entropy; the dinucleotide row is
  `E[len] × H_step`), `gene_marginal`, `bayes_net_dot` / `render_bayes_net` / **`render_dot`**.
- **Compare two models** `vdjtools.model.analyze`: **`compare_models(a, b, by="allele"|"gene")`**
  → per-event `status n_groups support_* tv tv_max jsd_bits` (union-aligned with zero fill,
  parent-marginal-weighted; JSD is primary because it stays finite on disjoint support; `by="gene"`
  bridges different germline namespaces; `status="schema_differs"` when an event is factorized
  differently). **`compare_usage(a, b, seg)`**, **`compare_net_dot(a, b)`** (bnlearn
  `compare_networks` style), `compare_entropy`.
- **Corpus build** `vdjtools.model.data`: `load_prepared(group, chain, label)` (arda-mapped TRA/TRB
  examples shipped in `tests/python/fixtures/model_reads/` as gzipped FASTA with the V/J/D calls in
  the header — offline, no arda; source tree only, not in the wheel), `write_prepared`,
  `build_model(chain, ...)`, **`build_all(chains, groups=, workers=)`** — the full FASTQ → arda-map
  → EM pipeline, parallel across chains. Ambiguous junction bases → `A` by default in both training
  entry points (`infer.sanitize_junctions`; `ambiguous=None` drops instead).
  NOTE: arda ≥2.19: stage-1 mapping is **`arda map`**, not `arda rnaseq map`.
- Tandem-D (D-D) supported throughout (`vdjtools.model.dd`).

### `vdjtools.stats` — diversity, spectratype, usage
- `diversity_stats` (all indices, per sample), **`diversity_cohort(cohort)`** (whole cohort in one
  streamed count-spectrum pass, bit-exact vs per-sample); individual: `observed_richness, chao1,
  chao_e, efron_thisted, shannon_wiener, normalized_shannon_wiener, inverse_simpson, d50`.
- Rarefaction: `rarefaction`, `inext` (Hill q=0/1/2, size+coverage), `inext_batch`,
  `rarefaction_batch`, `inext_coverage`, `asymptotic_diversity`, `coverage`, `sample_coverage`,
  `estimate_d`.
- `segment_usage`, `vj_usage`, `spectratype`, `vj_spectratype` — each takes **`by=["sample_id"]`**
  to compute the whole cohort in one fused `group_by` over a `scan_cohort` LazyFrame (stream-collect).

### `vdjtools.features` — CDR features
`physchem_profile` (region × property), `kmer_profile`, `v_kmer_c_profile`, `load_property_table`,
`DEFAULT_PROPERTIES`.

### `vdjtools.overlap` — overlap + TCRnet (delegates to vdjmatch/seqtree)
`overlap_metrics`, `overlap_pair`, `DEFAULT_KEY`; `fuzzy_overlap`, `fuzzy_overlap_metrics`;
`similarity_overlap`, `similarity_matrix`, `SimilarityMatrices` (TINA / Leinster-Cobbold);
`tcrnet`; `pairwise_distances`, `cluster_samples`; `track_clonotypes`.

### `vdjtools.dynamics` — longitudinal clonotype dynamics (within-donor, across timepoints)
Tests **frequency change across timepoints** within a subject (vs `biomarker`'s incidence across subjects).
- **`test_pair(a, b, *, neff="auto", key=, min_total=, alpha=)`** — per-clonotype paired test:
  two-step-sampling `N_eff` (Ayestaran 2024) + two-tailed Fisher; classes emergent/expanded/
  persistent/contracted/vanishing. Also `estimate_neff`.
- **`test_metaclonotypes(a, b, *, scope="1,0,0,1"|"1,1,1,1", match_v=, match_j=)`** — the same test on
  1-Hamming / 1-Levenshtein CDR3 groups (power for convergent expansions; delegates to `metaclonotypes`).
- **Recapture model** (VDJtrack; Pavlova, Zvyagin & Shugay 2024) `vdjtools.dynamics.capture`:
  `capture_rates` (size buckets singleton/doubleton/tripleton/large × recapture fraction + Beta CI),
  `capture_test` (log-linear `log(recapture) ~ size + group + log(div_ratio)`), `capture_paired_test`,
  `poisson_capture`, `size_class`, `SIZE_CLASSES`.
- **`expansion_test(a, b, *, dispersion="auto", log2fc=, alpha=)`** — edgeR NB-exact caller (TMM +
  qCML common dispersion + beta-binomial exact test; the paper's §2.5 complementary method).

### `vdjtools.preprocess`
`downsample`, `select_top`; `filter_functional`, `filter_frequency`, `filter_segment`,
`filter_by_sample`; `correct` (freq error-correction), `decontaminate`; `pool_samples`,
`join_samples`, `resolve_key`; `correct_vj_usage` (VJ batch-effect: `transform="location"` ComBat
default or `"sigmoid"` = Vlasova 2026 z-score + grand-mean-preserving sigmoid), `apply_vj_correction`
(rescale + roulette-wheel resample the clonotype table to the corrected usage).

### `vdjtools.biomarker`
Incidence contingency testing across a cohort (Emerson 2017 / Howie 2015 / De Witt 2018 / Vlasova 2026).
- `association(cohort, design, *, test=, level_col=, stratum_col=, key=, match=, min_incidence[_frac]=,
  candidates=, alternative=, features=)` — feature-vs-condition; `test` ∈ {`fisher`,`chi2`,`bayes_logodds`,
  `bayes_bf`,`permutation`} (str or list → long output w/ `test` col); category via `level_col` (one-vs-rest),
  paired via `stratum_col` (Cochran–Mantel–Haenszel). Match scope = `key` (`(junction_aa,)`/`+v`/`+v+j`) × `match`:
  - `exact` — the key itself.
  - **`fuzzy`** — 1mm **SEARCH** (Vlasova 2026): `incidence(c) = #subjects carrying ANY feature within `scope` of c`.
    Candidate KEEPS its identity and GAINS incidence; V/J in the key must match exactly. Delegates to
    `vdjmatch.cluster.overlap`. **This is what finds biomarkers.** `key=(junction_aa,v_call)` ≫ `junction_aa`
    alone (real cohort: donor q<0.01 7 → 78). NB `candidates=` is the QUERY set only — the universe stays the
    whole cohort (a candidate's neighbours usually aren't candidates). The search itself depends only
    on `(cohort,key,candidates,scope)`, never the design — `prepare_fuzzy_features(cohort, key,
    candidates=, scope=) -> FeatureFrame` builds it once, reuse via `association(..., features=)`
    across many designs against the same cohort (e.g. one call per HLA gene) instead of paying the
    collect+search cost per call.
  - `1mm` — **CLUSTERING** via `metaclonotypes`: MERGES candidates, tests the group. Different operation;
    belongs *downstream* of a biomarker list (Hamming graph / classifier), not to discovery.
- **Unit + null (the two things that go wrong):** the sampling unit is the **subject** — Emerson beat
  template-weighted abundance head-to-head; weighting a 2×2 by reads is pseudoreplication (Hurlbert 1984).
  If you count **rearrangements** instead (unique nt row = one recombination event), counts hit ~10⁷ → use a
  smooth test (conditional binomial / G-test), never factorials; and the null MUST be the **subject** ratio
  `n_pos/(n_pos+n_neg)`, not the row ratio. Depth differs by arm in real cohorts (FMBA controls are 1.4–1.5×
  deeper/donor), so the two nulls differ ~15–20% and any clonotype not scaling with depth gets exactly that
  much spurious enrichment — hyper-significant at large counts. Depth also biases the subject test the other
  way; for repeated samples of one donor, `preprocess.downsample` each pair to a common read count first.
- **HLA restriction:** measure it **per motif, within cases** (Fisher: carries motif × carries allele). Do NOT
  read it off per-stratum hit counts — the commonest allele wins on power alone (A\*02 is ~half a cohort and
  collects hits restricted by *other* alleles; in HIP-CMV, A\*01 gives the most hits and the weakest specificity
  because it has no dominant CMV epitope).
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
The `vdjtools` typer app. Model: `models`, `generate`, `pgen`, plus the **`vdjtools model <sub>`**
workshop sub-app — `list check template learn build extend rescale export net entropy diversity
compare compare-pgen loglik log`. A model is named as a directory **or** as
`LOCUS[:source[:organism]]` (`TRB`, `TRB:learned`, `TRA:arda:mouse`). **`model check` exits 1 on any
error-severity issue**, so it works as a build gate. Data: **`convert`** (any format →
canonical), **`downsample`**, **`filter`** (`--coding`/`--noncoding`/`--min-freq`/`--v`/`--j`),
**`pool`** (`--join`). Analytics: `diversity`, `overlap`, `segment-usage`, `spectratype`.
Longitudinal/enrichment: `dynamics`, `tcrnet`, `alice`. Inputs auto-detected; **`-o` is
format-aware** — `.parquet`/`.pq` → Parquet, else TSV (or stdout). The per-sample analytics
commands take **`--threads N`** (parallel over samples, `map_samples`) and **`--cohort DIR`** (one
streamed pass over a pre-ingested `scan_cohort` Parquet dataset).

### Notebooks (`pip install "vdjtools[examples]"` → `marimo edit examples/<name>.py`)
`model_explorer` (recombination Bayes net), **`model_workshop`** (custom germline → learn → check →
compare → BIC → diversity → extend → rescale), `biomarker_explorer` (Emerson association /
co-occurrence), **`vaccination_tracking`** (YFV/flu/TBEV clonotype dynamics + recapture model),
**`aging`** (cohort-streaming diversity / clone-size / spectratype vs age), **`ankspond_motif`**
(AS27 TRBV9 motif, disease-vs-B27-carriage). Local-first data (`./` → `~/hf/` → HuggingFace).

## Conventions
- **arda germline is the single source of germline truth** — resolve V/D/J germline + CDR3 anchors
  by allele name via `model.reference`. Never mix germline sources within a model.
- Delegate: overlap/TCRnet → **vdjmatch**; annotation/markup/germline → **arda**; search/e-value →
  **seqtree**. Only `arda`/`vdjmatch`/`seqtree` are dependencies; everything else under `~/vcs/code/`
  is read-only oracle/reference (never modify).
- Native C++ only for the hot loops (Pgen DP, generation sampler, EM E-step) via the single `_core`
  ext; everything else is polars.
