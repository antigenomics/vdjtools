# ROADMAP

vdjtools 2.x — phased Python + C++ rewrite. Gitflow: `master` (tagged releases) ← `dev`
(integration) ← `feature/*` (one per phase). Merge a feature to `dev` when its tests are
green; cut tags off `master` as phases land.

Legacy v1.x is preserved on the `legacy-1.x` branch and under all existing tags
(`v0.0.1`…`1.2.1`). The v2 history is an orphan root.

**Current release: `v2.1.0`** (2026-07) — model engine + analytics Phases 1–8, native Pgen
speedups, typer CLI, example notebooks. Published to PyPI.

The full approved v2.0 design lives in `~/.claude/plans/i-want-to-complenely-gleaming-snail.md`;
the Phase 6 biomarker sub-plan in `~/.claude/plans/biomarker-fisher-emerson.md`. `SOURCES.md`
tracks all data provenance.

## Shipped (v2.0.0-alpha → v2.1.0)

| # | Branch | Scope | Status |
|---|---|---|---|
| 0 | (master root) | Git surgery, GPLv3 relicense, repo scaffold, CI/docs/publish, AIRR+polars IO layer (vdjtools/AIRR/parquet readers) | **done** |
| 1 | `feature/model-engine`, `feature/model-dd-native` | Native V(D)J model: polars marginals (single-D **and** tandem D-D), arda scenarios, contig stitching, EM from OLGA-synthetic + real out-of-frame bootstrap, native Pgen/generation/E-step (pybind11), aa transfer-matrix Pgen, bundled 7-locus models (OLGA + learned), model diagnostics. Exact vs OLGA on all 7 loci. | **done** |
| 2 | `feature/repertoire-stats` | Diversity (Chao1/ChaoE/Shannon/Simpson/d50/Efron–Thisted, exact+resampled+rarefaction+quantile), spectratype, V/J/VJ usage | **done** |
| 3 | `feature/cdr-features` | CDR physicochemical profiles (regions × properties), k-mer / V+k-mer summaries | **done** |
| 4 | `feature/overlap-tcrnet` | Sample overlap + TCRnet via vdjmatch/seqtree; **similarity-aware** overlap (TINA / Leinster-Cobbold), pairwise-distance matrix, clustering/MDS, tracking | **done** |
| 5 | `feature/preprocess` | Downsampling, error-correction, decontaminate, filters, VJ-usage batch-effect correction, pool/join | **done** |
| 6 | `feature/biomarker-assoc` | Fisher incidence association (V/J-match + exact/1mm) vs HLA/condition; metaclonotype grouping; Emerson-2017 CMV/HLA benchmark validated vs VDJdb | **done** |
| 6b | `feature/biomarker-cooccurrence` | Full association framework (Fisher / χ² / Bayesian / Beta-Binomial BF / permutation; binary / per-HLA-allele / Cochran–Mantel–Haenszel stratified conditions) + feature-vs-feature **co-occurrence** (α-β pairing, same-chain co-specificity; θ-lift + Fisher + FDR + Poisson E-value, **depth-conditioned by default**); airr_hip / airr_covid19 / airr_covid19_vacc benchmarks. **v2.7.0** | **done** |
| 7 | `feature/singlecell-interop` | AIRR Cell / 10x paired-chain interop, single-cell metadata bridge, QC (doublet/mispairing), cluster evaluation, AnnData bridge | **done** (paired α/β Pgen residual → Phase 12) |
| 8 | `feature/cli-docs` | typer CLI (pgen/generate/diversity/spectratype/segment-usage/overlap/models), 6 example marimo notebooks, **v2.1.0** release | **done** (full Sphinx API docs residual → Phase 12) |

## Upcoming

**Released — `v2.2.1`** (patch, 2026-07-13): docs user guide with runnable per-module examples +
version in the navbar; MiXcr converter now reads the **C-gene / BCR isotype** hit; new **TRUST4**
(`read_trust4`) and **arda** (`read_arda`, delegates to `read_airr`) converters; plus the two
post-2.2.0 audit bug fixes (metadata-column collision, immunoSEQ count fallback).

**Released — `v2.2.0`** bundled Phases **11 → 12 → 13** (owner directive 2026-07-13):
AIRR `junction_nt`/`junction_aa` rename → close the v2.1 residual gaps → model-engine residuals.

The **BCR track — Phases 9 (MiGEC UMI correction) + 10 (BCR SHM/lineage) — is deferred** to a
later release: it is the larger, higher-risk body of work and Phase 10 depends on Phase 9. Their
designs are kept below so the track is ready to pick up. Phase numbers are stable (referenced
elsewhere) and sections stay in numeric order; each heading is tagged with its release target.

### Phase 9 — `feature/umi-migec` — MiGEC UMI error correction

**Status: deferred** to a post-v2.2 release (BCR track, item 1 of 2).

Re-implement the MiGEC UMI/molecular-barcode error-correction pipeline (originally Java;
Shugay et al. Nat. Methods 2014) as Python-first (native C++ for the consensus hot loop only),
emitting AIRR clonotype tables consumable by the rest of vdjtools. **Prerequisite for the BCR
lineage work** — the UMI-tagged BCR datasets need molecular error correction before SHM/lineage
inference, or PCR/sequencing errors masquerade as hypermutations.

Stages (mirror MiGEC):
- **Checkout** — UMI/barcode extraction + demultiplex from raw FASTQ (adapter/barcode patterns).
- **Histogram** — UMI coverage distribution; pick the over-sequencing threshold.
- **Assemble** — group reads by UMI into MIGs (molecular identifier groups), build per-MIG
  consensus, drop under-covered MIGs. *(hot loop → `_core` if profiling warrants)*
- **Correct** — collapse UMI/consensus errors (1-mismatch barcode collision, hot-spot correction).
- Output → AIRR reads/clonotypes for arda annotation + model/lineage.

Delegate downstream V(D)J markup to **arda**. Overlap with Phase 12 gap #2: that gap ships a
read-only **MiGec output-table converter** now (v2.2.0); when this phase lands it supersedes that
reader — this phase *is* the MiGEC processor, not just a parser of its output.

### Phase 10 — `feature/bcr-lineage` — BCR SHM & lineage/tree reconstruction

**Status: deferred** to a post-v2.2 release (BCR track, item 2 of 2; depends on Phase 9).

Clonal-lineage and somatic-hypermutation (SHM) reconstruction for BCR/IGH. Standard TCR-style
exact/1mm clustering fails under heavy SHM; this phase builds a Pgen- and mutation-model-aware
lineage engine. Two hard problems the design must confront:

1. **Heavy hypermutation** — clonally related variants can differ at many positions; naive
   distance clustering breaks.
2. **Unknown common ancestor** — parsimony on `junction_aa` is ill-posed because the naive
   ancestor of the junction (especially the non-templated N region) is not observed, so the
   tree cannot be rooted by assuming the ancestor.

Approach (owner's design):

- **(a) Pgen as SHM-load proxy / rooting signal.** Lower Pgen ⇒ more hypermutated IGH: a naive
  (unmutated) sequence sits in the generative model's high-Pgen region, while SHM drives the
  observed sequence away from any plausible germline rearrangement and lowers its model Pgen.
  Pgen therefore orders lineage members by mutational load and identifies the likely near-germline
  root (highest Pgen = closest to the naive ancestor). Uses the native IGH model in `vdjtools.model`.

- **(b) Clonal-relationship via non-templated N-stretch matching.** Two clonotypes sharing a long
  run of identical *non-templated* insertion nucleotides (e.g. a 5–7 nt consecutive random-base
  k-mer in the N1/N2 region) almost certainly descend from the same VDJ rearrangement — random
  insertions are a per-rearrangement fingerprint. Candidate lineage edges come from matching long
  shared N-stretch k-mers.
  - **Markup caution**: the match must be on the genuinely non-templated segment. Templated
    **D-gene** nucleotides are shared germline and would spuriously inflate shared-k-mer counts,
    merging unrelated clones. Use proper V(D)J markup — potentially an HMM / Pgen-based
    segmentation that **excludes D** (and V/J germline) so only true non-templated insertions form
    the fingerprint. Ties directly to the model engine's scenario markup: the same event structure
    that yields Pgen yields the N-region boundaries.

- **(c) Likelihood-ratio test for clonal relationship → edges, SHM profiles, trees.** For a
  candidate pair compare:
  - H0 (independent rearrangements): `P_gen(c1) · P_gen(c2)`
  - H1 (one rearrangement + SHM):    `P_gen(c1) · P_mut(c1 → c2)`

  where `P_mut` is a context-aware SHM mutation model (per-position / 5-mer context, S5F-style).
  If H1 ≫ H0 the pair is clonally related; orient the edge from the more-ancestral (higher-Pgen)
  member. Aggregated over a clone these give **SHM profiles** (mutation-frequency spectra, WRC/GYW
  hot/cold-spot context) and a weighted graph from which lineage **trees** are inferred (root by
  Pgen, edges by LRT).

Deliverables (`vdjtools/lineage/` or `vdjtools/bcr/`):
- context-aware SHM mutation model `P_mut` (learn/bootstrap; native for the hot loop),
- N-region fingerprint extraction from D-excluded model markup,
- clonal-relationship LRT scorer combining native Pgen + `P_mut`,
- lineage-graph builder + tree inference (Pgen rooting, LRT edges), SHM-profile summaries.
- Delegate germline alignment/markup to **arda**, fuzzy search to **seqtree**; Pgen/model from
  `vdjtools.model`.

**Test data** — verified via PubMed, pinned in `SOURCES.md` (§Phase 10); all owner (Shugay)-co-authored
BCR datasets. Paths TBD (fill `SOURCES.md` when each is actually fetched — do not fabricate a path first):
- **Mikelov et al. 2022, *eLife*** — [doi:10.7554/eLife.79254](https://doi.org/10.7554/eLife.79254) —
  longitudinal UMI-tagged memory-B/ASC BCR repertoires with clonal-lineage/phylogenetic analysis: the
  canonical Phase-9-MiGEC + Phase-10-lineage target. *(Owner said "Mikelov allergy paper"; the literal
  allergy match is Mikelov et al. 2025 Nat. Immunol. peanut-OIT, doi:10.1038/s41590-025-02323-3 — but
  that is a **TCR/scRNA** study without Shugay, so the eLife BCR paper is used here. Confirm intent.)*
- **TCGA** repertoires via **Bolotin et al. 2017, *Nat. Biotechnol.*** —
  [doi:10.1038/nbt.3979](https://doi.org/10.1038/nbt.3979) — RNA-seq-derived TCR/BCR (MiXCR method), TCGA IGH.
- **Grimsholm et al. 2020, *Cell Rep.*** — [doi:10.1016/j.celrep.2020.02.022](https://doi.org/10.1016/j.celrep.2020.02.022)
  — CD27-dull/bright memory-B VH repertoires + SHM (owner wrote "CD20-dull"; the paper is **CD27**-dull).
- **Grimsholm et al. 2023, *Cell Rep.*** — [doi:10.1016/j.celrep.2023.112446](https://doi.org/10.1016/j.celrep.2023.112446)
  — CVID peripheral B-cell selection, Ig-seq.

Risks: SHM-model identifiability; N-fingerprint false merges under short insertions; rooting when
no near-germline member was sampled; indel-bearing SHM (not just substitutions).

### Phase 11 — `feature/airr-junction` — AIRR field-name consistency

**Status: DONE (`v2.2.0`, item 1 of 3).** Canonical columns renamed `cdr3_aa`/`cdr3_nt` →
`junction_aa`/`junction_nt` across schema/io/model/stats/features/overlap/preprocess/biomarker/sc/cli
+ tests + docstrings (paired sc cols → `alpha_junction_aa`/`beta_junction_aa`). **Decision locked**:
owner's `junction_nt`/`junction_aa` internally + on output; readers stay liberal, accepting strict-AIRR
`junction`, legacy `cdr3_nt`/`cdr3_aa` (compat shim), IMGT `cdr3`, native `cdr3nt`/`cdr3aa`, and 10x
fields as input aliases. Sequence-arg params in `model.pgen`/`native`/`stitch` intentionally unchanged
(they name a sequence, not a column). 315 tests green; adversarial diff review clean.

Rename the legacy vdjtools clonotype columns `cdr3nt` / `cdr3aa` → **`junction_nt` / `junction_aa`**
throughout (schema, IO, model, stats, features, overlap, preprocess, biomarker, sc, cli, tests,
docs, notebooks). vdjtools' "CDR3" is anchor-inclusive (starts Cys104, ends Phe/Trp118), i.e. the
AIRR **junction** — so this is a correctness alignment, not just cosmetics.

- **Decision to lock**: strict AIRR uses `junction` (nt) + `junction_aa`; the owner specified the
  parallel `junction_nt` / `junction_aa`. Default to the owner's `junction_nt`/`junction_aa`; flag
  the strict-`junction` alternative at the top of the phase. (Ref: `immunogenomics-conventions` —
  junction includes both conserved anchors; IMGT CDR3 excludes them.)
- Breaking change → version bump (**v2.2.0** or **v3.0.0** depending on the deprecation policy);
  consider a compat shim reading old column names with a deprecation warning.

### Phase 12 — `feature/v2.1-residuals` — close the v2.1 gaps (the "1/2/3" gaps)

**Status: `v2.2.0` (next minor), item 2 of 3.**

1. **Full Sphinx API docs** — only `docs/index.rst` + `docs/api.rst` exist and just `model` has an
   `automodule` ref. Write per-subpackage API pages for `io, stats, features, overlap, preprocess,
   biomarker, sc, cli` (and `model`). Zero-warning `sphinx-build -W` gate (`sphinx-docs` rules).
2. **Legacy input-format converters** — **DONE**: `vdjtools.io.convert` reimplements the legacy
   Groovy parsers for MiXcr (v1/2 + v3/4 header dialects), MiGec, ImmunoSeq v1/v2, ImgtHighVQuest,
   Vidjil (JSON), RTCR → the canonical junction frame; `sniff_format`/`io.read(fmt="auto")` detect
   and dispatch them. Fixtures copied to `tests/python/fixtures/legacy/`; conformance oracles in
   `test_convert.py` (migmap is already read by `read_vdjtools`). Adaptive→IMGT gene conversion +
   bidirectional `translate()` ported verbatim from `CommonUtil`.
3. **Paired α/β Pgen** (Phase 7 residual) — **DONE**: `vdjtools.sc.paired_pgen` adds
   `pgen_alpha`/`pgen_beta`/`pgen_paired` (= `Pgen(α)·Pgen(β)`) to a paired-chain frame using the
   **native model** + bundled per-locus models (no `vdjmatch` dependency — supersedes the
   `evalue.paired` route). Loci inferred from V-call prefixes; V/J conditioning when the call
   matches a model allele, else marginalised. `test_sc_pgen.py`.

### Phase 13 — `feature/model-residuals` — model-engine residuals (the "phase-1 4/5" items)

**Status: `v2.2.0` (next minor), item 3 of 3.**

4. **arda full-length V/J germline helper** — **DONE** (P1c prerequisite closed):
   `model.reference.load_full_vj_germline(organism)` recovers full-length V/J germline from arda
   2.5.2's bundled V–J scaffolds (`alleles.fasta`, sliced at `v_sequence_end`/`j_sequence_start`
   via `arda.annotate.reference.load_reference`); `arda_full_germline(locus, organism)` returns the
   stitch-ready `{(seg, allele): (full_germline, anchor)}` with anchors derived self-consistently
   from the CDR3-region germline (verified: `full_V[anchor:]` == CDR3-region germline for all 895
   functional V / 132 J alleles). arda-native models can now stitch full contigs. `test_full_germline.py`.
5. **Native perf** — **DONE (the exact, safe win):** `native.pgen_aa_batch(model, seqs, v=, j=,
   mismatches=, threads=)` parallelizes Pgen / Hamming-1-ball across sequences (GIL released, disjoint
   workers → **bitwise-identical** to the per-sequence calls, thread-count-invariant); **11.3× exact /
   11.6× 1-mm on 16 cores**. Also releases the GIL on the single `pgen_aa`/`pgen_aa_hamming1` bindings
   so pure-Python threading works too. `test_native_pgen_batch.py`. This is the real-workload speedup
   (Pgen over many clonotypes — TCRnet / biomarker matching) and matches the engine ethos ("the clean
   exact speedup is read-parallelization, not gating").
   - *Set aside (not rushed):* the literal VJ / Hamming-1 **codon-boundary sweep** (collapse the L+1
     TM passes to ~1). Analysis of `pgen_aa_vj`: the V/J combine boundary **migrates with the delJ
     sum**, so a wildcarded codon is left-of-boundary for some delJ and right-of-boundary for others
     — there is no clean O(1)-per-codon leave-one-out. Forcing it is a high-risk rewrite of
     numerically-delicate code for a non-bottleneck; batch parallelization delivers the same
     end-goal (fast large-scale 1-mm Pgen) safely.
   - *Still low-priority:* (b) native generation sampler — Python is already fast. *(`estep_batch`
     read-parallelization was already done.)*

## Design principles

AIRR schema + polars everywhere, minimal OO. Python-first — native C++ (pybind11, single `_core`
ext) only for hot loops (Pgen DP, generation sampler, EM E-step, UMI consensus). Delegate
overlap/TCRnet to **vdjmatch**, annotation/markup/germline to **arda**, search/e-value to
**seqtree** rather than reimplementing. **arda's germline library is the single source of germline
truth** — all V/D/J germline + CDR3 anchors resolve from arda by allele name
(`model.reference.load_germline`), so annotation ↔ scenarios ↔ stitching ↔ Pgen ↔ lineage share one
coordinate frame. Never modify the dependency libraries (`arda`, `vdjmatch`, `seqtree`); everything
else under `~/vcs/code/` (mirpy, OLGA, IGoR, …) is read-only oracle/reference.
