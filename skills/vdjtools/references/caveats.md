# Workflow caveats

Read the relevant section before choosing defaults, interpreting results or changing
model/feature code. The linked documentation owns full signatures and measured results.

## Input, pairing and missing evidence

- Bulk estimators assume one row per clonotype. Preserve read/count semantics when
  converting, pooling or streaming; do not average chunk-local frequencies.
  `io.convert.read_mitcr(..., collapse=False)` preserves repeated source rows,
  missing calls and ambiguous calls when that is the intended input unit.
  See [clonotype premise](https://docs.isalgo.dev/vdjtools/preprocessing.html#one-row-per-clonotype-the-premise-every-estimator-rests-on)
  and [bounded batches](https://docs.isalgo.dev/vdjtools/usage.html#read-a-large-sample-in-bounded-batches).
- Productive rearrangement, functional germline gene and length filters are distinct.
  Missing SHM/isotype annotations are missing evidence, not evidence of absence.
  See [filtering axes](https://docs.isalgo.dev/vdjtools/preprocessing.html#filtering-three-axes-and-they-are-not-the-same-axis)
  and [SHM](https://docs.isalgo.dev/vdjtools/usage.html#somatic-hypermutation-b-cells).
- Start cell workflows with `sc.read_10x` or `sc.read_airr_cell`, then resolve and
  pair chains. `paired_pgen` defaults to V/J conditioning and representative-gene
  resolution. Unknown calls or missing junctions produce null chain scores; the
  paired score needs both chains. Investigate whole-locus null warnings.
  `condition_vj=False` deliberately marginalizes; `resolve_genes=False` disables
  representative selection for multi-allele models, while single-entry genes still
  accept suffix aliases. See [pairing and QC](https://docs.isalgo.dev/vdjtools/singlecell.html#pairing-and-qc).
- `sc.push_obs` requires one row per cell; aggregate multiple receptor pairs first.
  Match sample/barcode conventions across scRepertoire and expression data.
  See [object augmentation](https://docs.isalgo.dev/vdjtools/singlecell.html#augmenting-an-existing-object)
  and [scRepertoire](https://docs.isalgo.dev/vdjtools/singlecell.html#screpertoire-r).

## Feature identity, coverage and memory

- Preserve feature order, vocabulary, clone weights, coverage and `embedding_space`;
  retain `rotation_space` when comparing configured projections. Do not infer
  compatibility from column count. Reuse a reference bundle's saved vocabulary,
  including its declared loci. See [reference rotations](https://docs.isalgo.dev/vdjtools/features.html#reference-rotation-bundles).
- Legacy signature coverage is a runtime choice; separate cohort minima need not
  be comparable. Configured corpora store the measurement policy and default to
  its saved coverage target; `own` is refused. Use the same explicit shared target
  across compared tables. See [legacy coverage](https://docs.isalgo.dev/vdjtools/signature.html#the-coverage-target-is-a-runtime-argument)
  and [configured fitting](https://docs.isalgo.dev/vdjtools/signature.html#fitting-configurable-bases).
- Configured length bins include junction anchors; motif windows remove the C/F/W
  anchors. Gapped motifs have three residues over four positions and require valid
  amino acids even at the skipped position. Missing genes and uncalled isotype
  retain residual categories. See [configurable bases](https://docs.isalgo.dev/vdjtools/features.html#configurable-statistical-bases).
- Feature-table PCA and robust corpus transformation have different scaling and
  imputation contracts; use `TablePCA` and `Corpus` respectively. Per-locus corpus
  rotations preserve missing-locus holes. Selecting output columns does not avoid
  computing their source groups. See [column selection](https://docs.isalgo.dev/vdjtools/features.html#choosing-feature-columns)
  and [signature output selection](https://docs.isalgo.dev/vdjtools/signature.html#columns-selects-output-it-does-not-skip-work).
- Sparse randomized fitting is not constant-memory: sparse entries, rank workspaces
  and output size still grow. Already-loaded frames remain resident; prefer
  repeatable deferred readers. `fit_cohort(..., n_jobs=B, threads=B)` joins sample
  workers before BLAS fitting; zero fitting threads selects allocated cores,
  while the fitting default is one. See [resource limits](https://docs.isalgo.dev/vdjtools/features.html#resource-limits)
  and [parallelism](https://docs.isalgo.dev/vdjtools/signature.html#parallelism).
- Keep stored corpus source/fingerprint identities. A new bundled model default
  does not authorize corpus refitting or substituting weights. Do not cache computed
  feature tables or other analysis intermediates; keep requested final artifacts
  distinct from reusable raw inputs and frozen model parameters. See
  [building corpora](https://docs.isalgo.dev/vdjtools/signature.html#building-a-corpus).

## Model references, inference and interpretation

- New ARDA scaffolds select one usable functional major allele per gene: *01,
  otherwise the lowest usable numeric allele. `gene_level=False` retains multiple
  alleles. Between-gene ambiguity remains. Never mix ARDA and frozen OLGA germlines
  inside a model; check anchor coordinates when importing libraries. V and J each
  contribute at least one nucleotide. See [custom germlines](https://docs.isalgo.dev/vdjtools/model.html#custom-germline-libraries)
  and [OLGA caveats](https://docs.isalgo.dev/vdjtools/model.html#known-quirks-of-the-olga-models).
- First-insertion-base distributions are separate from Markov transitions; a
  one-base insertion supplies first-base evidence. Biological deletion counts can
  be negative for P nucleotides. Preserve the model's schema and transition
  orientation when importing/exporting: `(from_nt, to_nt, p)` means
  `P(next=to_nt | previous=from_nt)`, unlike OLGA's column-stochastic storage. See [custom germlines](https://docs.isalgo.dev/vdjtools/model.html#custom-germline-libraries)
  and [table interchange](https://docs.isalgo.dev/vdjtools/model.html#exporting-and-importing-tables).
- `d_history` normalizes complete single/tandem histories within J. A fully deleted
  single D is allowed; tandem histories retain at least one nucleotide from each D.
  TRB D1–J2 and D1–D2–J2 are valid. Recognized human TRD forward histories follow
  TRDD1–TRDD2–TRDD3; do not infer custom-species order from names or impose a universal
  10% tandem ceiling. See [joint histories](https://docs.isalgo.dev/vdjtools/model.html#joint-single-tandem-d-histories).
- `candidate_masks` is a finite ARDA V/J gene-proposal search, not exhaustive
  inference. It expands proposals to all functional model templates; missing V/J
  proposals remain unrestricted, as do D/D2. Candidate restriction and scenario
  pruning are separate approximations: `prune_ratio=0` does not undo restricted
  candidate support. See [corpus model building](https://docs.isalgo.dev/vdjtools/model.html#building-from-the-read-corpus).
- Matching germline bases do not prove zero insertion: trimming followed by
  coincident N bases remains a possible history. Do not fix longest alignments or
  add a hard reassignment-length cap. Constraints belong before fitting, not in
  a posthoc repair of fitted probabilities. See [learning](https://docs.isalgo.dev/vdjtools/model.html#learning-from-your-own-sequences).
- `infer_native(..., error_rate=0.01)` enables fitted retained-germline substitution
  likelihood; `None` uses exact matching. `report.loglik` averages scoreable inputs,
  not the full cohort when unsupported observations remain. Report both
  `n_scoreable` and `n_sequences`, plus `pruned_mass_mean`, `pruned_mass_max` and
  `n_pruned_observations`. Checkpoint resume restores recorded error/pruning policy;
  preserve candidate support or explicitly report changes. Resume availability is
  not permission to reuse intermediates where a fresh run is required. See
  [fit monitoring](https://docs.isalgo.dev/vdjtools/model.html#watching-and-surviving-a-long-fit).
- Training likelihood improvement is not independent validation. Compare supported
  populations and use held-out/oracle checks appropriate to the requested work;
  do not rerun completed analyses or add unrequested controls. See
  [checking models](https://docs.isalgo.dev/vdjtools/model.html#checking-a-model)
  and [likelihood/BIC](https://docs.isalgo.dev/vdjtools/model.html#likelihood-and-bic).
- Rank-one tandem accuracy is empirical, not a per-query guarantee or proof of
  unchanged downstream decisions. Exact NT, generation, inference and Hamming
  paths remain distinct from this opt-in AA approximation. See
  [approximation limits](https://docs.isalgo.dev/vdjtools/model.html#optional-approximation-of-tandem-histories).
- Selection calibration is model-, organism- and locus-specific. Calibration counts
  are independent recombination events, not sequencing reads. Modern ALICE uses
  model-bound Q; historical OLGA TRB Q must not transfer silently to other contexts.
  Optional selection profiles carry their own fitted shape and Q and do not change
  ordinary Pgen. See [selection mass](https://docs.isalgo.dev/vdjtools/model.html#explicit-selection-mass-scaling),
  [ALICE exposure](https://docs.isalgo.dev/vdjtools/model.html#alice-exposure-convention)
  and [selection profiles](https://docs.isalgo.dev/vdjtools/model.html#opt-in-coarse-selection-profiles).
