# Changelog

Notable changes to vdjtools v2. Releases before 3.0.0 are recorded in the git tags
(`v2.5.0` … `v2.9.0`) and their commit history.

## 3.3.0 — 2026-08-12

The recombination model becomes a workshop: buildable on your own reference, checkable, comparable,
scoreable and extendable. See the new [user guide](https://docs.isalgo.dev/vdjtools/model.html) and
`examples/model_workshop.py`.

### Added

- **Custom V(D)J reference libraries.** `model.io.from_germline(germline_df, locus=...)` builds a
  model on any germline library; `from_arda` is now a six-line wrapper over it (its output is
  byte-identical, verified table by table on TRA/TRB/TRG). `reference.read_germline_fasta(v, j, d,
  anchors=)` reads your own FASTA — segment comes from which argument a file is passed as, so no
  header convention is assumed — and `reference.validate_germline` audits it. The audit includes
  the two anchor-frame checks (V starts on a Cys codon, J ends on Phe/Trp) that catch the most
  damaging custom-library mistake: a CDR3 anchor one codon off shifts every deletion profile by a
  constant and nothing downstream complains.
- **`model.check.check_model`** — a consistency audit returning a tidy issue frame
  (`severity check event segment allele detail value`) rather than raising, so every problem in a
  model is visible at once. Covers normalization and probability range, alleles missing from the
  germline (or absent from the marginals), functional genes stuck at `P = 0`, unreachable deletion
  mass, incomplete dinucleotide tables, and VDJ/VJ event-set mismatches. Deletion reachability is
  derived from the Pgen DP (`ndel = len(cut_segment) − contributed − max_palindrome`) and reported
  as the **fraction** of each allele's mass that is lost, ranked — up to 25% for `IGKV3-20*01` in
  the bundled IGK model. `vdjtools model check` exits 1 on any error-severity issue.
- **`model.score`** — likelihood and diversity. `model_fit` reports log-likelihood, free parameters,
  AIC and BIC; it uses **nucleotide** Pgen, because `Σ Pgen_nt = 1` makes the log-likelihood proper
  (amino-acid Pgen sums only the in-frame, stop-free fiber, so its missing normalizing constant
  differs between models and it is a relative score only). A sequence the model cannot generate is
  counted in `n_scoreable`, never turned into `-inf`. `free_params` counts **occupied cells**, not
  rows, and drops undefined and unreachable conditional groups — the difference between ~700 and
  ~3,600 parameters for human TRB's `v_3_del` alone. Also `pgen_frame`, `compare_pgen` +
  `pgen_summary` (KS, Spearman, and the headline one-sided coverage counts), and `pgen_spectrum`.
- **Information content and total diversity.** `analyze.total_entropy` gives each recombination
  event's contribution to the scenario entropy (the dinucleotide term is `E[length] × H_step`), and
  `score.diversity` adds a Monte-Carlo sequence entropy with its standard error plus both Hill
  numbers — `2^H` and `1/E[Pgen]`. Human TRB: ~52 bits per rearrangement, ~45 bits per sequence,
  ~3·10¹³ effective sequences.
- **Model comparison.** `analyze.compare_models` reports per-event total variation, `tv_max` and
  Jensen-Shannon over the union of both models' realizations with zero fill, weighted by the
  parent's marginal; `by="gene"` bridges different germline namespaces, and an event factorized
  differently in the two models is flagged `schema_differs` rather than joined. Plus
  `compare_usage` and `compare_net_dot` (a bnlearn `compare_networks`-style graph).
- **Training log.** Every EM fit now appends a run to `model.training["runs"]`, persisted beside the
  model as a `training.json` sidecar and readable as a table with `infer.training_frame`. Both the
  field and the sidecar are optional, so every previously-saved model still loads (with
  `training is None`).
- **`infer.infer_frame`** fits from a clonotype frame, building the per-read V/J masks for you, and
  **`infer.extend_alleles`** adds alleles from a larger germline library.
- **`data.build_all`** runs the full corpus pipeline — fetch FASTQ, map with arda, collapse, EM —
  parallel across chains, exposed as `vdjtools model build`. Two arda-mapped clonotype examples
  (human TRA 34,238 and TRB 100,000 out-of-frame clonotypes, 2 MB total) ship in
  `tests/python/fixtures/model_reads/` as gzipped FASTA with the V/J/D calls in the header, and
  load offline with `data.load_prepared` — no network, no arda, no mmseqs2.
- **`vdjtools model` CLI sub-app** with 13 subcommands. A model is named as a directory or as
  `LOCUS[:source[:organism]]`.
- **Table export/import**: `marginals_frame` / `set_marginals`, and `save_model(..., fmt="tsv")`
  with format auto-detection on load, so a hand-edited TSV directory is a first-class model input.

### Changed

- **`pgen_nt` now releases the GIL**, the one Pgen binding that still held it after the Phase-13
  batch work. Threaded nucleotide Pgen is **11.7× faster** on this Mac and bitwise-identical to the
  serial result.
- `extend_alleles` **preserves each pre-existing gene's total usage**. Alleles of one gene are
  alternative versions of the same gene — a diploid carries at most two — so a richer library must
  split a gene's mass more finely, never multiply it. Seeding each new allele at its gene's average
  without this correction moved gene-level V usage on human TRB by up to 6 percentage points,
  silently reweighting every Pgen through those genes.

### Fixed

- **`collapse_alleles` could give a gene a germline it could not use — in the default path.**
  `load_bundled(..., collapse=True)` picks one representative allele per gene and averages the
  other alleles' conditionals onto it. The representative was chosen by usage alone, but IMGT ships
  some alleles with a **truncated** CDR3-region germline: human `IGKV3-20*02` is 11 nt against
  `*01`'s 30 and carried the higher learned usage, so it became the gene's germline — relabelled
  `*01`, which was doubly misleading — and 25% of the gene's own averaged deletion distribution
  landed on trims the 11-nt germline cannot reach. The Pgen DP never visits those, so that quarter
  of the probability vanished from every Pgen through IGKV3-20 instead of being redistributed.
  The representative is now chosen by **germline length first, usage second** (alleles of a gene are
  near-identical through the CDR3 region, so a large length gap means an incomplete database entry),
  and the collapsed deletion conditionals are **projected onto the representative's reachable
  support and renormalized**. Every bundled model is now clean at `collapse=True`.
- **A failed `arda` run reported nothing but an exit code.** `annotate_reads` passed
  `capture_output=True` and let `CalledProcessError` propagate, so arda's own message was swallowed
  — which is precisely how the CLI rename above stayed invisible. It now raises with arda's stderr
  and the installed `arda-mapper` version.
- **`annotate_reads` was calling an arda CLI that no longer exists.** It shelled out to
  `arda rnaseq map -o …`, but arda 2.19 turned `rnaseq` into the full map→assemble→correct preset
  with no stage positional and `-p/--out-prefix` in place of `-o`, so every real invocation exited
  2 — i.e. the whole model-training pipeline was broken against the installed arda. Stage-1 mapping
  is `arda map`; the pin is now `arda-mapper>=2.19.0`. Worth recording that a `--help` smoke test
  would *not* have caught this: typer short-circuits `--help` before argument parsing, so
  `arda rnaseq map --help` still exits 0.
- **Ambiguous junction bases crashed EM** with a bare `KeyError: 'N'` from inside the native
  encoder. `infer_frame` and `build_model` now substitute `A` by default and warn with the count
  (`ambiguous=None` drops the clonotype instead) — see `infer.sanitize_junctions`.
- `reference.read_fasta` wraps arda's FASTA parser with gzip support; arda's opens with plain
  `open()`, so a `.gz` reached it as mojibake and died on the first byte.
- `tests/python/test_io_hf.py::test_control_native_schema_capped` called `list_repo_files` outside
  the `hf` fixture's guard, so an offline run failed instead of skipping.

## 3.2.0 — 2026-08-09

### Fixed

- **Adaptive/immunoSEQ gene names were wrong for 100 of the 161 tokens** seen in the IMMREP25
  release + the pairSEQ mock cohort (22,058 of 44,000 gene calls), and **every one of those 100
  outputs is a gene name absent from the IMGT human reference**. `_adaptive_to_imgt` normalised
  Adaptive tokens with a global `re.sub(r"0([1-9])", r"\1", …)`, which always re-emits the trailing
  `-01` as an IMGT *subgroup*: `TCRAJ39-01 → "TRAJ39-1"` (no human TRAJ gene has a subgroup),
  `TCRBV09-01 → "TRBV9-1"`, `TCRBD01-01 → "TRBD1-1"`. Slash ties (`TCRBV03-01/03-02`), family-only
  calls (`TCRBV20-X`) and co-locus names (`TCRAV38-02`, IMGT `TRAV38-2/DV8`) were passed through
  verbatim. Any consumer resolving gene names against a germline reference silently lost the rows —
  tcrdist3 dropped 100 % of both cohorts.

  Whether a token's trailing group is a subgroup or an allele is a per-family fact (`TCRAV01-01` =
  `TRAV1-1`, but `TCRAV22-01` = `TRAV22`), so no regex can decide it. `read_immunoseq` now resolves
  V, D **and** J calls through a shipped CDR-validated table, `resources/adaptive_imgt_map.tsv`
  (163 tokens; the choice among candidates is decided by exact-matching germline CDR1+CDR2 —
  provenance in `SOURCES.md`, rationale in `appendix/adaptive_imgt_map.md`). Tokens outside the
  table fall back to the legacy rewrite, so unknown input behaves exactly as before.
  The legacy Groovy `CommonUtil.extractVDJImmunoSeq` has the same defect — a v1 bug inherited by
  v2, not a porting regression.
- **`read_mixcr` accepts every MiXcr count spelling** (`cloneCount`, `readCount`,
  `uniqueTagCountMolecule`) — a v4 `-readCount` export used to raise outright.

## 3.1.2 — 2026-07-30

### Changed

- Bumped the `vdjmatch` floor to `>=0.1.2` — a fresh install of 3.1.1 could resolve `vdjmatch`
  0.1.1, whose `cluster.overlap()` raised `SchemaError` on any query with zero fuzzy matches
  (`a_idx`/`b_idx` defaulted to a `Null` dtype with no hits to infer from, then failed to join
  against the `Int64`-typed lookup frames); this broke `vdjtools.overlap.fuzzy.fuzzy_overlap` and
  CI's `test_overlap_fuzzy.py::test_fuzzy_no_match_empty_and_zero_metrics`. Fixed upstream in
  vdjmatch 0.1.2.

## 3.1.1 — 2026-07-30

### Fixed

- **`biomarker.association(match="fuzzy")` re-ran the full-cohort search on every call**, even
  when a caller tests many phenotype designs (e.g. one per HLA gene) against the same
  cohort/key/candidates/scope — the search depends only on those, never on the design. At
  full-corpus scale (~50k donors) this drove peak memory from 48G to 256-350G per SLURM task on
  diverse BCR light chains, entirely from redundant `collect()`/`.to_list()` work repeated once
  per design instead of once. Added `prepare_fuzzy_features(cohort, key, candidates=, scope=) ->
  FeatureFrame` + `association(..., features=)` so the search is opt-in-cacheable: build it once,
  reuse across every design against the same cohort.
- Documented that `level_col`'s memory cost is **multiplicative, not additive**: `association()`'s
  feature-join duplicates every matched row once per design level (correct behaviour — each level
  needs its own incidence table — but easy to miss from the prior wording).

### Changed

- Bumped the `seqtree` floor to `>=0.6.1` — fixes a corrupted Miyazawa–Jernigan A–N contact
  energy in `structural()` (0.6.0) and names the offending sequence/index in `gapblock_matrix`'s
  alphabet error instead of just the bad symbol (0.6.1).

## 3.1.0 — 2026-07-28

### Added

- **The bundled `arda` model set is reachable.** `_bundled/arda/` has shipped 9 EM-refit models —
  the 7 human loci **plus mouse TRA/TRB** — in every wheel, but no public call could reach them:
  `SOURCES` listed only `("olga", "learned")`, and `load_bundled` keyed
  `_bundled/<source>/<LOCUS>` while the arda directories are `<organism>_<LOCUS>`, so
  `list_bundled()` reported them as absent. (`from_arda` is not a substitute — it returns the
  *placeholder* marginals meant to be refit by `infer_native`, not these refit ones.)

  `load_bundled` gains a keyword-only `organism=` (default `"human"`) and derives the directory key
  per set. This is the only bundled set in the **arda IMGT allele namespace** — the frame that
  arda-annotated pipelines such as `mirpy`'s prototypes and baked germline distances live in — and
  the only bundled set covering a non-human organism.

### Changed

- `load_bundled` now **raises** when a non-human `organism` is asked of the human-only `olga` /
  `learned` sets, instead of silently handing back the human model; `FileNotFoundError` lists the
  keys that *are* available for that set. Every existing caller passes `(locus, source)`
  positionally and is unaffected.

### Repository

- Consolidated the example notebooks: the old `notebooks/` directory was merged into
  **`examples/`**, so every marimo explorer now lives under `examples/` (docs / README / skills
  updated to match). Examples are not shipped in the wheel — this is a repository-layout change only.
- Merged the two aging notebooks into one **`examples/aging.py`**: it now covers the
  cohort-streaming stats (`diversity_cohort` / clone-size / spectratype), the coverage-standardized
  iNEXT diversity + rarefaction, and the pairwise-overlap→MDS divergence — the union of the old
  `aging.py` and `aging_airr_benchmark.py` (both removed, along with the now-unused
  `aging_manifest.json`). Fixes a latent `cdr3_aa`-vs-`junction_aa` key bug in the old benchmark
  notebook (stale since the v2.2.0 junction rename).
- **Every example is now a marimo notebook** — converted the two remaining plain scripts
  (`emerson_cmv_hla.py`, `scale_cohort.py`).
- **No user-specific absolute paths in the examples/tests** — each notebook resolves data from a
  gitignored **`./data_dump/`** directory first (symlink your copies there), then falls back to
  HuggingFace. VDJdb is fetched from the latest `antigenomics/vdjdb-db` release (cached to
  `./data_dump/`) instead of a hardcoded local checkout.

## 3.0.0

### Added — longitudinal clonotype dynamics (`vdjtools.dynamics`)

- **Recapture model** (`dynamics.capture`) — the VDJtrack size-bucket model: clonotypes binned
  singleton / doubleton / tripleton / large, Poisson capture probability `P = 1 − exp(−f·R)`,
  `Beta(captured, missing)` credible intervals, and the group-effect test — a log-linear
  `log(recapture) ~ size + group + log(div_ratio)` plus a per-bucket paired t-test across donors.
  Python port of Pavlova, Zvyagin & Shugay, *Front Immunol* 2024
  ([10.3389/fimmu.2024.1321603](https://doi.org/10.3389/fimmu.2024.1321603)).
- **Metaclonotype-grouped testing** (`dynamics.test_metaclonotypes`) — collapse a 1-Hamming
  (`scope="1,0,0,1"`) or 1-Levenshtein (`"1,1,1,1"`) CDR3 ball into one feature before the paired
  test, for power on convergent expansions.
- **edgeR NB-exact caller** (`dynamics.expansion_test`) — TMM normalization + qCML common
  dispersion + the negative-binomial exact test (as a Beta-Binomial conditional); the paper's §2.5
  complementary per-clone caller.

  (Complements the existing per-clonotype `dynamics.test_pair`, Ayestaran 2024.)

### Added — cohort-streaming summary statistics

- `io.map_samples(fn, items, *, workers=)` — thread-parallel per-sample reduce, `O(workers)`-sample
  peak memory, results in input order.
- `stats.diversity_cohort(cohort)` — the whole cohort's diversity table in one streamed
  count-spectrum pass, bit-exact vs the per-sample path.
- A `by=["sample_id"]` group-prefix on `spectratype` / `vj_spectratype` / `segment_usage` /
  `vj_usage` / `kmer_profile` / `v_kmer_c_profile` / `physchem_profile` — the whole cohort in one
  fused `group_by` over a `scan_cohort` LazyFrame.
- CLI `diversity` / `spectratype` / `segment-usage` / `overlap` gain `--threads N` (parallel over
  samples) and `--cohort DIR` (one streamed pass over a pre-ingested Parquet cohort); the `overlap`
  command now pre-aggregates each sample once.

### Added — CLI & packaging

- New `vdjtools` subcommands: **`convert`** (read any supported format — native / AIRR / Parquet /
  MiXcr / MiGec / MiTCR / immunoSEQ / IMGT / Vidjil / RTCR / TRUST4 / arda — and write the canonical
  table), **`downsample`**, **`filter`** (coding / non-coding / frequency / V-J segment), and
  **`pool`** (flat pool or incidence `--join`).
- Every command's `-o` is now **format-aware**: a `.parquet` / `.pq` path writes Parquet, anything
  else (or stdout) writes TSV.
- Development switched to **uv** — one repo-local `.venv`, no conda. `setup.sh` is rewritten to be
  uv-first (with a `python -m venv` fallback) and **portable across bash and zsh**. `environment.yml`
  is now optional, needed only for MMseqs2 (arda's aligner) + the slow arda round-trip tests.

### Added — notebooks (marimo, `[examples]` extra)

- `examples/vaccination_tracking.py` — clonotype tracking + recapture model across YFV / influenza
  / TBE vaccination time courses.
- `examples/aging.py` — cohort-streaming diversity / clone-size / spectratype across the Britanova
  ageing cohort.
- `examples/ankspond_motif.py` — the ankylosing-spondylitis TRBV9 "AS27" motif (disease vs HLA-B27
  carriage; Komech 2018).

### Fixed

- **C++ CI version drift** — the native `version()` is now single-sourced from `pyproject.toml`
  (parsed by CMake into the `VDJTOOLS_VERSION` compile definition), and both the C++ and Python
  version tests assert *agreement* rather than a hand-copied literal. A release bump can no longer
  redden CI (as the 2.9.0 bump did, leaving `tests/cpp/test_core.cpp` asserting `"2.8.0"`).
