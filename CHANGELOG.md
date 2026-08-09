# Changelog

Notable changes to vdjtools v2. Releases before 3.0.0 are recorded in the git tags
(`v2.5.0` … `v2.9.0`) and their commit history.

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
