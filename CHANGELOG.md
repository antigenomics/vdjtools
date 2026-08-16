# Changelog

Notable changes to vdjtools v2. Releases before 3.0.0 are recorded in the git tags
(`v2.5.0` … `v2.9.0`) and their commit history.

## 3.9.3 — 2026-08-16

### Changed — the release gate no longer runs the ten slowest tests

The publish gate added **~25 minutes of wall time to every release**, and the comment introducing
it in 3.9.2 asserted the opposite: *"running in parallel with the wheel builds, costs no extra wall
time."* Measured on the 3.9.2 release run (33m03s end to end), that was wrong by a wide margin:

| job | duration |
|---|--:|
| build-sdist | 0m16s |
| build-wheels macos | 2m42s |
| build-wheels ubuntu | 3m07s |
| build-wheels windows | 6m06s |
| **test (the gate)** | **31m18s** |
| publish | 0m43s |

Wheels finished at +6m06s; `publish` started at +32m19s. For comparison 3.9.1, before the gate
existed, took 15m26s. Inside the gate the split is **51 s** of install-and-compile against
**30m21s** of pytest — roughly 1:180, so the C++ build was never the problem.

The obvious fix does not work: `slow` is already deselected by `addopts`, so those 19 tests **never
run in CI at all** and the 30 minutes is entirely non-`slow` tests. Hence a second marker.

**`heavy`** — ten tests, measured at **318 s of a 322 s** three-file run:

| test | measured | why |
|---|--:|---|
| `test_gene_prior.py` (all 4) | 162.9 s | four full `infer_native` EM runs over the 89-allele bundled TRB locus; the C++ E-step is thread-parallel, so it degrades further on a 4-vCPU runner |
| `test_viterbi.py` `[TRB]` (4) + chosen-D | 130.2 s | the pure-Python reference scenario DP, documented at ~600× slower than native; the `[TRA]` twins cost 0.4 s because VJ has no D enumeration |
| `test_dynamics_paired.py::test_pvalues_are_calibrated…` | 25.4 s | a 200k/2M-read replicate pair plus a mandatory negative control |

`publish.yml` now runs `-m "not slow and not heavy"`. **Nothing loses coverage:** `ci.yml` still
runs the full suite on every push and PR at the same SHA — only the release path is trimmed.

Unlike `slow`, `heavy` is **not** in `addopts`, so a plain `pytest` locally still runs all ten.

Measured locally (M3, `HF_HUB_OFFLINE=1`):

| suite | result | wall |
|---|---|--:|
| full (`-m 'not slow'`, what the gate ran) | 1150 passed, 12 skipped, 19 deselected | 387.23 s |
| release path (`-m 'not slow and not heavy'`) | **1140 passed, 12 skipped, 29 deselected** | **72.65 s** |

**5.3× faster**, ten more tests deselected, same pass count otherwise. CI runs ≈4.7× slower than
this machine, so the projected release-gate job was **~6 min against the measured 31m18s** — at which
point it really would finish alongside the 6m06s Windows wheel build, as 3.9.2 wrongly claimed it
already did.

**Measured on this release run, and it beat the projection:**

| release | gate job | wall |
|---|---|--:|
| 3.9.2 | `Test before publishing` | 31m18s |
| 3.9.3 | `Test before publishing` | **4m04s** |

**7.7× on CI**, against 5.3× locally and the ~6 min projected. The runner gains more than the
arithmetic predicted because `test_gene_prior.py`'s thread-parallel C++ E-step degrades worst on a
4-vCPU runner, and that is exactly what `heavy` removes. The gate now finishes well inside the
6m06s Windows wheel build, so it genuinely costs no extra wall time — which is what 3.9.2's comment
asserted without measuring.

Coverage is unchanged: `ci.yml` ran the **full** suite green on this same SHA across ubuntu and
macOS × Python 3.10/3.12 (run 31970256603).


### Added — `Manifest.builder_version`

Which vdjtools built a model, set by `data.build_model` and persisted in `manifest.json`. A germline
defect lives in the **builder**, not the schema, so `model_version` could not answer *was this built
before or after the fix* — for the 3.9.1 J-anchor defect the answer had to be reconstructed by
comparing a shipped model's posterior against an unaffected reference fit. Now it is a field lookup.

Backward compatible: a manifest written earlier reads `""`, which means *predates 3.9.2*, not
*missing*. The bundled models are unchanged in this release and so all read `""`. A new
`test_collapse.py` check asserts the seven `learned` loci always report the **same** builder — a
*partial* regeneration is the dangerous state, because half the set would silently answer a
different question from the other half.

### Fixed — the PyPI upload was not gated on a green test suite

3.9.1 uploaded to PyPI while CI was still running its Test step. Same commit, so the code was
covered, but nothing enforced the **order**, and a CI-only failure could have landed on an
already-published version. `needs:` cannot reference another workflow, so `publish.yml` now carries
its own `test` job — same extras as `ci.yml`, so the OLGA oracle suite actually runs instead of
silently skipping — and `publish` needs it. It runs alongside the wheel builds, so it costs no
extra wall time.

### Note — the bundled models are NOT regenerated, deliberately

3.9.1 fixed `from_olga(derive_orf=True)`, which had reconstructed J CDR3-region germlines as
`full[anchor:]` (the **V** convention) instead of `full[:anchor + 3]`. Fixing the builder does not
fix already-built parquet, so the bundled `learned` human TRA model still carries 11 J germlines
from the wrong side of the anchor, and `TRAJ35` — the one functional allele among them — reads
**~774× low** against the unaffected `arda` fit.

That is left in place on purpose. **0** of VDJdb's 30,937 human TRA records (27,272 unique
junctions, 49 TRAJ genes) use any of the 11 affected alleles — verified as real absence, not a
name-matching artefact, since neighbouring `TRAJ34` (891 records) and `TRAJ36` (494) are well
covered. Regenerating the set is a multi-hour all-loci job, and there is no measured consumer of
the difference. The affected alleles stay pinned by name in `test_collapse.py`.

**If you use human TRA and care about `TRAJ35` usage, use `load_bundled("TRA", "arda")`** — that set
is built from arda germline, where `derive_orf` never runs, so it cannot carry this defect.

Measurements: `bench/results/vdjtools_germline_pgen_shift.md` in the benchmark repo, which also
records what the **3.9.1** collapse fix already moved — 6,676 of 41,322 VDJdb human TRB junctions
(16.2%) went from `Pgen` exactly `0.0` to positive, all of them `TRBJ2-7`, with max
`|Δlog10| = 0.000000` across every junction that was already non-zero.

## 3.9.1 — 2026-08-16

### Fixed — a collapsed gene could be represented by a **non-functional** allele, making Pgen a silent zero

`collapse_alleles` picks one allele per gene, keeps its germline and relabels it `gene*01`. It
ranked candidates by CDR3-region germline length, then usage, then the allele **name** — but usage
was only ever passed for **V**. For J and D the key therefore degenerated to the name after the
length tie, and `max` took the lexicographically last allele, which was then relabelled `*01`.

Human `TRBJ2-7` is the worst case. `*01` is functional and templates `SYEQYF`; `*02` is an IMGT
**ORF** and templates `SYEQYV`. Both are 19 nt, so length did not separate them and `*02` won on
the name — the collapsed model shipped `*02`'s germline under `*01`'s label:

```
Pgen("CASSIRSSYEQYF" | TRBV19*01, TRBJ2-7*01), bundled human TRB `olga` model
             collapse=True   collapse=False
  3.9.0        0.0             4.7889e-08
  3.9.1        4.9986e-08      4.7889e-08
```

**Exactly zero, with no error raised.** Measured on 864 `TRBJ2-7` nucleotide junctions drawn from
the model itself (`generate(seed=7)`, n=4000): **864 of 864 scored `pgen_nt == 0`** against the
3.9.0 collapsed model, 3 of 864 against 3.9.1 — and those 3 are genuine `*02`-germline draws, which
a single-germline collapsed model cannot represent by construction. Summed over the 864: `0.0` at
3.9.0, `3.357e-07` at 3.9.1, against `3.725e-07` uncollapsed (mean `|Δlog10|` 0.157). This affects
**every released version that shipped these bundled models**, on the default `collapse=True` path.
`load_bundled(..., collapse=False)` reproduces the old numbers for anyone who needs them.

The representative is now ranked by **length → IMGT functionality (`F` > `ORF` > `P`, read from
arda's `cdr3_anchors.tsv`) → usage → prefer `*01` → name**, and the J and D usage marginals are
actually passed. Length still leads: over all 23 bundled models the two orders differ on exactly one
gene, human `TRBV23/OR9-2`, where the only non-pseudogene allele has an **empty** CDR3 germline and
leading with functionality would install it — trading one silent zero for another.

**65 genes change their collapsed germline** (every change length-preserving, 63 of them onto
`*01`): 10 in `olga`/human, 9 in `learned`/human, 29 in `arda`/human, 17 in `arda`/mouse. Beyond
`TRBJ2-7`: `IGKJ2`, `IGKJ4`, `TRAJ47` in all three sources; `TRAJ24/32/37/41`, `IGHJ4`, `IGHJ6`,
`IGKV1-39`, `IGKV3D-20`, `IGHV2-70`, `IGLV1-41` and 7 `IGHD` genes in `arda`; 12 mouse `TRAV` plus
`TRBJ1-1`/`TRBJ1-5`. **Any Pgen, generation or scoring output on a collapsed model can move.**

### Fixed — `from_olga(derive_orf=True)` rebuilt J germlines from the wrong side of the anchor

The CDR3 region lies on opposite sides of the conserved-codon anchor per segment: V runs
`full[anchor:]`, J runs `full[:anchor + 3]`. `_genomic_table` used the V slice for both, so every
ORF/P J allele it reconstructed got the framework *downstream* of Phe118 — still a plausible
in-locus sequence, hence silent. 11 alleles in the bundled `learned` human TRA model carry it
(`TRAJ1*01`, `TRAJ2*01`, `TRAJ19*01`, `TRAJ25*01`, `TRAJ35*01`, `TRAJ51*01`, `TRAJ55*01`,
`TRAJ58*01`, `TRAJ59*01`, `TRAJ60*01`, `TRAJ61*01`); the builder is fixed, but clearing the shipped
parquet needs a model rebuild, so the new test pins those 11 as a named exception.

### Added — permanent anchor tests over every shipped model

`tests/python/test_collapse.py` now asserts, across all 23 bundled models (`olga`/`learned`/`arda`
× human, `arda` × mouse TRA/TRB) and both `collapse=True` and `False`, against arda's per-allele
`functionality`/`status`/`templated_aa` rather than a guessed convention:

- every functional J germline translates, in its anchor frame, to a terminal **F or W**. A
  collapsed row is held to its *gene's* standard, not to whichever allele supplied the germline —
  judging it by that allele would have exempted `TRBJ2-7*02` as "an ORF" and waved the defect
  through. Named exceptions: `TRAJ35*01` (arda records `templated_aa=IGFGNVLHC` at `status=ok`;
  IMGT still calls it F) and `IGHJ6*02` in OLGA's namespace only (OLGA ships it 1 nt short of the
  Trp118 codon; arda's own copy templates `YYYYYGMDVW` and passes);
- no gene is represented by an ORF/pseudogene allele where a functional one exists — 18 genes
  failed this before the fix;
- the CDR3-region germline sits on the documented side of the anchor;
- regression pin: `Pgen("CASSIRSSYEQYF" | TRBJ2-7*01) > 0` in all three sources, with collapsed and
  uncollapsed agreeing to within 0.3 in `log10`.

## 3.9.0 — 2026-08-16

### Added — `Pgen` of a degenerate motif

`native.pgen_aa_degenerate(model, allowed, v=None, j=None)` and `pgen_aa_degenerate_batch` expose
the masked transfer-matrix DP that already backed `pgen_aa` and `pgen_aa_hamming1`. `allowed` is one
entry per position, each a string of permitted residues; `""` or `"X"` means any residue.

This makes the total generation probability of a V/J/length-pinned motif — a VDJdb cluster PWM, say
— a single exact call, with no enumeration and no inclusion–exclusion. Motivated by epitope
precursor-frequency estimation (`appendix/pgen_motif.md`).

Note `pgen_aa` itself still scores an `X` as **0.0**, because `mask_for_aa` matches the genetic code
by exact character. Use `pgen_aa_degenerate` when a position is meant to be a wildcard.

## 3.8.0 — 2026-08-15

Single-cell interop: vdjtools now sits inside the downstream single-cell ecosystem instead of
ending at its own frame.

### Fixed — `paired_pgen` returned nothing but nulls on real CellRanger data

The bug that mattered most here, and it was silent. CellRanger reports **gene**-level V/J calls
(`TRBV10-3`); the model is keyed by **allele**. `native.pgen_aa` raises on a gene name on
purpose — the old `-1` fallback meant *marginalise over every allele* and once returned a Pgen
**2.38x too high** with no error — but `sc.pgen._chain_pgen` caught that with a bare
`except Exception: return None`. Net effect: on the single most common real input, `pgen_alpha`,
`pgen_beta` and `pgen_paired` were **100% null**, with nothing to indicate why. Measured on the
public dCODE donor-4 run: **27,268 of 27,268 receptors null**.

`paired_pgen` now resolves a gene to its representative allele (`*01` where the model has it)
before scoring — deliberately and documented, *not* by falling back to marginalising. Same
dataset: **24,325 of 27,268 now scored** (median paired Pgen 2.1e-19). `resolve_genes=False`
restores exact-allele-only matching, and an all-null locus now emits a `UserWarning` instead of
shipping a silent column. The `except` is narrowed to `(KeyError, ValueError)` with a
non-`str` junction guard, so unrelated failures stop being swallowed.

### Fixed — a barcoded AIRR table was silently collapsed into a bulk repertoire

`io.sniff_format` had no `cell_id` branch, so CellRanger's `airr_rearrangement.tsv` (or any
barcoded AIRR table) sniffed as `"airr"`/`"arda"` and `read_airr` pooled reads **across cells**,
dropping the barcode with no error. It now sniffs as `"airr_cell"` and `io.read` refuses it,
naming `sc.read_airr_cell` instead; `fmt="airr"` still pools on purpose.

### Added — one interchange format, four ecosystems

scirpy, dandelion and scRepertoire all read the same thing: a flat AIRR Rearrangement table with
`sequence_id` + `cell_id`. So `vdjtools/sc/airr.py` is one emitter (`to_airr`) and one inverse
(`from_airr`), and each bridge is a thin adapter — `write_airr`, plus `write_screpertoire`
(`format="airr"|"10x"`). It reconciles the two spellings that otherwise bite: AIRR says
`junction` where vdjtools says `junction_nt`, and scRepertoire's parser reads `consensus_count`
where scirpy and dandelion prefer `umi_count`, so both are emitted.

- **scirpy / scverse** — `to_scirpy` (scirpy's `obsm["airr"]` awkward layout, `index_chains` run
  by default; `gex=` gives a `MuData`) and `from_scirpy`. Writing **delegates** to
  `scirpy.io.read_airr` so no copy of their schema can drift here; reading is ours and needs only
  `awkward`, so consuming someone else's AnnData costs no scirpy install.
- **dandelion** — `to_dandelion` / `from_dandelion`, plus `read_h5ddl`: `.h5ddl` is plain HDF5, so
  a dandelion result opens with `h5py` alone.
- **`push_obs`** — attach vdjtools-computed columns (`pgen_paired`, mispairing flags) to an
  `AnnData.obs` or `Dandelion.metadata` you did not build. Refuses a multi-pair frame rather than
  silently picking one row per cell.

### Added — ingestion

`read_10x` now accepts `filtered_contig_annotations.csv` as well as `all_contig_annotations.csv`
(one CellRanger writer, one layout) and tolerates version drift — `fwr*`/`cdr1`/`cdr2` are CR6+,
`exact_subclonotype_id` CR4+, `sample` only under `cellranger multi`, and `raw_consensus_id` is
used when present rather than required. `read_arda_cells` reads `arda cells` output
(`.contigs.airr.tsv` + `.chains.tsv`), surfacing arda's own per-chain verdict as `arda_status`
**without acting on it** — arda's call and `resolve_chains`' call are independent answers to the
same question. `productive` joins `SC_COLUMNS` so the emitted AIRR table is schema-valid.

### Added — CLI, docs, example

`vdjtools sc` — `convert`, `pair`, `qc`, `pgen`, and
`export --to airr|scirpy|dandelion|screpertoire|screpertoire-10x|airr-cell`, each exposing the
matching library options: `--fmt`, `--require-cell`, `--require-high-conf`, `--consensus`,
`--locus-pair`, `--resolve`, `--flag-mispairing`, `--max-slaves-per-master`, `--drop-mispaired`,
`--source`, `--condition-vj`, `--resolve-genes`, `--alpha-locus`, `--beta-locus`,
`--index-chains`, `--repertoire-id`. (GEX pairing stays library-only --
`to_scirpy(cells, gex=...)` returns the MuData; a CLI flag for it only added a MuData
serialisation step, which failed in CI environments we could not reproduce.) The input format is sniffed from the
**header**, not the filename — a renamed export still works and a bulk table is refused by name
rather than mis-parsed. `sc pgen` reports `scored N/M receptors`, so a naming mismatch is a
number on screen rather than a column of nulls to notice later. A dedicated
`docs/singlecell.rst` (the `usage.rst` section is now a pointer), and
`examples/single_cell_interop.py`, a marimo notebook running the whole path on dCODE donor 4 —
which is what surfaced the Pgen bug above.

`[sc]` gains `awkward` + `mudata`; a new **test-only** `[interop]` extra carries `scirpy` and
`sc-dandelion` (PyPI name; imports as `dandelion`). CI installs it best-effort and reports
whether the round-trip tests ran or skipped, since their dep chains break on new matplotlib.
The format contract itself (`test_sc_airr.py`) has no optional deps and never skips.

## 3.7.3 — 2026-08-15

Housekeeping. The first PyPI release since 3.7.0, so it carries 3.7.1 and 3.7.2 with it.

### Fixed — the iNEXT bootstrap carried a fallback that could never be taken

`stats/inext.py` guarded its `_core` import in a `try/except` and dispatched between the native
bootstrap and the numpy reference at call time. `_core` is a build-time dependency — an install
without it does not exist — so the `except` branch and `_bootstrap_se_dispatch` were dead, and
`inext_batch` raised its own "requires the native _core extension" for a state that cannot occur.
Both are gone; the import is now local to the two functions that need it, which keeps
`import vdjtools.stats` as light as the guard made it. The numpy `_bootstrap_se` stays, unchanged
— it is the reference the tests compare the native kernel against.

### Removed — three dev-notes files the changelog had already absorbed

`NOTES.md`, `ROADMAP.md` and `SUGGESTED_EDITS.md` recorded the phase narrative from before the
changelog existed, and had been drifting from it since. `CHANGELOG.md` is the release-by-release
record; `CLAUDE.md`'s "Open loops" is what is in flight. `CLAUDE.md`, `README.md`,
`docs/index.rst` and the sdist exclude list no longer point at the deleted files.

## 3.7.2 — 2026-08-14

Documentation accuracy. No code change.

### Fixed — the preset-ranking corpus was described as larger than it is

`vdjtools.signature.presets` and `docs/signature.rst` both said the rankings come from "several
hundred study groups, tens of thousands of samples". Counted, the sweep panel is **182 study
groups over 198 accessions, 14,553 samples** — which the same page already stated correctly two
sections earlier ("14,553 samples × 1,369 columns, 182 studies"), so the file disagreed with
itself. Both places now carry the counted figure. The accession list is published in the analysis
repo's `heldout/signature_studies.tsv`, so a reader can check the claim rather than take it.

## 3.7.1 — 2026-08-14

Audit pass. Two fixes, both cases where a feature was reachable from Python and not from the
command line that ships it.

### Fixed — `keep=` stopped at the readers, so the SHM block could never be computed from the CLI

3.7.0 added `keep=` to `read_airr` / `read_vdjtools` / `read_parquet` so `v_identity` — the one
field the signature needs that the canonical eight columns do not carry — could reach
`vsig:shm:IGH:mean_v_identity`. The dispatcher `io.read` and the batch mapper `io.map_samples`
did not take the argument, and those are what the CLI uses: `vdjtools signature` on a file
carrying `v_identity` reported `mask:IGH:shm = 0` and `mean_v_identity = nan`. Both now take
`keep=`, and the `signature` command passes `("v_identity",)`.

The column now populates on any input that has the field. Nothing else changes: `keep=()` is the
default everywhere, and the legacy converters, which narrow to the canonical schema, ignore it.

### Fixed — a coverage level of exactly 1.0 warned its way to the right answer

`mir.signature` passes `cstar = 1.0` deliberately, as an "unreachable" sentinel, for a locus
where no coverage level could be established — the diversity block is then supposed to fail its
own estimability check and mask out. It did, but `_invert_coverage` got there by evaluating
`log(1 - 1.0)` and doing inf arithmetic, emitting three `RuntimeWarning: divide by zero` per
call into the user's terminal. It now returns `inf` directly. Same `m`, same method, same mask —
without the noise.

## 3.7.0 — 2026-08-14

### Added — `vdjtools signature` on the CLI, with the help text as the primary documentation

The command a collaborator actually runs. No Python:

```bash
vdjtools signature --preset classify -m metadata.txt --base-dir samples/ -o sig.tsv
vdjtools signature --preset compact a.tsv b.tsv.gz -o sig.tsv
vdjtools signature --preset classify --describe     # the columns, reading no input
vdjtools presets                                    # the named feature sets, ranked
```

`--help` on both commands carries worked examples, the three `recommended` presets and when each
applies, the pointer to `mir signature` for the geometry half, and the CDR3-vs-junction trap
(a file carrying only IMGT `cdr3_aa` is two residues short everywhere, which shifts the length,
k-mer and Pgen features). `docs/signature.rst` opens with the same quickstart, and the README and
`examples/README.md` lead with `--preset classify` rather than a `specific`-ranked set.

Because that help text is written for a terminal — indented example blocks, which are not valid
reStructuredText — `signature` and `presets` are excluded from the `vdjtools.cli` autodoc, with a
note on the API page saying where to read them instead.

### Fixed — 221 signature tests were invisible to CI

The seven new test files landed in `tests/` while `testpaths = ["tests/python"]`, so a plain
`pytest` collected 789 of 1010 and the CI job (`pytest tests/python -q`) never ran one of them.
Moved into `tests/python/` with the rest; default collection is 1004 passed / 6 skipped.

### Added — `vdjtools.signature`: VSIG, the statistics half of a portable repertoire signature

One repertoire in, a **fixed, named, positional** feature vector out — the object you hand a
collaborator so their matrix and yours are the same coordinate system. The geometry half lives in
`mir.signature`; the shared column contract lives here, because mirpy depends on vdjtools and not
the reverse, and two copies of a contract are not a contract.

```python
from vdjtools.signature import vsig, vsig_cohort, columns, describe
v = vsig({"TRB": df}, tier="standard")
describe("standard")            # column, sig, block, locus, feature, tier, transform, flags
```

Four modules. `layout` is the contract — loci, the `core ⊂ standard ⊂ full` tiers as exact
**index subsets** of one frozen column order, and a per-feature (not per-block) transform
declaration, because a clonality block legitimately mixes a CLR-transformed composition with a
logit-transformed proportion. `transform` is the variance-stabilising layer. `blocks` computes.
`assemble` puts them in order.

Every transform choice is **denominator-aware**, because the alternative silently lies about
shallow samples: Haldane–Anscombe `logit` so `0/3` and `0/500` are different numbers, Anscombe
`arcsine` so a share is defined at exactly zero, and `clr` over the *whole* composition before any
coordinate is selected — shipping *k−1* parts, since all *k* are linearly dependent and would put
a guaranteed zero eigenvalue in any PCA.

Diversity is compared at a **frozen coverage level**, and `estimable()` **refuses** rather than
extrapolates. Real repertoires attain Good–Turing coverage 0.24–0.58, so a textbook `C* = 0.95`
puts every sample into extrapolation, where the same statistic inflates roughly tenfold. A hole a
model can see beats a confident wrong number. For the same reason `clonality` is rebuilt from the
coverage-standardised Hill numbers, `1 − ln(¹D)/ln(⁰D)`: the observed Pielou evenness it replaced
drifted 0.510 over a 667× depth range, against 0.023 for the standardised form.

### Fixed — the CLR zero replacement could consume the composition it was correcting

The textbook multiplicative replacement puts `delta = 0.5/m` on each zero part and scales the rest
by `1 − n_zero·delta`. On a *shallow* composition that is bigger than the composition: three
parts, one observed, `m = 1` gives two replacements of 0.5 and scales the one real part to exactly
zero, whose log is `-inf` — a value that then propagates through every downstream reduction. Found
while emitting a real 4,000-sample corpus. The replaced mass is now capped below half the smallest
*observed* part, which is the only property a replacement needs; the cap is inactive whenever `m`
exceeds the number of parts, i.e. everywhere outside that tail.

### Fixed — `pgen_block` reloaded the recombination model on every call

Loading and collapsing a bundled model costs 0.4–1.8 s; the Pgen batch that follows costs ~0.15 s.
A corpus emission therefore spent 80–95% of its time re-reading seven files it had already read,
and a seven-locus sample paid it seven times over. Memoised per locus: **1.54 s → 0.01 s** on the
second call.

Also worth knowing when emitting a corpus: `pgen_block`/`vsig` default to `threads=0`, meaning
*all cores*, which inside your own process pool means every worker claims the whole machine. On a
16-core box, 14 workers took the load average to 227. Pass `threads=1` there.

## 3.6.1 — 2026-08-14

Audit pass. No library behaviour changes.

- **`examples/emerson_cmv_hla.py` did not parse on Python 3.10/3.11.** A backslash inside an
  f-string *expression* (`f"…{meta['hla'].str.contains(r'HLA-A\*02').sum()}…"`) is 3.12-only syntax,
  but the package declares `requires-python = ">=3.10"`. The regex is hoisted to a local.
- Unused imports and multi-import lines cleaned out of `examples/` (`ruff --fix`).
- `[tool.ruff.lint]` now ignores `E702`/`E741`/`E731` — the paired-short-statement style and `l`/`O`
  loop scalars are deliberate throughout the examples and OLGA oracle shims. `ruff check .` is
  green, so a real finding is visible again instead of being buried in 30 style hits.
- Repo cleanup: 1.7 GB of regenerable artifacts removed (`examples/.data` notebook caches,
  `docs/_build`, `examples/__marimo__`, tool caches, `.DS_Store`), plus four worktrees whose
  branches were already fully merged into `master` (`feature/biomarker-cooccurrence`,
  `feature/dynamics`, `chore/vdjmatch-pin`, `claude/trusting-torvalds-65e7b7`). `feature/cdr3-viterbi`
  and `signature` still carry unmerged commits and were left in place.

Verified at this commit: `pytest tests/python` 783 passed / 6 skipped; `sphinx-build -W` clean.

## 3.6.0 — 2026-08-12

### Added — `infer_nt`: the nucleotide CDR3 behind an amino-acid one

A VDJdb record carries `(V, J, CDR3aa)` and no nucleotides, so none of the boundary markup a
repertoire analysis wants is there. `infer_nt` reconstructs all of it:

```python
from vdjtools.model import infer_nt
sc = infer_nt(model, "CASSLGQAYEQYF", v="TRBV5-1*01", j="TRBJ2-3*01")
sc.cdr3_nt, sc.v_end, sc.d_call, sc.d_start, sc.d_end, sc.j_start, sc.pgen, sc.margin
```

Two stages. A codon-constrained **max-product DP** over every scenario `pgen_aa` sums — germline
positions pinned to their segment, each free N-region position taking the nucleotide that maximises
`P(nt₁)·∏P(nt_k | nt_{k−1})` under the VD/DJ/VJ dinucleotide model; then a `pgen_nt` re-score of the
survivors, because stage 1 maximises the *joint* `P(nt, scenario)` while the contract is about the
*marginal* `P(nt)`. `pgen` on the result is a real `pgen_nt`.

Stage 1 is **native**: the same Murugan/OLGA `Pi_L·Pi_R` transfer matrix `pgen_aa` already uses,
with `max` in place of the sums and the winning `(V, delV)` / `(J, delJ)` carried through the state
(`native.best_aa_scenarios`). It returns *scenarios*, not nucleotide paths — once the scenario is
fixed, recovering the nt string is one cheap DP over the two insertion blocks, so the sweep carries
no back-pointers. A first cut enumerated the scenarios in Python instead and cost 1.7 s per TRB
CDR3; that implementation is kept as the reference the native one is tested against.

Per-sequence cost, 25 generated productive draws per locus. Leaving the calls out is nearly free,
because the DP sweeps V and J either way:

| locus | ms/seq, V/J known | ms/seq, V/J free |
|---|---|---|
| IGK | **0.41** | **0.48** |
| IGL | 0.46 | 0.50 |
| TRG | 0.78 | 0.68 |
| TRA | 0.86 | 2.43 |
| TRB | 3.02 | 3.20 |
| TRD | 7.07 | 7.78 |
| IGH | 87.14 | 85.88 |

So all 80k VDJdb records — TRA and TRB — take **about 2-4 minutes**. IGH is the outlier at 87 ms:
it carries 60+ D alleles and the D placement loop scales with that library, so plan minutes per
10k there rather than per 80k.

**Measured against the exponential oracle** (generated productive draws, V/J pinned, 10 residues):

| variant | TRG | TRA | ms/seq |
|---|---|---|---|
| one scenario, best codon per residue | 9/25 | 4/19 | 0.06 |
| stage 1 only (`keep=1, n_best=1`) | 21/25 | 15/19 | 0.3 |
| stage 1 + marginal re-score (defaults) | **25/25** | **19/19** | 1.2 |

The cheap shortcut fails because a trim chosen before the codons pins a codon the true optimum
would have trimmed away — the same unsoundness as pinning the germline flanks.

**Three call-input modes**, because annotation tables have all three: one allele (the normal mode);
several, as the comma-separated string an ambiguous AIRR `v_call` carries or as a list; or nothing
at all, where the DP marginalizes over every gene — 0.26 ms against 0.23 ms with both pinned, so
unknown calls cost essentially nothing.

Tandem-D is not enumerated in stage 1 (a single D trimmed to zero length already reaches every
middle, so D-D can only reorder candidates, not add them); the stage-2 `pgen_nt` counts it in full.

### Fixed — EM could relearn a genomically impossible D–J pair

A D can only recombine with a J lying 3′ of it, and in TRB the clusters interleave
(TRBD1·TRBJ1·TRBD2·TRBJ2), so `P(TRBD2 | TRBJ1-*)` must be zero. The learned TRB model had it at
**0.0909**, and OLGA's own TRB gives `P(TRBD2*01 | TRBJ1-6*01) = 0.333`.

The constraint is now applied in **both M-steps before normalization** — a post-hoc patch would be
undone by the very next iteration. `reference.forbidden_dj_pairs` derives the forbidden set from
IMGT cluster numbering, `infer.enforce_dj_order` repairs an existing model, and `check_model` gains
an `impossible_dj_pair` check (`warn` for faithful OLGA imports, which must stay byte-exact for the
Pgen invariant; `error` otherwise).

The bundled `learned` TRB is rebuilt with the constraint: **9 iterations instead of 11**, final
log-likelihood **−33.7575** against −33.7604 — the constraint *improves* the fit.

### Changed

- No emoji anywhere in the repository; the `⛔`/`⚠` markers are now `WARNING:` / `NOTE:`.

## 3.5.0 — 2026-08-12

### Fixed — `generate(seed=)` was not reproducible across processes

`collapse_alleles` used unordered polars `group_by`, so a collapsed model's table row order varied
between processes and the *same* seed drew a different allele. `collapse=True` is the default load
path, so this affected ordinary use. The collapse now maintains order throughout. Expectations
recorded from `generate()` before this release are stale; the bundled models are unaffected, since
they ship uncollapsed and are collapsed at load.


### Added — `vdjtools.model.viterbi`: the argmax side of the Pgen DP

`pgen_nt` **sums** over every recombination that could produce a nucleotide CDR3.
**`best_scenario(model, cdr3_nt, v=, j=)`** takes the **maximum** over the same loops and returns
the single most likely one — which *is* the V/D/J boundary markup:

```python
from vdjtools.model import best_scenario
sc = best_scenario(model, cdr3_nt, v=v_call, j=j_call)
sc.v_end, sc.d_call, sc.d_start, sc.d_end, sc.j_start   # 0-based, half-open, CDR3-nt space
```

It re-derives nothing: every probability comes from `prepare()`'s tables, and the D placement is a
max-product mirror of `pgen._d_middle` over the same `P(D|J)·P(delD|D)·Pins(VD)·Pins(DJ)` terms.

WARNING: **The D therefore obeys `P(D|J)`.** TRBD2 lies 3′ of the whole TRBJ1 cluster, so deletional
joining can never produce a TRBD2–TRBJ1 pair and the model encodes that as a zero. An earlier draft
here chose D by longest exact substring and ignored `j` entirely — it would have called the
impossible pair. There is a regression test.

Validated on 200–500 generated draws per locus: the V span **is** the V germline, the J span the J
germline, the D span the D germline, `scenario_p ≤ pgen_nt` (a maximum cannot exceed the sum it is
taken over), and `scenario_p` recomputes exactly from the reported path's own table entries.

`infer_nt_bruteforce` is an exact but exponential **oracle for tests**.

NOTE: Both functions assume an **in-frame** CDR3 (`len(nt) == 3 × len(aa)`). Real productive receptors
satisfy this; an out-of-frame draw does not (measured: 8 aa against 25 nt).

### Changed

- `arda-mapper` pinned to **>= 2.19.0** (was >= 2.5.5).

### Fixed — `generate(model, n, seed=)` was not reproducible across processes

Same seed, same wheel, a **new interpreter → a different draw**. Within one process it looked
perfect, which is why no existing test caught it.

Root cause: **`collapse_alleles`** — which `load_bundled(..., collapse=True)` runs by default —
built its tables with unordered polars `group_by().agg()`. `group_by` is a multithreaded hash
aggregation, so the collapsed table's **row order varied per process**; `_cum` then assigned the
same cumulative interval to a different allele, and the same `rng.random()` drew a different one.
`_pick` and `default_rng` were correct throughout — the ordering beneath them was not.

NOTE: **Not hash randomisation.** `PYTHONHASHSEED=0` did not help, which is what ruled it out and
pointed at the aggregation. Same class as the nondeterminism recorded against arda's `correct`
stage.

Fixed by `maintain_order=True` on all 12 `group_by` calls in `collapse.py` and all 7 in
`generate.py`. Verified identical across 5 separate processes on TRA, TRB and IGH.

WARNING: **This changes generated output** for a given seed — it has to, since the old order was
arbitrary. Any recorded expectation from `generate()` predating 3.3.0 must be re-derived.

New tests run the sampler in a **subprocess**, because an in-process test agrees even with the bug
present; they fail 4 of 5 without the fix.

## 3.4.0 — 2026-08-12

Follows 3.3.0's model workshop with the defects that workshop then found, and retrains every
bundled model.

### Added

- **Live EM progress, checkpointing and exact resume.** `infer`/`infer_native` take
  `progress=callable(iter, loglik, rel_change, n_scoreable)` — `infer.print_progress()` is a
  ready-made one — so a long fit is visibly converging rather than merely running; the relative
  change it reports is exactly what is compared against `tol`. They also take
  `checkpoint=DIR`/`checkpoint_every=N`, saving the model after each iteration (written to a
  sibling directory and swapped in, so a kill mid-write leaves the previous checkpoint loadable),
  and `infer.resume(DIR, seqs)` continues from one. **Resuming is exact**: 3 iterations plus a
  resumed 4 give the same log-likelihood *and* the same tables as an uninterrupted 7, and the
  training log spans every attempt. Exposed as
  `vdjtools model learn -v --checkpoint DIR --resume DIR`; `vdjtools model build -v` additionally
  stops swallowing arda's mapping output. This matters because IGH's EM enumerates ~1,225 D pairs
  per read against TRB's 9 — 12 minutes per iteration on a 112-core node, and more than 110 minutes
  without finishing one on a laptop.

### Changed

- **All seven bundled `learned` models retrained** on the full non-functional read corpus (every
  available read; out-of-frame *and* stop-codon, since both escaped selection). Every one converged
  on tolerance with a step-by-step monotone log-likelihood, and `check_model` reports zero errors:

  | locus | clonotypes | iters | log-likelihood |
  |---|---|---|---|
  | TRA | 34,238 | 10 | −21.52 → −20.01 |
  | TRB | 122,703 | 11 | −37.01 → −33.76 |
  | TRG | 14,305 | 8 | −21.90 → −20.33 |
  | TRD | 10,915 | 7 | −40.94 → −35.79 |
  | IGH | 141,607 | 9 | −57.16 → −49.88 |
  | IGK | 256,347 | 6 | −16.18 → −14.28 |
  | IGL | 23,469 | 12 | −6.63 → −6.00 |

  Each now ships its training log, so a model states what it was fitted on. IGH was trained on a
  112-core cluster node; the rest fit comfortably on a laptop.
- **Ambiguous junction bases are substituted, not dropped.** They previously crashed EM with a bare
  `KeyError` from inside the native encoder. Both training entry points now substitute `A` by
  default and warn with the count (`ambiguous=None` drops instead) — it affects ~0.01% of these
  reads, so dropping cost sample size for nothing. It is a substitution, not a marginalization.

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
  — precisely how the CLI rename in 3.3.0 stayed invisible. It now raises with arda's stderr and the
  installed `arda-mapper` version.
- `reference.read_fasta` wraps arda's FASTA parser with gzip support; arda's opens with plain
  `open()`, so a `.gz` reached it as mojibake and died on the first byte.

### Documented

- **Known quirks of the OLGA models** (`docs/model.rst`). The bundled `olga` set is a bit-faithful
  import, and that includes its defects. Verified against OLGA's raw `model_marginals.txt` with
  OLGA's own parser: deletion mass on unreachable trims is **OLGA's**, with identical fractions to
  4 dp (`IGHV4-30-4*01` 100% — its Pgen is identically zero in OLGA too — `IGKJ4*02` 80.9%,
  `TRAV20*03` 54.7%), and our Pgen matches olga-pip exactly (ratio 1.000000) on every sequence OLGA
  will score. Correcting it would break the exact-OLGA-Pgen invariant, so `check_model` reports it
  at `warn` rather than `error` for an OLGA-sourced model. Also covers the empty-germline ORF genes,
  protocol-specific V/J usage, and OLGA's refusal of out-of-frame input.

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
