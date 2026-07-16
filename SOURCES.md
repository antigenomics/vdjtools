# SOURCES

Provenance of every dataset used or produced by vdjtools. Never guess a source — record it here.

## Germline reference (canonical)

**arda's germline library is the single source of germline truth** — every V/D/J germline
sequence and CDR3 anchor resolves from it by allele name (`vdjtools.model.reference.load_germline`),
so annotation ↔ scenarios ↔ stitching ↔ Pgen share one coordinate frame. arda's anchor
convention is byte-identical to OLGA's (0-based Cys104/[FW]118 offset into the full germline).

| Dataset | Origin | Format | Notes |
|---|---|---|---|
| V/J germline + CDR3 anchors | arda `database/vdj/<org>/cdr3_anchors.tsv` (via `arda.cdr3fix.load_anchors`) | TSV, per-allele | CDR3-region germline + anchor + functionality; **full-length V/J germline not shipped** (build-time only) — a P1c/stitching prerequisite |
| D germline | arda `database/vdj/<org>/d_germlines.fasta` | FASTA `>LOCUS\|allele` | full D germline, no anchor |

**OLGA/IGoR model files contribute recombination probabilities only** (and, for bootstrap
models, their own IMGT-vintage germline is kept for exact-Pgen fidelity — see the P1a note).

## Phase 1 — model engine (bootstrap)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| OLGA default models (7 human loci + mouse) | `mirpy/mir/resources/olga/default_models/{human_T_alpha,beta,gamma,delta, human_B_heavy,kappa,lambda}/` (reference only — do not copy mirpy code) | `model_params.txt`, `model_marginals.txt`, `V/J_gene_CDR3_anchors.csv` | read from the local mirpy checkout | published OLGA models (Sethna et al.); **derived** generative-model parameters |
| Synthetic out-of-frame training seqs (~100k/chain) | OLGA generation via mirpy tooling; sample at `mirpy/tests/assets/olga_humanTRB_1000.txt.gz` | TSV (V call, J call, CDR3 nt) | generate with OLGA from the models above | **computed** (Monte-Carlo draws); **no tandem D by design** — see plan Phase 1 note |
| Real AIRR reads (model training, all 7 loci) | HF dataset [`isalgo/airr_model_read`](https://huggingface.co/datasets/isalgo/airr_model_read) (owner, cc-by-nc-nd) | raw 5'-RACE **FASTQ** `{group}/{CHAIN}.{functional,nonfunctional}.fq.gz`; groups = `human`, `human_fetal` (TdT-low), `mouse` | `vdjtools.model.data.prepare(group, chain, label)` = fetch (lazy `huggingface_hub`) → `arda rnaseq map` (V/D/**D2**/J + junction + productivity) → `unique_clonotypes` (dedup to `(v,j,junction)`) | **experimental**; out-of-frame reads train the models (EM), functional reads are the selection test set. D-D annotated per-read via `d2_call`. `yields.tsv` records per-bucket read counts |

## Phase 7 — single-cell (dCODE)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| 10x dCODE dextramer donors (1–4) | HF dataset `isalgo/airr_benchmark`, path `dcode/vdj_v1_hs_aggregated_donor{1..4}_{all_contig_annotations,consensus_annotations,binarized_matrix}.csv.gz` | gzipped CSV (10x CellRanger VDJ all-contig / consensus annotations; CITE-seq binarized dextramer matrix: `barcode` + 50 `*_binder` boolean cols) | `huggingface_hub.hf_hub_download(repo_id="isalgo/airr_benchmark", repo_type="dataset", filename=...)` (fetched at test time into the HF cache; skips cleanly offline) | **experimental** — 10x Genomics "A New Way of Exploring Immunity" dCODE dextramer single-cell TCR + surface-marker dataset; antigen labels are the single-True `*_binder` column per cell |

## Phase 6 — biomarker association (Emerson benchmark + validation)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| Emerson HIP cohort (786 subjects) | HF dataset [`isalgo/airr_hip`](https://huggingface.co/datasets/isalgo/airr_hip) — redistributed from Adaptive immuneACCESS **Emerson-2017-NatGen** | per-subject VDJtools tables `corr/HIP#####.txt.gz` (`count freq cdr3nt cdr3aa v d j VEnd…`); `metadata.txt` (TAB-sep: `file_name sample_id age race sex cmv hla`) | `examples/emerson_cmv_hla.py` → `huggingface_hub.snapshot_download(repo_type="dataset", allow_patterns=["corr/{sample}.txt.gz"])`; ingest with `io.ingest_cohort(fmt="vdjtools")` | **experimental** TCRβ repertoires + phenotypes (Emerson et al., *Nat Genet* 2017, doi:10.1038/ng.3822). `cmv` ∈ {`+`,`-`,`NA`}; `hla` = 2-digit HLA-A/B only (`HLA-A*02`); ⚠ `race` contains commas — split on TAB. No discovery/validation split column |
| VDJdb (CMV validation target) | local checkout `/Users/mikesh/vcs/code/vdjdb-db/database/vdjdb.slim.txt` (canonical `antigenomics/vdjdb-db`, 2024-06 release; 2-digit HLA matches airr_hip) | TSV, 16 cols: `gene cdr3 species antigen.epitope antigen.gene antigen.species … v.segm j.segm … mhc.a mhc.b mhc.class … vdjdb.score` | read with polars; filter `gene==TRB & species==HomoSapiens & antigen.species~CMV` | **curated** TCR↔epitope database; newer 4-digit dump at `/Users/mikesh/vcs/code/vdjdb-iedb-concordance/vdjdb_dump_2026/vdjdb.slim.txt`. Canonical fetch: `antigenomics/vdjdb-db` GitHub releases |
| FMBA covid19 (TRA+TRB) — Phase 6b benchmark | aldan3 `/projects/fmba_covid` (clonotype tables `COV_V_usage_adjustment_v3/FMBA_functional/*.clonotypes.TRB.txt` + `data/*.clonotypes.{TRA,TRB}.pool.aa.table.txt`) ↔ private HF `isalgo/airr_covid19`; phenotype = `metadata_fmba_full.txt`/`desc_fmba_not_nan_hla.csv` (`COVID_status`, 4-digit `HLA-*`), join by 12-digit `id` | legacy VDJtools tables; TAB metadata | `scripts/biomarker_bench.sbatch -- covid19` on aldan3 (reads `/projects/fmba_covid` directly) | **experimental** deep TCRα/β covid repertoires (Vlasova et al., *Genome Medicine* 2026;18:20). COVID-association + α-β co-occurrence benchmark |
| covid19 biomarker oracle | HF `isalgo/airr_covid19` → `covid_associated_clonotypes.csv` (`cdr3,cluster,has_covid_association,chain,v,j`) | CSV | validation target for the covid19 association hits | **curated** — the study's published COVID-associated clonotype/cluster list (Vlasova 2026) |
| FMBA covid19_vacc (TRA+TRB) — Phase 6b benchmark | aldan3 `/projects/fmba_covid/vaccine/corr_func/*.clonotypes.{TRA,TRB}.txt` ↔ private HF `isalgo/airr_covid19_vacc`; phenotype = `vaccine/processed_metadata.tsv` / `vaccine/metadata.csv` (`timepoint` ∈ before/20d-after, `vaccine` ∈ GamCOVIDVac/CoviVac), join by 12-digit `id` | legacy VDJtools tables; CSV/TSV metadata | `scripts/biomarker_bench.sbatch -- covid19_vacc` on aldan3 | **experimental** pre/post-vaccination TCRα/β repertoires (Vlasova 2026). Timepoint/vaccine association benchmark |
## Phase 5 — preprocess (VJ batch-correction validation)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| FMBA covid Cohort I (TCRβ) | aldan3 HPC project `fmba_covid`: clonotype tables `COV_V_usage_adjustment_v3/FMBA_functional/*.clonotypes.TRB.txt`; raw usage `COV_V_usage_adjustment_v3/V_usage_FMBA.tsv`; paper-adjusted `…_adjusted.tsv`; batch labels = `unprocessed_fmba_metadata.csv` col `folder` (9 NovaSeq runs = the paper's 9 batches), join by 12-digit id in col `name` | legacy VDJtools tables (`count freq cdr3nt cdr3aa v d j VEnd…`); usage matrices = sample × V-gene/family TSV | `aldan3 pull /projects/fmba_covid/COV_V_usage_adjustment_v3/…` (needs the `aldan3` client + `fmba_covid` group). Validation driver: `appendix/validate_batch_covid.py` | **experimental** — deep TCRβ repertoires (median ~3.3M reads/sample) from Vlasova et al., *Genome Medicine* 2026;18:20 (10.1186/s13073-025-01589-4). Used to validate `preprocess.correct_vj_usage(transform="sigmoid")` + `apply_vj_correction`: batch η² 0.109→0.002, grand mean preserved, reads preserved |

## Phase 10 — BCR SHM & lineage (test data, planned; not yet fetched)

Verified via PubMed (2026-07); all four are owner (Shugay)-co-authored BCR datasets. Data paths are
**TBD** — fill the `How to obtain` details when each dataset is actually sourced; do not fabricate a
path before then. All are **experimental** repertoire data.

| Dataset | Reference (verified) | DOI | Relevance |
|---|---|---|---|
| Longitudinal memory-B / ASC BCR repertoires | Mikelov et al., *eLife* 2022;11:e79254 | [10.7554/eLife.79254](https://doi.org/10.7554/eLife.79254) | UMI-tagged BCR-seq + clonal-lineage/phylogeny — primary Phase-9 (MiGEC) + Phase-10 (lineage) target |
| TCGA tumor repertoires from RNA-seq | Bolotin et al., *Nat. Biotechnol.* 2017;35(10):908–911 | [10.1038/nbt.3979](https://doi.org/10.1038/nbt.3979) | RNA-seq-derived TCR/BCR (MiXCR method); TCGA IGH |
| CD27-dull/bright memory-B VH repertoires | Grimsholm et al., *Cell Rep.* 2020;30(9):2963–2977.e6 | [10.1016/j.celrep.2020.02.022](https://doi.org/10.1016/j.celrep.2020.02.022) | memory-B subsets, VH usage + SHM frequency (owner wrote "CD20-dull"; the paper is **CD27**-dull) |
| CVID peripheral B-cell selection | Grimsholm et al., *Cell Rep.* 2023;42(5):112446 | [10.1016/j.celrep.2023.112446](https://doi.org/10.1016/j.celrep.2023.112446) | CVID Ig-seq, peripheral B-cell selection |

**Flag (unresolved):** the owner's "Mikelov *allergy* paper" literally matches Mikelov et al.,
*Nat. Immunol.* 2025;26(12):2328–2342 ([10.1038/s41590-025-02323-3](https://doi.org/10.1038/s41590-025-02323-3),
peanut oral-immunotherapy) — but that study is **TCR / single-cell** and does **not** list Shugay
as an author, so the Shugay-co-authored **eLife 2022 BCR** paper is recorded above instead. Confirm
which Mikelov dataset is intended before use. (References verified via PubMed; DOIs copied verbatim.)

## Bundled precomputed models (shipped in the wheel)

Ship under `python/vdjtools/model/_bundled/<source>/<LOCUS>/` (parquet marginals + `manifest.json`);
loaded with `vdjtools.model.load_bundled(locus, source)`. ~0.4 MB total (30–150 KB/model).

| Model set | Origin | How to rebuild | Provenance |
|---|---|---|---|
| `olga` (7 loci) | `from_olga` on the OLGA default models above | `python appendix/build_bundled_models.py` (the OLGA part is inline in the precompute) | **derived** — OLGA generative parameters converted to the polars schema (keeps OLGA germline; exact-Pgen bootstrap, single-D) |
| `learned` (7 loci) | native EM (`infer_native`) on real out-of-frame HF reads (2 000 unique clonotypes/locus, 12 fixed iters, arda V/J[/D]-masked) seeded from the `olga` model | `python appendix/build_bundled_models.py` | **computed from experimental data** — real-repertoire gene-usage/trim/insertion marginals. D-bearing loci (IGH/TRD/TRB) carry an **arda-anchored tandem-D** event: a read may be `n_D=2` only where arda called a second D (`d2_call`), which counters the tandem-vs-long-insertion identifiability that inflates unregularized D-D EM (TRB 0.28→**0.00**, TRD 0.18→**0.006**, IGH **0.009**). `EM_SINGLE_D=1` / `ND_PRIOR` are alternative regularizers. |

## Golden fixtures (tests)

| Dataset | Origin | Use |
|---|---|---|
| IGoR model files (human/mouse, all loci) | `IGoR-models/` (local) | model-loader + Pgen oracle fixtures; canonical `model_parms`/`model_marginals`/anchors format |
| OLGA `default_models/` (7 human loci) | **vendored in-repo** at `tests/python/fixtures/olga/default_models/` — copied from `antigenomics/mirpy` (`main`/`legacy-v2`, `mir/resources/olga/default_models/`; byte-identical to OLGA's published models, Sethna et al.). Re-fetch: `git -C ../mirpy archive main mir/resources/olga/default_models \| tar -x`. Mouse loci not vendored (unused). | model-loader + native-vs-OLGA Pgen oracle (`from_olga`). Tests resolve them here (override with `VDJTOOLS_OLGA_MODELS`); the `olga` pip package (`[oracle]` / test extra) supplies the *runtime* Pgen comparison. Makes the suite self-contained — no external mirpy checkout needed. |
| Legacy input-format samples | `tests/python/fixtures/legacy/*.txt.gz` — copied from legacy vdjtools `src/test/resources/samples/` (on the `legacy-1.x` branch; re-fetch with `git show legacy-1.x:src/test/resources/samples/<f>.txt.gz`) | format-conversion conformance for `vdjtools.io.convert` (MiXcr v1/2+v3, MiGec, ImmunoSeq v1/v2, ImgtHighVQuest, Vidjil, RTCR). Oracles (input row → canonical values) derived from the legacy Groovy parsers; see `tests/python/test_convert.py` |
| TRUST4 report sample | `tests/python/fixtures/legacy/trust4.txt.gz` — **synthetic** (not from legacy vdjtools), hand-built from the literal TRUST4 `*_report.tsv` header (`#count frequency CDR3nt CDR3aa V D J C cid cid_full_length`) in [`liulab-dfci/TRUST4`](https://github.com/liulab-dfci/TRUST4) `trust-simplerep.pl`; CDR3 nt/aa reuse verified junctions from the MiXcr/RTCR fixtures | conformance for `read_trust4` (`tests/python/test_convert.py`), incl. `partial`/gene-less row skipping. arda AIRR output is exercised inline (arda writes standard AIRR; `read_arda` delegates to `read_airr`) |

## Carried-over resource data (from legacy vdjtools, `legacy-1.x` branch)

| File | Use |
|---|---|
| `resources/profile/aa_property_table.txt` | amino-acid physicochemical properties (Phase 3) |
| `resources/profile/cdr3contact.txt` | CDR3 contact-probability estimate (Phase 3) |
