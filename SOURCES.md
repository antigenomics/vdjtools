# SOURCES

Provenance of every dataset used or produced by vdjtools. Never guess a source — record it here.

## Phase 1 — model engine (bootstrap)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| OLGA default models (7 human loci + mouse) | `mirpy/mir/resources/olga/default_models/{human_T_alpha,beta,gamma,delta, human_B_heavy,kappa,lambda}/` (reference only — do not copy mirpy code) | `model_params.txt`, `model_marginals.txt`, `V/J_gene_CDR3_anchors.csv` | read from the local mirpy checkout | published OLGA models (Sethna et al.); **derived** generative-model parameters |
| Synthetic out-of-frame training seqs (~100k/chain) | OLGA generation via mirpy tooling; sample at `mirpy/tests/assets/olga_humanTRB_1000.txt.gz` | TSV (V call, J call, CDR3 nt) | generate with OLGA from the models above | **computed** (Monte-Carlo draws); **no tandem D by design** — see plan Phase 1 note |
| Real non-functional seqs (all 7 loci) | **owner-provided, ships later** | TBD | TBD | **experimental**; the data that will actually exercise D-D — update this row on arrival |

## Phase 7 — single-cell (dCODE)

| Dataset | Origin | Format | How to obtain | Provenance |
|---|---|---|---|---|
| 10x dCODE dextramer donors (1–4) | HF dataset `isalgo/airr_benchmark`, path `dcode/vdj_v1_hs_aggregated_donor{1..4}_{all_contig_annotations,consensus_annotations,binarized_matrix}.csv.gz` | gzipped CSV (10x CellRanger VDJ all-contig / consensus annotations; CITE-seq binarized dextramer matrix: `barcode` + 50 `*_binder` boolean cols) | `huggingface_hub.hf_hub_download(repo_id="isalgo/airr_benchmark", repo_type="dataset", filename=...)` (fetched at test time into the HF cache; skips cleanly offline) | **experimental** — 10x Genomics "A New Way of Exploring Immunity" dCODE dextramer single-cell TCR + surface-marker dataset; antigen labels are the single-True `*_binder` column per cell |

## Golden fixtures (tests)

| Dataset | Origin | Use |
|---|---|---|
| IGoR model files (human/mouse, all loci) | `IGoR-models/` (local) | model-loader + Pgen oracle fixtures; canonical `model_parms`/`model_marginals`/anchors format |
| OLGA `default_models/` | OLGA pip package | Pgen oracle (`pip install olga`), as seqtree CI does |
| Legacy input-format samples | legacy vdjtools `src/test/resources/samples/*.txt.gz` (on the `legacy-1.x` branch) | format-conversion conformance (MiXcr, MiGec, ImmunoSeq v1/v2, ImgtHighVQuest, Vidjil, RTCR, …) |

## Carried-over resource data (from legacy vdjtools, `legacy-1.x` branch)

| File | Use |
|---|---|
| `resources/profile/aa_property_table.txt` | amino-acid physicochemical properties (Phase 3) |
| `resources/profile/cdr3contact.txt` | CDR3 contact-probability estimate (Phase 3) |
