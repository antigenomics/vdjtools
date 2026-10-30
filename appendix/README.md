# Model-engine appendix

[Read the compiled appendix](murugan_model.pdf) · [LaTeX source](murugan_model.tex)

The appendix describes the model engine in vdjtools 5.4.0: exact nucleotide and
amino-acid Pgen, joint single/tandem-D histories, sequence generation,
error-aware inference, pruning and model diagnostics. It uses AIRR **junctions**
including the conserved anchors; IMGT CDR3 excludes them.
The current API guide is [Model workshop](https://docs.isalgo.dev/vdjtools/model.html).

## Build

```sh
make -C appendix
```

Requires LuaLaTeX, latexmk, BibTeX and the OldStandard text/math fonts.
The two network PDFs are frozen historical illustrations and declared build inputs.
`make` only typesets the document; it does not fit models or run benchmarks.

## Library builders and validators

Run scripts from the repository root in its configured environment.
Full fitting and concordance workloads belong on scheduled compute nodes.

| File | Purpose |
|---|---|
| `build_bundled_models.py` | Build twelve isalgo contexts: seven human loci and mouse TRA/TRB/IGH/IGK/IGL. Reannotates raw reads; native EM fits a retained-germline substitution rate and the supported history distribution. Defaults and environment options are documented in the script. |
| `compare_models.py` | Separate the identical-raw-OLGA engine oracle from independent functional-source software validation. Uses the accepted first-20000-reads/first-300-eligible-junction protocol. Functional-source likelihood is selection-shaped, not unbiased held-out recombination likelihood. |
| `concordance.py` | Historical seven-human-locus OLGA engine comparison and bounded synonymous-codon checks. Does not train models or test tandem inference. |
| `build_motif_vocab.py` | Build frozen motif/V dictionaries for statistical feature bases. |
| `build_adaptive_imgt_map.py` | Build the Adaptive-to-IMGT mapping using reference annotations. |
| `bench_junction_pipeline.py` | Compare reconstructed junctions to observed nucleotides and preliminary annotation evidence; D assignments are not independent biological truth. |

The performance producer `bench/bench_pgen.py` and historical inference diagnostic
`bench/bench_em.py` live in the separate `2026-vdjtools-benchmark` project.
They are consumers, not part of this library directory.
The appendix preserves their recorded timing and entropy values with their
historical parameter scope. It does not present the old training-report
likelihood as a held-out result or the old aligned-call fit as a current isalgo fit.

## References checked during reconciliation

The bibliography entries were checked against PubMed records
[22988065](https://pubmed.ncbi.nlm.nih.gov/22988065/),
[30657870](https://pubmed.ncbi.nlm.nih.gov/30657870/),
[29422654](https://pubmed.ncbi.nlm.nih.gov/29422654/) and
[30323809](https://pubmed.ncbi.nlm.nih.gov/30323809/).
The latent-history account also follows the Murugan 2012 Supporting Information,
particularly its scenario enumeration and likelihood updates.
IGoR's [`GenModel.cpp`](https://github.com/qmarcou/IGoR/blob/master/igor_src/GenModel.cpp)
documents its traversal-relative probability cutoff; vdjtools uses its own
factorised, global-best-cell rule and records actual omitted mass.
