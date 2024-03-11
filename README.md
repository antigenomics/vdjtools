> VDJtools is free for academic use. his is a legacy version, the software is no longer supported, if you need similar software solutions please visit [MiLaboratory](https://milaboratory.ru/).

[![JitPack](https://jitpack.io/v/mikessh/vdjtools.svg)](https://jitpack.io/#mikessh/vdjtools)

## VDJtools

A comprehensive analysis framework for T-cell and B-cell repertoire sequencing data. Compiled binaries are available from [here](https://github.com/mikessh/vdjtools/releases/latest). You can download them and execute as

```bash
java -jar vdjtools.jar ...
```

Make sure that you've specified the full/correct path to jar file. In case of Java Heap Space exception, you can increase the JVM memory limit by adding ``-Xmx20G`` (for extra 20G) after the ``-jar`` argument.

The software is cross-platform and requires Java v1.8+ to run and R to perform plotting.

Easy installation on **MacOS/Linux** via [Homebrew](http://brew.sh/) or [Linuxbrew](http://linuxbrew.sh/):
```bash
brew tap homebrew/science
brew tap mikessh/repseq
brew install vdjtools
vdjtools CalcBasicStats ...
```
See [homebrew-repseq](https://github.com/mikessh/homebrew-repseq) for other RepSeq analysis software Homebrew installers.

List of features and detailed documentation can be found at [ReadTheDocs](http://vdjtools-doc.readthedocs.io).

Example datasets and shell scripts are provided in a separate [repository](https://github.com/mikessh/vdjtools-examples).

### Please cite VDJtools as:

Shugay M et al. VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires. [PLoS Comp Biol 2015; 11(11):e1004503-e1004503](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004503).

### Some recent publications where VDJtools was used:

- Feng Y et al. A mechanism for expansion of regulatory T-cell repertoire and its role in self-tolerance. [Nature 2015; doi:10.1038/nature16141](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature16141.html)
- Britanova OV et al. Dynamics of Individual T Cell Repertoires: From Cord Blood to Centenarians. [J Immunol 2016; doi:10.4049/jimmunol.1600005](http://www.jimmunol.org/content/196/12/5005.short)
- Joachims ML et al. Single-cell analysis of glandular T cell receptors in Sjögren’s syndrome. [JCI Insight 2016; doi:10.1172/jci.insight.85609](https://insight.jci.org/articles/view/85609)
- Plitas G et al. Regulatory T cells exhibit distinct features in human breast cancer. [Immunity 2017; doi.org:10.1016/j.immuni.2016.10.032](http://www.sciencedirect.com/science/article/pii/S1074761316304435)
- Izraelson M et al. Comparative Analysis of Murine T Cell Receptor Repertoires. [Immunology 2017; doi.org:10.1111/imm.12857](http://onlinelibrary.wiley.com/doi/10.1111/imm.12857/full)
- Bolotin DA et al. Antigen receptor repertoire profiling from RNA-seq data. [Nat Biotech 2017; doi:10.1038/nbt.3979](https://www.nature.com/articles/nbt.3979)
- Meng W et al. An atlas of B-cell clonal distribution in the human body. [Nat Biotech 2017; doi:10.1038/nbt.3942](https://www.nature.com/articles/nbt.3942)