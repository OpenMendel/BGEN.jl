# BGEN.jl
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://OpenMendel.github.io/BGEN.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://OpenMendel.github.io/BGEN.jl/stable)
[![codecov](https://codecov.io/gh/OpenMendel/BGEN.jl/branch/main/graph/badge.svg?token=W28QPREGC7)](https://codecov.io/gh/OpenMendel/BGEN.jl)
[![build Actions Status](https://github.com/OpenMendel/BGEN.jl/workflows/CI/badge.svg)](https://github.com/OpenMendel/BGEN.jl/actions)

Routines for reading compressed storage of genotyped or imputed markers

[*Genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) data with imputed markers are often saved in the [**BGEN format**](https://www.well.ox.ac.uk/~gav/bgen_format/) or `.bgen` file.
It can store both hard calls and imputed data, unphased genotypes and phased haplotypes. Each variant is compressed separately to make indexing simple. An index file (`.bgen.bgi`) may be provided to access each variant easily. [UK Biobank](https://www.ukbiobank.ac.uk/) uses this format for genome-wide imputed genotypes.

## Installation

This package requires Julia v1.0 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.


This package is registered in the default Julia package registry, and can be installed through standard package installation procedure: e.g., running the following code in Julia REPL.
```julia
using Pkg
pkg"add BGEN"
```

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. OPENMENDEL: a cooperative programming project for statistical genetics. Hum Genet. 2020 Jan;139(1):61-71. doi: 10.1007/s00439-019-02001-z. Epub 2019 Mar 26. PMID: 30915546; PMCID: [PMC6763373](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6763373/).*

## Acknowledgments

Current implementation incorporates ideas in a [bgen Python package](https://github.com/jeremymcrae/bgen).

This project has been supported by the National Institutes of Health under awards R01GM053275, R01HG006139, R25GM103774, and 1R25HG011845.
