# BGEN.jl

Routines for reading compressed storage of genotyped or imputed markers

[*Genome-wide association studies (GWAS)*](https://en.wikipedia.org/wiki/Genome-wide_association_study) data with imputed markers are often saved in the [**BGEN format**](https://www.well.ox.ac.uk/~gav/bgen_format/) or `.bgen` file.
It can store both hard calls and imputed data, unphased genotypes and phased haplotypes. Each variant is compressed separately to make indexing simple. An index file (`.bgen.bgi`) may be provided to access each variant easily. [UK Biobank](https://www.ukbiobank.ac.uk/) uses this format for genome-wide imputed genotypes.

## Installation

This package requires Julia v1.0 or later, which can be obtained from
https://julialang.org/downloads/ or by building Julia from the sources in the
https://github.com/JuliaLang/julia repository.

The package has not yet been registered and must be installed using the repository location. 
It can be done with the following Julia code:
```julia
using Pkg
pkg"https://github.com/OpenMendel/BGEN.jl"
```

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*OPENMENDEL: a cooperative programming project for statistical genetics. Zhou H, Sinsheimer JS, Bates DM, Chu BB, German CA, Ji SS, Keys KL, Kim J, Ko S, Mosher GD, Papp JC, Sobel EM, Zhai J, Zhou JJ, Lange K. Hum Genet. 2019 Mar 26. doi: 10.1007/s00439-019-02001-z. [Epub ahead of print] PMID: 30915546*

## Acknowledgments

Current implementation incorporates ideas in a [bgen Python package](https://github.com/jeremymcrae/bgen).

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
