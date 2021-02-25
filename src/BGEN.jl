module BGEN
import Base: length, getindex, setindex, firstindex, lastindex, eltype, size,
            iterate, close
import Tables: columntable
import Statistics: mean
export Bgen, Samples, Variant, Genotypes, Index
export io, fsize, samples, n_samples, n_variants, compression
export varid, rsid, chrom, pos, n_alleles, alleles, minor_allele, major_allele
export phased, min_ploidy, max_ploidy, ploidy, bit_depth, missings
export parse_variants, iterator, probabilities!, minor_allele_dosage!, clear!
export select_region, variant_by_rsid, variant_by_pos, variant_by_index
export rsids, chroms, positions
export VariantIteratorFromStart, VariantIteratorFromOffsets
using CodecZlib, CodecZstd, SQLite
include("structs.jl")
include("iterator.jl")
include("header.jl")
include("minor_certain.jl")
include("sample.jl")
include("variant.jl")
include("bgen_ftns.jl")
include("genotypes.jl")
include("index.jl")
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
