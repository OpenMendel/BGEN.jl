module BGEN
import Base: length, getindex, setindex, firstindex, lastindex, eltype, size,
            iterate, close, Iterators.filter
import Tables: columntable
import Statistics: mean
import SpecialFunctions: gamma_inc
import TranscodingStreams: initialize, finalize, buffermem, process, Buffer, Error
export Bgen, Samples, Variant, Genotypes, Index
export io, fsize, samples, n_samples, n_variants, compression
export varid, rsid, chrom, pos, n_alleles, alleles, minor_allele, major_allele
export phased, min_ploidy, max_ploidy, ploidy, bit_depth, missings
export parse_variants, iterator, probabilities!, minor_allele_dosage!
export check_decompressed_length
export ref_allele_dosage!, clear!
export select_region, variant_by_rsid, variant_by_pos, variant_by_index
export rsids, chroms, positions
export hwe, maf, info_score, counts!
export VariantIteratorFromStart, VariantIteratorFromOffsets
using CodecZlib, CodecZstd, SQLite, SIMD
include("structs.jl")
include("iterator.jl")
include("header.jl")
include("minor_certain.jl")
include("sample.jl")
include("variant.jl")
include("bgen_ftns.jl")
include("genotypes.jl")
include("index.jl")
include("utils.jl")
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
