module BGEN
import Base: length, getindex, setindex, firstindex, lastindex, eltype, size,
            iterate
import Tables: columntable
export Bgen, Samples, Variant, Genotypes, Index
export parse_variants, iterator, probabilities!, minor_allele_dosage!, clear!
export select_region, variant_by_rsid, variant_by_pos, rsids, chroms, positions
export VariantIteratorFromStart, VariantIteratorFromOffsets
using CodecZlib, CodecZstd, SQLite
include("structs.jl")
include("iterator.jl")
include("header.jl")
include("minor_certain.jl")
include("sample.jl")
include("variant.jl")
include("bgen.jl")
include("genotypes.jl")
include("index.jl")
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
