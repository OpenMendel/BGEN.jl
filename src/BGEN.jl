module BGEN
import Base: length, getindex, setindex, firstindex, lastindex
import Tables: columntable
export Bgen, Samples, Variant, Genotypes, Index
export parse_probs!, parse_dosage!, clear!
export rsids, chroms, positions
using CodecZlib, CodecZstd, SQLite
include("structs.jl")
include("header.jl")
include("minor_certain.jl")
include("sample.jl")
include("variant.jl")
include("bgen.jl")
include("genotypes.jl")
include("index.jl")
datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
