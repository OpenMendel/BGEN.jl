module BGEN
import Base: length, getindex, setindex, firstindex, lastindex
export Bgen, parse_probs!, parse_dosage!, clear!
using CodecZlib, CodecZstd
include("header.jl")
include("minor_certain.jl")
include("genotype_structs.jl")
include("sample.jl")
include("variant.jl")
include("bgen.jl")
include("genotypes.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
