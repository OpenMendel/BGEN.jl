module BGEN
import Base: length, getindex, setindex, firstindex, lastindex
using CodecZlib, CodecZstd
include("header.jl")
include("sample.jl")
include("variant.jl")
include("bgen.jl")
include("genotypes.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
