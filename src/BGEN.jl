module BGEN
import Base: length, getindex, setindex, firstindex, lastindex
include("header.jl")
include("sample.jl")
include("variant.jl")
include("bgen.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
