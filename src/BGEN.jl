module BGEN

include("bgen.jl")
include("header.jl")
include("sample.jl")
include("variant.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
