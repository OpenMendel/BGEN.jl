module Bgen

include("header.jl")
include("bgen.jl")

datadir(parts...) = joinpath(@__DIR__, "..", "data", parts...)

end
