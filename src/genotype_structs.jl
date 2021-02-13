struct Preamble
    n_samples::UInt32
    n_alleles::UInt16
    phased::UInt8
    min_ploidy::UInt8
    max_ploidy::UInt8
    ploidy::Vector{UInt8}
    bit_depth::UInt8
    max_probs::Int
    missings::Vector{Int}
end

struct Genotypes{T}
    preamble::Preamble # Once parsed, it will not be destroyed unless Genotypes is destroyed
    decompressed::Vector{UInt8}
    probs::Vector{T}
    minor_idx::Vector{UInt8} # index of minor allele
    dose::Vector{T}
end

function Genotypes{T}(p::Preamble, d::Vector{UInt8}) where T <: AbstractFloat
    Genotypes{T}(p, d, T[], UInt8[0], T[])
end
