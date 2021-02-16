struct Header
    offset::UInt32
    header_length::UInt32
    n_variants::UInt32
    n_samples::UInt32
    compression::UInt8
    layout::UInt8
    has_sample_ids::Bool
end

struct Samples
    samples :: Vector{String}
end

struct Index
    path::String
    db::SQLite.DB
    offsets::Vector{UInt32}
    rsids::Vector{String}
    chroms::Vector{String}
    positions::Vector{UInt32}
end

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

struct Variant
    offset::UInt64
    geno_offset::UInt64 # to the start of genotype block
    next_var_offset::UInt64
    geno_block_size::UInt32
    n_samples::UInt32
    varid::String
    rsid::String
    chrom::String
    pos::UInt32
    n_alleles::UInt16
    alleles::Vector{String}
    # length-1 for parsed one, empty array for not parsed yet or destroyed,
    genotypes::Vector{Genotypes}
end

struct Bgen
    io::IOStream
    fsize::UInt64

    header::Header
    samples::Samples
    idx::Union{Index, Nothing}
    # note: detailed information of variables stored in Variable struct,
    # accessed thru VariantIterator
end
