struct Header
    offset::UInt32
    header_length::UInt32
    n_variants::UInt32
    n_samples::UInt32
    compression::UInt8
    layout::UInt8
    has_sample_ids::Bool
end

const Samples = Vector{String}

struct Index
    path::String
    db::SQLite.DB
    offsets::Vector{UInt64}
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
    ploidy::Union{UInt8, Vector{UInt8}}
    bit_depth::UInt8
    max_probs::Int
    missings::Vector{Int}
end

mutable struct Genotypes{T}
    preamble::Preamble # Once parsed, it will not be destroyed unless Genotypes is destroyed
    decompressed::Union{Nothing,Vector{UInt8}}
    probs::Union{Nothing, Vector{T}}
    minor_idx::UInt8 # index of minor allele
    dose::Union{Nothing, Vector{T}}
    dose_mean_imputed::Bool
    minor_allele_dosage::Bool
end

mutable struct BgenVariant <: GeneticVariantBase.Variant
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
    genotypes::Union{Nothing, Genotypes}
end

struct Bgen <: GeneticVariantBase.GeneticData
    io::IOStream
    fsize::UInt64

    header::Header
    samples::Samples
    idx::Union{Index, Nothing}
    # note: detailed information of variables stored in Variable struct,
    # accessed thru VariantIterator
    ref_first::Bool
end
