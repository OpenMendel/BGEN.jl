struct Bgen
    io::IOStream
    fsize::UInt64
    #offset::UInt64

    header::Header
    samples::Samples
    variants::Union{Vector{Variant}, Nothing}
    # note: detailed information of variables stored in Variable struct
end

"""
    Bgen(path; sample_path=nothing, delay_parsing=false)

Read in the Bgen file information: header, list of samples.
Variants and genotypes are read separately.

- `parse_variants`: option to load variant information in memory
when a Bgen file is opened.
"""
function Bgen(path::AbstractString;
    sample_path::Union{AbstractString,Nothing}=nothing,
    parse_variants=true
    )
    io = open(path)
    fsize = filesize(path)
    # read header
    header = Header(io)

    # read samples
    if header.has_sample_ids
        samples = Samples(io, header.n_samples)
    elseif sample_path !== nothing
        samples = Samples(sample_path, header.n_samples)
    else
        samples = Samples(header.n_samples)
    end

    offset = header.offset + 4 # location of the first variant_ids

    if parse_variants
        variants = parse_all_variants(io, fsize, header)
    else
        variants = nothing
    end
    Bgen(io, fsize, header, samples, variants)
end

"""
    offset_first_variant(x)
returns the start of the first variant
"""
@inline function offset_first_variant(x::Bgen)
    return x.header.offset + 4
end

"""
    get_variant(io, fsize, h, offset)
    get_variant(x::Bgen, offset::Integer)
Parse a variant starting at position `offset`.
"""
function get_variant(io::IOStream, fsize::Integer, h::Header, offset::Integer)
    if eof(io) || offset >= fsize
        @error "reached end of file"
    end
    Variant(io, offset, h.compression, h.layout, h.n_samples)
end

function get_variant(x::Bgen, offset::Integer)
    get_variant(x.io, x.fsize, x.header, offset)
end

"""
    parse_all_variants(io, fsize, h)
    parse_all_variants(x::Bgen)
Parse all variants (but not genotypes) of the file
"""
function parse_all_variants(io::IOStream, fsize::Integer, h::Header)
    offset = offset_first_variant(h)
    variants = Vector{Variant}(undef, h.n_variants)
    for i in 1:h.n_variants
        variants[i] = get_variant(io, fsize, h, offset)
        offset = variants[i].next_var_offset
    end
    variants
end

function parse_all_variants(x::Bgen)
    parse_all_variants(x.io, x.fsize, x.header)
end
