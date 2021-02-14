"""
    Bgen(path; sample_path=nothing, delay_parsing=false)

Read in the Bgen file information: header, list of samples.
Variants and genotypes are read separately.

- `sample_path`: path to  ".sample" file, if applicable
- `idx_path`: path to ".bgi" file, if applicable
- `parse_variants`: option to load variant information in memory
when a Bgen file is opened.
"""
function Bgen(path::AbstractString;
    sample_path = nothing,
    idx_path = path * ".bgi",
    parse_variants = true
    )
    io = open(path)
    fsize = filesize(path)
    # read header
    header = Header(io)
    # read samples
    if sample_path !== nothing
        samples = Samples(sample_path, header.n_samples)
        # disregard sample names in the header if .sample file is provided
        if header.has_sample_ids
            header_sample_length = read(io, UInt32)
            read(io, header_sample_length)
        end
    elseif header.has_sample_ids
        samples = Samples(io, header.n_samples)
    else
        samples = Samples(header.n_samples)
    end

    offset = header.offset + 4 # location of the first variant_ids

    if idx_path !== nothing
        if isfile(idx_path)
            idx = Index(idx_path)
        else
            idx = nothing
        end
    else
        idx = nothing
    end

    if parse_variants
        variants = parse_all_variants(io, fsize, header)
    else
        variants = Variant[]
    end
    Bgen(io, fsize, header, samples, variants, idx)
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

"""
    rsids(bgen)
Get rsid list of all variants
"""
function rsids(b::Bgen)
    if b.idx !== nothing
        return rsids(b.idx)
    elseif length(b.variants) != 0
        return [v.rsid for v in b.variants]
    else
        @error "need either .bgi file or parsed all variants"
    end
end

"""
    chroms(bgen)
Get chromosome list of all variants
"""
function chroms(b::Bgen)
    if b.idx !== nothing
        return chroms(b.idx)
    elseif length(b.variants) != 0
        return [v.chrom for v in b.variants]
    else
        @error "need either .bgi file or parsed all variants"
    end
end

"""
    positions(idx)
Get base pair positions of all variants
"""
function positions(b::Bgen)
    if b.idx !== nothing
        return positions(b.idx)
    elseif length(b.variants) != 0
        return [v.pos for v in b.variants]
    else
        @error "need either .bgi file or parsed all variants"
    end
end
