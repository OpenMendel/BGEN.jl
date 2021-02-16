"""
    Bgen(path; sample_path=nothing, delay_parsing=false)

Read in the Bgen file information: header, list of samples.
Variants and genotypes are read separately.

- `path`: path to the ".bgen" file.
- `sample_path`: path to  ".sample" file, if applicable.
- `idx_path`: path to ".bgi" file, defaults to `path * ".bgi`.
"""
function Bgen(path::AbstractString;
    sample_path = nothing,
    idx_path = isfile(path * ".bgi") ? path * ".bgi" : nothing
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
            @error "$idx_path is not a file"
        end
    else
        idx = nothing
    end

    Bgen(io, fsize, header, samples, idx)
end

"""
    iterator(b::Bgen; offsets=nothing)
Retrieve a variant iterator for `b`. If offsets === nothing, it returns:

- a `VariantIteratorFromStart` if `.bgi` file was not provided, or
`from_bgen_start` is `true`.
- `VariantIteratorFromOffsets` containing the offsets of each variant if `.bgi`
was provided.
- If `offsets` is provided, it returns a `VariantIteratorFromOffsets`.

"""
function iterator(b::Bgen; offsets=nothing, from_bgen_start=false)
    if offsets === nothing
        if b.idx === nothing || from_bgen_start
            return VariantIteratorFromStart(b)
        else
            return VariantIteratorFromOffsets(b, BGEN.offsets(b.idx))
        end
    else
        return VariantIteratorFromOffsets(b, offsets)
    end
end

"""
    offset_first_variant(x)
returns the offset of the first variant
"""
@inline function offset_first_variant(x::Bgen)
    return x.header.offset + 4
end

"""
    parse_variants(b::Bgen; offsets=offsets)
Parse variants of the file.
"""
function parse_variants(b::Bgen; offsets=nothing, from_bgen_start=false)
    collect(iterator(b; offsets=offsets, from_bgen_start=from_bgen_start))
end

"""
    rsids(b; vi=nothing, offsets=nothing)
Get rsid list of all variants. If `vi === nothing`, it will return based on
iterator(b; offsets=offsets).

Arguments:

- `bgen`: `Bgen` object
- `vi`: a collection of `Variant`s
- `offsets`: offset of each variant to be returned
"""
function rsids(b::Bgen; vi=nothing, offsets=nothing, from_bgen_start=false)
    if vi === nothing
        if b.idx !== nothing && offsets === nothing
            return rsids(b.idx)
        else
            vi = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        end
    else
        collect(v.rsid for v in vi)
    end
end

"""
    chroms(bgen; vi=nothing)
Get chromosome list of all variants. If `vi === nothing`, it will return based on
iterator(b; offsets=offsets).

Arguments:

- `bgen`: `Bgen` object
- `vi`: a collection of `Variant`s
- `offsets`: offset of each variant to be returned
"""
function chroms(b::Bgen; vi=nothing, offsets=nothing, from_bgen_start=false)
    if vi === nothing
        if b.idx !== nothing && offsets === nothing
            return chroms(b.idx)
        else
            vi = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        end
    else
        collect(v.chrom for v in vi)
    end
end

"""
    positions(bgen; vi=nothing)
Get base pair positions of all variants. If `vi === nothing`, it will return based on
iterator(b; offsets=offsets).

Arguments:

- `bgen`: `Bgen` object
- `vi`: a collection of `Variant`s
- `offsets`: offset of each variant to be returned
"""
function positions(b::Bgen; vi=nothing, offsets=nothing,
    from_bgen_start=false)::Vector{UInt32}

    if vi === nothing
        if b.idx !== nothing && offsets === nothing
            return positions(b.idx)
        else
            vi = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        end
    else
        collect(v.pos for v in vi)
    end
end
