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
    idx_path = isfile(path * ".bgi") ? path * ".bgi" : nothing,
    ref_first = true
    )
    io = open(path)
    fsize = filesize(path)
    # read header
    header = Header(io)
    # read samples
    if sample_path !== nothing
        samples = get_samples(sample_path, header.n_samples)
        # disregard sample names in the header if .sample file is provided
        if header.has_sample_ids
            header_sample_length = read(io, UInt32)
            read(io, header_sample_length)
        end
    elseif header.has_sample_ids
        samples = get_samples(io, header.n_samples)
    else
        samples = get_samples(header.n_samples)
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

    Bgen(io, fsize, header, samples, idx, ref_first)
end

@inline io(b::Bgen) = b.io
Base.close(b::Bgen) = close(b.io)
@inline fsize(b::Bgen)::Int = b.fsize
@inline samples(b::Bgen) = b.samples
@inline n_samples(b::Bgen)::Int = b.header.n_samples
@inline n_variants(b::Bgen)::Int = b.header.n_variants
const compression_modes = ["None", "Zlib", "Zstd"]
@inline function compression(b::Bgen)
    compression_modes[b.header.compression + 1]
end

"""
    iterator(b::Bgen; offsets=nothing, from_bgen_start=nothing)
Retrieve a variant iterator for `b`.

- If `offsets` is provided, or `.bgen.bgi` is provided and
`from_bgen_start` is `false`, it returns a `VariantIteratorFromOffsets`,
iterating over the list of offsets.
- Otherwise, it returns a `VariantIteratorFromStart`, iterating from the start
of bgen file to the end of it sequentially.
"""
function iterator(b::Bgen; offsets=nothing, from_bgen_start=false)
    if offsets === nothing
        if b.idx === nothing || from_bgen_start
            return BgenVariantIteratorFromStart(b)
        else
            return BgenVariantIteratorFromOffsets(b, BGEN.offsets(b.idx))
        end
    else
        return BgenVariantIteratorFromOffsets(b, offsets)
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
function parse_variants(v::BgenVariantIterator)
    collect(v)
end

"""
    rsids(vi)
    rsids(b; offsets=nothing, from_bgen_start=false)
Get rsid list of all variants.

Arguments:
- `vi`: a collection of `Variant`s
- `bgen`: `Bgen` object
- `offsets`: offset of each variant to be returned
"""
function rsids(b::Bgen; offsets=nothing, from_bgen_start=false)
    if b.idx !== nothing && offsets === nothing
        return rsids(b.idx)
    else
        vi = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        rsids(vi)
    end
end
function rsids(vi::BgenVariantIterator)
    collect(v.rsid for v in vi)
end

"""
    chroms(vi)
    chroms(bgen; offsets=nothing)
Get chromosome list of all variants.

Arguments:
- `vi`: a collection of `Variant`s
- `bgen`: `Bgen` object
- `offsets`: offset of each variant to be returned
"""
function chroms(b::Bgen; vi=nothing, offsets=nothing, from_bgen_start=false)
    if b.idx !== nothing && offsets === nothing
        return chroms(b.idx)
    else
        vi = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        chroms(vi)
    end
end
function chroms(vi::BgenVariantIterator)
    collect(v.chrom for v in vi)
end

function chrom(b::Bgen, v::BgenVariant)
    chrom(v)
end

"""
    positions(vi)
    positions(bgen; offsets=nothing)
Get base pair positions of all variants.

Arguments:
- `vi`: a collection of `Variant`s
- `bgen`: `Bgen` object
- `offsets`: offset of each variant to be returned
"""
function positions(b::Bgen; offsets=nothing,
    from_bgen_start=false)::Vector{Int}
    if b.idx !== nothing && offsets === nothing
        return positions(b.idx)
    else
        vi = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        positions(vi)
    end
end
function positions(vi::BgenVariantIterator)
    collect(v.pos for v in vi)
end

function pos(b::Bgen, v::BgenVariant)
    pos(v)
end

function rsid(b::Bgen, v::BgenVariant)
    rsid(v)
end

function alleles(b::Bgen, v::BgenVariant)
    alleles(v)
end

function alt_allele(b::Bgen, v::BgenVariant)
    allele_list = alleles(v)
    b.ref_first ? allele_list[2] : allele_list[1]
end

function ref_allele(b::Bgen, v::BgenVariant)
    allele_list = alleles(v)
    b.ref_first ? allele_list[1] : allele_list[2]
end
