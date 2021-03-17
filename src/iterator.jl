abstract type VariantIterator end

@inline function Base.eltype(vi::VariantIterator)
    Variant
end

"""
    VariantIteratorFromStart(b::Bgen)
Variant iterator that iterates from the beginning of the Bgen file
"""
struct VariantIteratorFromStart <: VariantIterator
    b::Bgen
end

function Base.iterate(vi::VariantIteratorFromStart,
                        state=offset_first_variant(vi.b))
    if state >= vi.b.fsize
        return nothing
    else
        v = Variant(vi.b, state)
        nextstate = v.next_var_offset
        return (v, nextstate)
    end
end

@inline function Base.length(vi::VariantIteratorFromStart)
    vi.b.header.n_variants
end

@inline function Base.size(vi::VariantIteratorFromStart)
    (vi.b.header.n_variants, )
end

"""
    VariantIteratorFromOffsets(b::Bgen, offsets::Vector{UInt})
Variant iterator that iterates over a vector of offsets
"""
struct VariantIteratorFromOffsets <: VariantIterator
    b::Bgen
    offsets::Vector{UInt}
end

function Base.iterate(vi::VariantIteratorFromOffsets, state=1)
    state > length(vi.offsets) ? nothing :
        (Variant(vi.b, vi.offsets[state]), state + 1)
end

@inline function Base.length(vi::VariantIteratorFromOffsets)
    length(vi.offsets)
end

@inline function Base.size(vi::VariantIteratorFromOffsets)
    size(vi, offsets)
end

struct Filter{F,I}
    itr::I
    min_maf::AbstractFloat
    min_hwe_pval::AbstractFloat
    min_success_rate_per_variant::AbstractFloat
    rmask::BitVector
    cmask::Union{Nothing,BitVector}
    decompressed::Union{Nothing, Vector{UInt8}}
end

"""
    BGEN.filter(itr; min_maf=NaN, min_hwe_pval=NaN, min_success_rate_per_variant=NaN, 
        cmask=trues(n_variants(itr.b)), rmask=trues(n_variants(itr.b)))
"Filtered" iterator for variants based on min_maf, min_hwe_pval, min_success_rate_per_variant, 
cmask, and rmask.
"""
filter(itr::VariantIterator; 
    min_maf=NaN, min_hwe_pval=NaN, min_success_rate_per_variant=NaN, 
    cmask = trues(n_variants(itr.b)), 
    rmask = nothing,
    decompressed = nothing) = 
    Filter(itr, min_maf, min_hwe_pval, min_success_rate_per_variant, 
        cmask, rmask, decompressed)

function iterate(f::Filter, state...)
    y = iterate(f.itr, state...)
    cnt = 1
    while y !== nothing
        v = y[1]
        if length(v.genotypes) == 0 || length(v.genotypes[1].decompressed) == 0
            decompressed = decompress(io, v, h; decompressed=f.decompressed)
        else
            decompressed = v.genotypes[1].decompressed
        end
        idx = [1]
        if length(v.genotypes) == 0
            p = parse_preamble!(decompressed, idx, h, v)
            push!(v.genotypes, Genotypes{T}(p, decompressed))
        else
            p = v.genotypes[1].preamble
            if h.layout == 2
                idx[1] += 10 + h.n_samples
            end
        end
        passed = true
        if !cmask[cnt]
            passed = false
        end
        if passed && !isnan(f.min_hwe_pval) 
            hwe_pval = hwe(p, d, idx, p.layout, rmask)
            if hwe_pval < f.min_hwe_pval 
                passed = false
            end
        end
        if passed && !isnan(f.min_maf)
            maf = maf(p, d, idx, p.layout, rmask)
            if maf < f.min_maf && passed
                passed = false
            end
        end
        if passed && !isnan(min_success_rate_per_variant)
            successes = length(intersect(p.missings, (1:n_samples(f.itr.b)[rmask])))
            success_rate = successes / count(rmask)
            if success_rate < f.min_success_rate_per_variant && passed 
                passed = false
            end
        end       
        if passed
            return y
        end
        y = iterate(f.itr, y[2])
        cnt += 1
    end
    nothing
end

eltype(::Type{Filter{F,I}}) where {F,I} = Variant
IteratorEltype(::Type{Filter{F,I}}) where {F,I} = IteratorEltype(I)
IteratorSize(::Type{<:Filter}) = SizeUnknown()

reverse(f::Filter) = Filter(reverse(f.itr), min_maf, min_hwe_pval, 
    min_success_rate_per_varinat, rmask, cmask, decompressed)
