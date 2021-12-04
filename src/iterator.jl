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
    size(vi.offsets)
end

struct Filter{I, T} <: VariantIterator
    itr::I
    min_maf::AbstractFloat
    min_hwe_pval::AbstractFloat
    min_info_score::AbstractFloat
    min_success_rate_per_variant::AbstractFloat
    rmask::Union{Nothing,BitVector}
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
    T=Float32,
    min_maf=NaN, min_hwe_pval=NaN, min_info_score=NaN,
    min_success_rate_per_variant=NaN, 
    cmask = trues(n_variants(itr.b)), 
    rmask = nothing,
    decompressed = nothing) = 
    Filter{typeof(itr), T}(itr, min_maf, min_hwe_pval, min_info_score, min_success_rate_per_variant, 
        rmask, cmask, decompressed)

function Base.iterate(f::Filter{I,T}, state...) where {I,T}
    io, h = f.itr.b.io, f.itr.b.header
    if state !== ()
        y = Base.iterate(f.itr, state[1][2:end]...)
        cnt = state[1][1]
    else
        y = Base.iterate(f.itr)
        cnt = 1
    end
    while y !== nothing
        v, s = y
        passed = true
        if !f.cmask[cnt]
            passed = false
        end
        if  passed && (v.genotypes === nothing || v.genotypes.decompressed === nothing)
            decompressed = decompress(io, v, h; decompressed=f.decompressed)
        elseif passed
            decompressed = v.genotypes.decompressed
        end
        startidx = 1
        if passed && v.genotypes === nothing
            p = parse_preamble(decompressed, h, v)
            v.genotypes = Genotypes{T}(p, decompressed)
        elseif passed
            p = v.genotypes.preamble
        end
        if passed && h.layout == 2
            startidx += 10 + h.n_samples
        end
        if passed && !isnan(f.min_maf)
            current_maf = maf(p, v.genotypes.decompressed, startidx, h.layout, f.rmask)
            if current_maf < f.min_maf && passed
                passed = false
            end
        end
        if passed && !isnan(f.min_hwe_pval) 
            hwe_pval = hwe(p, v.genotypes.decompressed, startidx, h.layout, f.rmask)
            if hwe_pval < f.min_hwe_pval 
                passed = false
            end
        end
        if passed && !isnan(f.min_info_score)
            current_info_score = info_score(p, v.genotypes.decompressed, startidx, h.layout, f.rmask)
            if current_info_score < f.min_info_score && passed
                passed = false
            end
        end
        if passed && !isnan(f.min_success_rate_per_variant)
            successes = length(intersect(p.missings, (1:n_samples(f.itr.b)[f.rmask])))
            success_rate = successes / count(f.rmask)
            if success_rate < f.min_success_rate_per_variant && passed 
                passed = false
            end
        end    
        cnt += 1
        if passed
            return v, (cnt, s...)
        end
        y = iterate(f.itr, s...)
    end
    nothing
end

eltype(::Type{Filter{I,T}}) where {I,T} = Variant
IteratorEltype(::Type{Filter{I,T}}) where {I,T} = IteratorEltype(I)
IteratorSize(::Type{<:Filter}) = SizeUnknown()

reverse(f::Filter) = Filter(reverse(f.itr), f.min_maf, f.min_hwe_pval, 
    f.min_info_score, f.min_success_rate_per_varinat, f.rmask, f.cmask, f.decompressed)
