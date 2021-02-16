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
