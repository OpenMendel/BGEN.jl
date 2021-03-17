"""
    hwe(b::Bgen, v::Variant; T=Float32, decompressed=nothing)
Hardy-Weinberg equilibrium test
"""
function hwe(b::Bgen, v::Variant; T=Float32, decompressed=nothing, rmask=nothing)
    io, h = b.io, b.header
    if length(v.genotypes) == 0 || length(v.genotypes[1].decompressed) == 0
        decompressed = decompress(io, v, h; decompressed=decompressed)
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

    hwe(p, decompressed, idx, h.layout, rmask)
end

function hwe(p::Preamble, d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8, 
    rmask::Union{Nothing, Vector{UInt16}})
    @assert layout == 2 "hwe only supported for layout 2"
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy 
    @assert length(p.missings) == 0 "current implementation does not allow missingness"
    idx1 = idx[1]

    # "counts" times 255.
    n00 = 0
    n01 = 0
    n11 = 0
    if p.n_samples >= 16
        @inbounds for n in 1:16:(p.n_samples - p.n_samples % 16)
            idx_base = idx1 + ((n-1) >> 1) << 2
            if rmask !== nothing
                rs = vload(Vec{16,UInt16}, rmask, n)
                if sum(rs) == 0
                    continue
                end
            end
            r = reinterpret(Vec{16, UInt16}, vload(Vec{32, UInt8}, d, idx_base))
            first  = (r & mask_odd)
            second = (r & mask_even) >> 8
            third  = 0x00ff - first - second
            if rmask !== nothing
                first = first * rs
                second = second * rs
                third = third * rs
            end
            n00 += sum(first)
            n01 += sum(second)
            n11 += sum(third)
        end
    end
    rem = p.n_samples % 16
    if rem != 0
        @inbounds for n in ((p.n_samples - rem) + 1) : p.n_samples
            rmask !== nothing && rmask[n] == 0 && continue
            idx_base = idx1 + ((n - 1) << 1)
            n00 += d[idx_base]
            n01 += d[idx_base + 1]
            n11 += 255 - d[idx_base] - d[idx_base + 1]
        end
    end
    n00 *= one_255th
    n01 *= one_255th
    n11 *= one_255th
    return hwe(n00, n01, n11)
end

@inline ccdf_chisq_1(x) = gamma_inc(convert(typeof(x), 1/2), x/2, 0)[2]
"""
    hwe(n00, n01, n11)
Hardy-Weinberg equilibrium test. `n00`, `n01`, `n11` are counts of homozygotes 
and heterozygoes respectively. Output is the p-value of type Float64.
"""
function hwe(n00::Real, n01::Real, n11::Real)
    n = n00 + n01 + n11
    n == 0 && return 1.0
    p0 = (n01 + 2n00) / 2n
    (p0 ≤ 0.0 || p0 ≥ 1.0) && return 1.0
    p1 = 1 - p0
    # Pearson's Chi-squared test
    e00 = n * p0 * p0
    e01 = 2n * p0 * p1
    e11 = n * p1 * p1
    ts = (n00 - e00)^2 / e00 + (n01 - e01)^2 / e01 + (n11 - e11)^2 / e11
    #pval = ccdf(Chisq(1), ts)
    pval = ccdf_chisq_1(ts)

    # TODO Fisher exact test
    return pval
end

function maf(b::Bgen, v::Variant; T=Float32, decompressed=nothing, rmask=nothing)
    io, h = b.io, b.header
    if length(v.genotypes) == 0 || length(v.genotypes[1].decompressed) == 0
        decompressed = decompress(io, v, h; decompressed=decompressed)
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

    maf(p, decompressed, idx, h.layout, rmask)
end

function maf(p::Preamble, d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8, 
    rmask::Union{Nothing, Vector{UInt16}})
    @assert layout == 2 "hwe only supported for layout 2"
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy 
    @assert length(p.missings) == 0 "current implementation does not allow missingness"
    idx1 = idx[1]
    # "counts" times 255.
    dosage_total = 0
    if p.n_samples >= 16
        @inbounds for n in 1:16:(p.n_samples - p.n_samples % 16)
            idx_base = idx1 + ((n-1) >> 1) << 2
            if rmask !== nothing
                rs = vload(Vec{16,UInt16}, rmask, n)
                if sum(rs) == 0
                    continue
                end
            end
            r = reinterpret(Vec{16, UInt16}, vload(Vec{32, UInt8}, d, idx_base))
            first  = (r & mask_odd)  << 1
            second = (r & mask_even) >> 8
            s = first + second
            if rmask !== nothing
                s = s * rs
            end
            dosage_total += sum(s)
        end
    end
    rem = p.n_samples % 16
    if rem != 0
        @inbounds for n in ((p.n_samples - rem) + 1) : p.n_samples
            rmask !== nothing && rmask[n] == 0 && continue
            idx_base = idx1 + ((n - 1) << 1)
            dosage_total += 2 * d[idx_base] + d[idx_base + 1]
        end
    end
    dosage_total *= one_255th
    dosage_total /= (rmask !== nothing ? sum(rmask) : p.n_samples)
    dosage_total < 1.0 ? dosage_total / 2 : 1 - dosage_total / 2
end
