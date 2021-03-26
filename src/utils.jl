"""
    hwe(b::Bgen, v::Variant; T=Float32, decompressed=nothing)
    hwe(p::Preamble, d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8, 
        rmask::Union{Nothing, Vector{UInt16}})
Hardy-Weinberg equilibrium test for diploid biallelic case
"""
function hwe(p::Preamble, d::Vector{UInt8}, startidx::Integer, layout::UInt8, 
    rmask::Union{Nothing, Vector{UInt16}})
    @assert layout == 2 "hwe only supported for layout 2"
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy 
    @assert length(p.missings) == 0 "current implementation does not allow missingness"
    idx1 = startidx

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

"""
    maf(b::Bgen, v::Variant; T=Float32, decompressed=nothing)
    maf(p::Preamble, d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8, 
        rmask::Union{Nothing, Vector{UInt16}})
Minor-allele frequency for diploid biallelic case
"""
function maf(p::Preamble, d::Vector{UInt8}, startidx::Integer, layout::UInt8, 
    rmask::Union{Nothing, Vector{UInt16}})
    @assert layout == 2 "maf only supported for layout 2"
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy 
    @assert length(p.missings) == 0 "current implementation does not allow missingness"
    idx1 = startidx
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

"""
    info_score(b::Bgen, v::Variant; T=Float32, decompressed=nothing)
    info_score(p::Preamble, d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8, 
        rmask::Union{Nothing, Vector{UInt16}})
Information score of the variant.
"""
function info_score(p::Preamble, d::Vector{UInt8}, startidx::Integer, layout::UInt8, 
    rmask::Union{Nothing, Vector{UInt16}})
    @assert layout == 2 "info_score only supported for layout 2"
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy 
    @assert length(p.missings) == 0 "current implementation does not allow missingness"
    idx1 = startidx
    # "counts" times 255.
    samples_cum = 0
    mean_cum = 0.0
    sumsq_cum = 0.0
    dosage_sum = 0.0
    dosage_sumsq = 0.0
    rs = Vec{16,UInt16}((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
    rs_float = Vec{16, Float32}((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
    rs_float_ones = rs_float
    rs_floatarr = ones(Float32, 16)
    if p.n_samples >= 16
        @inbounds for n in 1:16:(p.n_samples - p.n_samples % 16)
            idx_base = idx1 + ((n-1) >> 1) << 2
            if rmask !== nothing
                rs = vload(Vec{16,UInt16}, rmask, n)
                if sum(rs) == 0
                    continue
                elseif sum(rs) == 16
                    rs_float = rs_float_ones
                else
                    rs_floatarr .= @view(rmask[n:n+15])
                    rs_float = vload(Vec{16, Float32}, rs_floatarr, 1)
                end
            end
            r = reinterpret(Vec{16, UInt16}, vload(Vec{32, UInt8}, d, idx_base))
            first  = (r & mask_odd)  << 1
            second = (r & mask_even) >> 8
            s = first + second
            if rmask !== nothing
                s = s * rs
            end
            dosage_float = one_255th * convert(
                Vec{16, Float32}, s)
            samples_new = sum(rs)
            sum_new = sum(dosage_float)
            samples_prev = samples_cum
            samples_cum = samples_cum + samples_new
            mean_new = sum_new / samples_new
            diff = dosage_float - mean_new
            if rmask !== nothing
                diff *= rs_float
            end
            sumsq_new = sum(diff ^ 2)
            delta = mean_new - mean_cum
            mean_cum += delta * samples_new / samples_cum
            sumsq_cum += sumsq_new + delta ^ 2 * samples_prev * samples_new / samples_cum
        end
    end
    rem = p.n_samples % 16
    if rem != 0
        @inbounds for n in ((p.n_samples - rem) + 1) : p.n_samples
            rmask !== nothing && rmask[n] == 0 && continue
            idx_base = idx1 + ((n - 1) << 1)
            dosage = one_255th * (2 * d[idx_base] + d[idx_base + 1])
            delta = dosage - mean_cum
            samples_cum += 1
            mean_cum += delta / samples_cum
            delta2 = dosage - mean_cum
            sumsq_cum += delta * delta2
        end
    end
    n_samples = rmask !== nothing ? sum(rmask) : p.n_samples
    p = mean_cum / 2
    v = sumsq_cum / (n_samples - 1)
    v / (2p * (1-p))
end

function counts!(p::Preamble, d::Vector{UInt8}, startidx::Integer, layout::UInt8, 
    rmask::Union{Nothing, Vector{UInt16}}; r::Union{Nothing,Vector{<:Integer}}=nothing, dosage::Bool=true)
    if dosage
        @assert layout == 2 "hwe only supported for layout 2"
        @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy 
        @assert length(p.missings) == 0 "current implementation does not allow missingness"
        idx1 = startidx
        if r !== nothing
            @assert length(r) == 512 
            fill!(r, 0)
        else
            r = zeros(UInt, 512)
        end
        if p.n_samples >= 16
            @inbounds for n in 1:16:(p.n_samples - p.n_samples % 16)
                idx_base = idx1 + ((n-1) >> 1) << 2
                if rmask !== nothing
                    rs = vload(Vec{16,UInt16}, rmask, n)
                    if sum(rs) == 0
                        continue
                    end
                end
                q = reinterpret(Vec{16, UInt16}, vload(Vec{32, UInt8}, d, idx_base))
                first  = (q & mask_odd)  << 1
                second = (q & mask_even) >> 8
                s = first + second
                if rmask !== nothing
                    s = s * rs
                end
                @inbounds for i in 1:16
                    if rs[i] != 0
                        ss = s[i]
                        r[ss + 1] += 1
                    end
                end
            end
        end
        rem = p.n_samples % 16
        if rem != 0
            @inbounds for n in ((p.n_samples - rem) + 1) : p.n_samples
                rmask !== nothing && rmask[n] == 0 && continue
                idx_base = idx1 + ((n - 1) << 1)
                r[2 * d[idx_base] + d[idx_base + 1]] += 1
            end
        end
    else
        @error "counts for non-dosage not implemented yet"
    end
    r
end

for ftn in [:maf, :hwe, :info_score, :counts!]
    @eval begin
        function $(ftn)(b::Bgen, v::Variant; T=Float32, decompressed=nothing, 
            is_decompressed=false, rmask=nothing, kwargs...)
            io, h = b.io, b.header
            if (decompressed !== nothing && !is_decompressed) || 
                (decompressed === nothing && (v.genotypes === nothing || 
                v.genotypes.decompressed === nothing))
                decompressed = decompress(io, v, h; decompressed=decompressed)
            else
                decompressed = v.genotypes.decompressed
            end
            startidx = 1
            if v.genotypes === nothing
                p = parse_preamble(decompressed, h, v)
                v.genotypes = Genotypes{T}(p, decompressed)
            else
                p = v.genotypes.preamble
            end
            if h.layout == 2
                startidx += 10 + h.n_samples
            end

            $(ftn)(p, decompressed, startidx, h.layout, rmask; kwargs...)
        end
    end
end