const lookup = [i / 255 for i in 0:510]
using SIMD
@inline function unsafe_load_UInt64(v::Vector{UInt8}, i::Integer)
    p = convert(Ptr{UInt64}, pointer(v, i))
    unsafe_load(p)
end

"""
    Genotypes{T}(p::Preamble, d::Vector{UInt8}) where T <: AbstractFloat
Create `Genotypes` struct from the preamble and decompressed data string.
"""
function Genotypes{T}(p::Preamble, d::Vector{UInt8}) where T <: AbstractFloat
    Genotypes{T}(p, d, T[], UInt8[0], T[])
end

"""
    decompress(io, v, h; decompressed=nothing)
Decompress the compressed byte string for genotypes.
"""
function decompress(io::IOStream, v::Variant, h::Header;
    decompressed::Union{Nothing, AbstractVector{UInt8}}=nothing
    )
    compression = h.compression
    seek(io, v.geno_offset)
    decompressed_field = 0
    if h.compression != 0
        if h.layout == 1
            decompressed_length = 6 * h.n_samples
        elseif h.layout == 2
            decompressed_field = 4
            decompressed_length = read(io, UInt32)
        end
    end
    if decompressed !== nothing
        @assert length(decompressed) ==
            decompressed_length "decompressed length mismatch"
    end
    compressed_length = v.next_var_offset - v.geno_offset - decompressed_field
    buffer_compressed = read(io, compressed_length)
    if compression == 0
        decompressed = buffer_compressed
    elseif compression == 1
        decompressed = transcode(ZlibDecompressor,
            (buffer_compressed))
    elseif compression == 2
        decompressed = transcode(ZstdDecompressor,
            (buffer_compressed))
    else
        @error "invalid compression"
    end
    return decompressed
end

"""
    parse_ploidy!(ploidy, d, idx, n_samples)
Parse ploidy part of the preamble.
"""
function parse_ploidy!(ploidy::AbstractVector{UInt8}, d::AbstractVector{UInt8},
    idx::Vector{<:Integer}, n_samples::Integer)
    missings = Int[]
    mask = 0x3f # 63 in UInt8
    mask_8 = 0x8080808080808080 # UInt64, mask for missingness
    if ploidy[1] != 0 # if constant ploidy, just scan for missingness
        # check eight samples at a time
        idx1 = idx[1]
        if n_samples >= 8
            @inbounds for i in 0:8:(n_samples - (n_samples % 8) - 1)
                if mask_8 & unsafe_load_UInt64(d, idx1 + i) != 0
                    for j in (i+1):(i+8)
                        if d[idx[1] + j - 1] & 0x80 != 0
                            push!(missings, j)
                        end
                    end
                end

            end
        end
        # remainder not in multiple of 8
        @inbounds for j in (n_samples - (n_samples % 8) + 1):n_samples
            if d[idx[1] + j - 1] & 0x80 != 0
                push!(missings, j)
            end
        end
    else
        @inbounds for j in 1:n_samples
            ploidy[j] = mask & d[idx[1] + j - 1]
            if d[idx[1] + j - 1] & 0x80 != 0
                push!(missings, j)
            end
        end
    end
    idx[1] += n_samples
    return missings
end

@inline function get_max_probs(max_ploidy, n_alleles, phased)
    phased == 1 ?
        n_alleles : binomial(max_ploidy + n_alleles - 1, n_alleles - 1)
end

"""
    parse_preamble!(d, idx, h, v)
Parse preamble of genotypes.
"""
function parse_preamble!(d::AbstractVector{UInt8}, idx::Vector{<:Integer},
        h::Header, v::Variant)
    if h.layout == 1
        n_samples = h.n_samples
        n_alleles = 2
        phased = false
        min_ploidy = 2
        max_ploidy = 2
        bit_depth = 16
    elseif h.layout == 2
        n_samples = reinterpret(UInt32, @view(d[idx[1]:idx[1]+3]))[1]
        idx[1] += 4
        @assert n_samples == h.n_samples "invalid number of samples"
        n_alleles = reinterpret(UInt16, @view(d[idx[1]:idx[1]+1]))[1]
        idx[1] += 2
        @assert n_alleles == v.n_alleles "invalid number of alleles"
        min_ploidy = d[idx[1]]
        idx[1] += 1
        max_ploidy = d[idx[1]]
        idx[1] += 1
    else
        @error "invalid layout"
    end
    constant_ploidy = (min_ploidy == max_ploidy)
    ploidy = Vector{UInt8}(undef, n_samples)
    if constant_ploidy
        fill!(ploidy, max_ploidy)
    else
        fill!(ploidy, 0)
    end
    missings = []
    if h.layout == 2
        # this function also parses missingness.
        missings = parse_ploidy!(ploidy, d, idx, n_samples)
        phased = d[idx[1]]
        idx[1] += 1
        bit_depth = d[idx[1]]
        idx[1] += 1
    end
    max_probs = get_max_probs(max_ploidy, n_alleles, phased)
    Preamble(n_samples, n_alleles, phased, min_ploidy, max_ploidy, ploidy,
        bit_depth, max_probs, missings)
end

"""
    parse_layout1!(data, p, d, idx)
Parse probabilities from layout 1.
"""
function parse_layout1!(data::AbstractArray{<:AbstractFloat},
    p::Preamble, d::AbstractArray{UInt8}, idx::Vector{<:Integer}
    )
    @assert length(data) == p.n_samples * p.max_probs
    factor = 1.0 / 32768
    idx1 = idx[1]
    @fastmath @inbounds for i in 1:p.max_probs:(p.n_samples * p.max_probs)
        j = idx1
        data[i]     = reinterpret(UInt16, @view(d[j     : j + 1]))[1] * factor
        data[i + 1] = reinterpret(UInt16, @view(d[j + 2 : j + 3]))[1] * factor
        data[i + 2] = reinterpret(UInt16, @view(d[j + 4 : j + 5]))[1] * factor
        idx1 += 6

        # triple zero denotes missing for layout1
        if data[i] == 0.0 && data[i+1] == 0.0 && data[i+2] == 0.0
            data[i:i+2] .= NaN
            push!(p.missings, (i-1) ÷ p.max_probs + 1)
        end
    end
    return data
end

"""
   parse_layout2!(data, p, d, idx)
Parse probabilities from layout 2.
"""
function parse_layout2!(data::AbstractArray{<:AbstractFloat},
    p::Preamble, d::AbstractArray{UInt8}, idx::Vector{<:Integer}
    )
    constant_ploidy = p.max_ploidy == p.min_ploidy
    if p.phased == 0
        nrows = p.n_samples
    else
        if constant_ploidy
            nrows = p.n_samples * p.max_ploidy
        else
            nrows = sum(p.ploidy)
        end
    end
    @assert length(data) == p.max_probs * nrows

    max_less_1 = p.max_probs - 1
    prob = 0.0
    factor = 1.0 / (2 ^ p.bit_depth - 1)
    # mask for depth not multiple of 8
    probs_mask = 0xFFFFFFFFFFFFFFFF >> (64 - p.bit_depth)
    bit_idx = 0

    if constant_ploidy && p.max_probs == 3 && p.bit_depth == 8
        # fast path for unphased, ploidy==2, 8 bits per prob.
        idx1 = idx[1]
        @inbounds for offset in 1:3:(3 * nrows)
            idx2 = 2 * ((offset-1) ÷ 3)
            first = d[idx1 + idx2]
            second = d[idx1 + idx2 + 1]
            data[offset] = lookup[first + 1]
            data[offset + 1] = lookup[second + 1]
            data[offset + 2] = lookup[256 - first - second]
        end
    else
        idx1 = idx[1]
        @inbounds for offset in 1:p.max_probs:(nrows * p.max_probs)
            # number of probabilities to be read per row
            if constant_ploidy
                n_probs = max_less_1
            elseif p.phased == 1
                n_probs = p.n_alleles - 1
            elseif p.ploidy[offset ÷ p.max_probs + 1] == 2 && p.n_alleles == 2
                n_probs = 2
            else
                n_probs = binomial(p.ploidy[offset ÷ p.max_probs + 1] +
                    p.n_alleles - 1, p.n_alleles - 1) - 1
            end

            remainder = 1.0
            @inbounds for i in 1:n_probs
                j = idx1 + bit_idx ÷ 8
                @inbounds prob = (unsafe_load_UInt64(d, j) >> (bit_idx % 8) &
                    probs_mask) * factor
                bit_idx += p.bit_depth
                remainder -= prob
                data[offset + i - 1] = prob
            end
            data[offset + n_probs] = remainder
            if n_probs + 1 < p.max_probs
                data[(offset + n_probs + 1):(offset + p.max_probs - 1)] .= NaN
            end
        end
    end
    for m in p.missings
        offset = p.max_probs * (m - 1) + 1
        data[offset:(offset + p.max_probs - 1)] .= NaN
    end
    return data
end

"""
    ref_dosage_fast!(data, p, d, idx, layout)
Dosage retrieval for 8-bit biallele case, no floating-point operations!
"""
function ref_dosage_fast!(data::Vector{<:AbstractFloat}, p::Preamble,
    d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8
    )
    @assert layout == 2
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy
    idx1 = idx[1]
    if p.n_samples >= 2
        @inbounds for n in 1:2:(p.n_samples - p.n_samples % 2)
            idx_base = idx1 + ((n-1) >> 1) << 2
            data[n] = lookup[d[idx_base] * 2 + 
                                d[idx_base + 1] + 1]
            data[n+1] = lookup[d[idx_base + 2] * 2 + 
                                d[idx_base + 3] + 1]
        end
    end
    if p.n_samples % 2 == 1
        idx_base = idx1 + ((p.n_samples - 1) << 1)
        data[p.n_samples] = lookup[d[idx_base] * 2 +
                            d[idx_base + 1] + 1]
    end
    return data
end

const  one_255th = 1.0f0 / 255.0f0
function ref_dosage_fast2!(data::Vector{Float32}, p::Preamble,
    d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8
    )
    @assert layout == 2
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy
    idx1 = idx[1]
    mask_odd = reinterpret(Vec{16, UInt16}, Vec{32, UInt8}(
        tuple(repeat([0xff, 0x00], 16)...)))
    mask_even = reinterpret(Vec{16, UInt16}, Vec{32, UInt8}(
        tuple(repeat([0x00, 0xff], 16)...)))

    if p.n_samples >= 16
        @inbounds for n in 1:16:(p.n_samples - p.n_samples % 8)
            idx_base = idx1 + ((n-1) >> 1) << 2
            r = reinterpret(Vec{16, UInt16}, vload(Vec{32, UInt8}, d, idx_base))
            second = (r & mask_even) >> 8
            first  = (r & mask_odd) << 1
            dosage_level = first + second
            dosage_level_float = one_255th * convert(
                Vec{16, Float32}, dosage_level)
            vstore(dosage_level_float, data, n)
        end
    end
    rem = p.n_samples % 16
    if p.n_samples % 16 != 0
        @inbounds for n in ((p.n_samples - rem) + 1) : p.n_samples
            idx_base = idx1 + ((n - 1) << 1)
            data[n] = lookup[d[idx_base] * 2 +
                                d[idx_base + 1] + 1]
        end
    end
    return data
end
"""
    ref_dosage_slow!(data, p, d, idx, layout)
Dosage computation for general case.
"""
function ref_dosage_slow!(data::Vector{<:AbstractFloat}, p::Preamble,
    d::Vector{UInt8}, idx::Vector{<:Integer}, layout::UInt8
    )
    @assert length(data) == p.n_samples
    ploidy = p.max_ploidy
    half_ploidy = ploidy / 2

    maxval = 2 ^ p.bit_depth - 1
    factor = layout == 2 ? 1.0 / maxval : 1.0 / 32768
    probs_mask = 0xFFFFFFFFFFFFFFFF >> (64 - p.bit_depth)
    bit_idx = 0
    for n = 1:p.n_samples
        if p.max_ploidy != p.min_ploidy
            ploidy = ploidy[n]
            half_ploidy = ploidy ÷ 2
        end
        j = idx[1] + bit_idx ÷ 8
        hom = (unsafe_load_UInt64(d, j) >> (bit_idx % 8)) & probs_mask
        bit_idx += p.bit_depth
        j = idx[1] + bit_idx ÷ 8
        het = (unsafe_load_UInt64(d, j) >> (bit_idx % 8)) & probs_mask
        bit_idx += p.bit_depth
        data[n] = ((hom * ploidy) + (het * half_ploidy)) * factor
        if layout == 1
            # layout 1 also stores hom_alt probability, and it indicates missing by
            # triple zero.
            j = idx[1] + bit_idx ÷ 8
            hom_alt = (unsafe_load_UInt64(d, j) >> (bit_idx % 8)) & probs_mask
            bit_idx += p.bit_depth
            if hom == 0 && het == 0 && hom_alt == 0
                push!(p.missings, n)
            end
        end
    end
    return data
end

"""
    alt_dosage(data, p)
Switch ref allele dosage `data` to alt allele dosage.
"""
function alt_dosage!(data::Vector{<:AbstractFloat}, p::Preamble)
    if p.n_samples >= 8
        @inbounds for n in 1:8:(p.n_samples - p.n_samples % 8)
            data[n]     = 2.0 - data[n]
            data[n + 1] = 2.0 - data[n + 1]
            data[n + 2] = 2.0 - data[n + 2]
            data[n + 3] = 2.0 - data[n + 3]
            data[n + 4] = 2.0 - data[n + 4]
            data[n + 5] = 2.0 - data[n + 5]
            data[n + 6] = 2.0 - data[n + 6]
            data[n + 7] = 2.0 - data[n + 7]
        end
    end
    @inbounds for n in (p.n_samples - p.n_samples % 8 + 1):p.n_samples
        data[n] = 2.0 - data[n]
    end
end

"""
    find_minor_allele(data, p)
Find minor allele index, returns 1 (ref) or 2 (alt)
"""
function find_minor_allele(data::Vector{<:AbstractFloat}, p::Preamble)
    batchsize = 100
    increment = max(p.n_samples ÷ batchsize, 1)
    total = 0.0
    freq = 0.0
    cnt = 0
    for idx2 in 1:increment
        for n in idx2:increment:p.n_samples
            cnt += 1
            total += data[n]
        end
        freq = total / (cnt * 2)
        @assert 0 <= freq <= 1
        if minor_certain(freq, batchsize * idx2, 5.0)
            break
        end
    end
    if freq <= 0.5
        return 1
    else
        return 2
    end
end

@inline function get_data_size(p::Preamble, layout::Integer)
    if layout == 1
        return p.n_samples * p.max_probs
    else
        constant_ploidy = p.max_ploidy == p.min_ploidy
        if p.phased == 0
            nrows = p.n_samples
        else
            if constant_ploidy
                nrows = p.n_samples * p.max_ploidy
            else
                nrows = sum(p.ploidy)
            end
        end
        return p.max_probs * nrows
    end
end

function _get_prob_matrix(d::Vector{T}, p::Preamble) where T <: AbstractFloat
    reshaped = reshape(d, p.max_probs, :)
    if p.phased == 1
        current = 1
        ragged = Matrix{T}(undef, p.max_ploidy * p.max_probs, p.n_samples)
        fill!(ragged, NaN)
        for (i, v) in enumerate(p.ploidy)
            for j in 1:v
                first = (j-1) * p.max_probs + 1
                last = j * p.max_probs
                ragged[first:last, i] = reshaped[:, current]
                current += 1
            end
        end
        return ragged
    else
        return reshaped
    end
end

"""
    probabilities!(b::Bgen, v::Variant; T=Float32, clear_decompressed=false)
Given a `Bgen` struct and a `Variant`, compute probabilities.
The result is stored inside `v.genotypes[1].probs`, which can be cleared using
`clear!(v)`.

- T: type for the resutls
- `clear_decompressed`: clears decompressed byte string after execution if set `true`
"""
function probabilities!(b::Bgen, v::Variant;
        T=Float32, clear_decompressed=false)
    io, h = b.io, b.header
    if length(v.genotypes) == 0 || length(v.genotypes[1].decompressed) == 0
        decompressed = decompress(io, v, h)
    else
        decompressed = v.genotypes[1].decompressed
    end
    idx = [1]
    if length(v.genotypes) == 0
        p = parse_preamble!(decompressed, idx, h, v)
        push!(v.genotypes, Genotypes{T}(p, decompressed))
    else
        if h.layout == 2
            idx[1] += 10 + h.n_samples
        end
        p = v.genotypes[1].preamble
    end

    genotypes = v.genotypes[1]
    data_size = get_data_size(p, h.layout)

    # skip parsing if already parsed
    if length(genotypes.probs) >= data_size
        _get_prob_matrix(genotypes.probs, p)
    end
    resize!(genotypes.probs, data_size)
    if h.layout == 1
        parse_layout1!(genotypes.probs, p, decompressed, idx)
    elseif h.layout == 2
        parse_layout2!(genotypes.probs, p, decompressed, idx)
    end
    if clear_decompressed
        clear_decompressed!(genotypes)
    end
    return _get_prob_matrix(genotypes.probs, p)
end

"""
    minor_allele_dosage!(b::Bgen, v::Variant; T=Float32,
    mean_impute=false, clear_decompressed=false)
Given a `Bgen` struct and a `Variant`, compute minor allele dosage.
The result is stored inside `v.genotypes[1].dose`, which can be cleared using
`clear!(v)`.

- `T`: type for the results
- `mean_impute`: impute missing values with the mean of nonmissing values
- `clear_decompressed`: clears decompressed byte string after execution if set `true`
"""
function minor_allele_dosage!(b::Bgen, v::Variant;
        T=Float32, mean_impute=false, clear_decompressed=false)
    io, h = b.io, b.header
    # just return it if already computed
    if length(v.genotypes) == 1 && length(v.genotypes[1].dose) == h.n_samples
        genotypes = v.genotypes[1]
        p = genotypes.preamble
        genotypes.dose[p.missings] .= NaN
        if mean_impute
            genotypes.dose[p.missings] .= mean(filter(!isnan, genotypes.dose))
        end
        return v.genotypes[1].dose
    end
    if length(v.genotypes) == 0 || length(v.genotypes[1].decompressed) == 0
        decompressed = decompress(io, v, h)
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

    @assert p.phased == 0
    @assert p.n_alleles == 2 "allele dosages are available for non-biallelic var"

    genotypes = v.genotypes[1]

    resize!(genotypes.dose, h.n_samples)
    if p.max_ploidy == p.min_ploidy && p.max_probs == 3 && p.bit_depth == 8 &&
            b.header.layout == 2
        ref_dosage_fast2!(genotypes.dose, p, decompressed, idx, h.layout)
    else
        ref_dosage_slow!(genotypes.dose, p, decompressed, idx, h.layout)
    end

    genotypes.minor_idx[1] = find_minor_allele(genotypes.dose, p)
    if genotypes.minor_idx[1] != 1
        alt_dosage!(genotypes.dose, p)
    end

    genotypes.dose[p.missings] .= NaN
    if mean_impute
        genotypes.dose[p.missings] .= mean(filter(!isnan, genotypes.dose))
    end
    if clear_decompressed
        clear_decompressed!(genotypes)
    end
    return genotypes.dose
end

"""
    clear!(g::Genotypes)
    clear!(v::Variant)
Clears cached decompressed byte representation, probabilities, and dose.
If `Variant` is given, it removes the corresponding `.genotypes` altogether.
"""
function clear!(g::Genotypes)
    resize!(g.decompressed, 0)
    resize!(g.probs, 0)
    resize!(g.dose, 0)
    return
end

"""
    clear_decompressed!(g::Genotypes)
Clears cached decompressed byte representation.
"""
function clear_decompressed!(g::Genotypes)
    resize!(g.decompressed, 0)
    return
end
