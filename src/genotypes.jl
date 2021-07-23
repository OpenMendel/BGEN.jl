const lookup = [i / 255 for i in 0:510]
@inline function unsafe_load_UInt64(v::Vector{UInt8}, i::Integer)
    p = convert(Ptr{UInt64}, pointer(v, i))
    unsafe_load(p)
end

"""
    Genotypes{T}(p::Preamble, d::Vector{UInt8}) where T <: AbstractFloat
Create `Genotypes` struct from the preamble and decompressed data string.
"""
function Genotypes{T}(p::Preamble, d::Vector{UInt8}) where T <: AbstractFloat
    Genotypes{T}(p, d, nothing, 0, nothing, false, false)
end

const zlib = ZlibDecompressor()

@inline function zstd_uncompress!(input::Vector{UInt8}, output::Vector{UInt8})
    r = ccall((:ZSTD_decompress, CodecZstd.libzstd),
        Csize_t, (Ptr{Cchar}, Cint, Ptr{Cchar}, Cint),
        pointer(output), length(output), pointer(input), length(input))
    @assert r == length(output) "zstd decompression returned data of wrong length"
end

function check_decompressed_length(io, v, h)
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
    return decompressed_length, decompressed_field
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
    decompressed_length, decompressed_field = check_decompressed_length(io, v, h)
    if decompressed !== nothing
        @assert length(decompressed) ==
            decompressed_length "decompressed length mismatch"
    else
        decompressed = Vector{UInt8}(undef, decompressed_length)
    end
    compressed_length = v.next_var_offset - v.geno_offset - decompressed_field
    buffer_compressed = read(io, compressed_length)
    if compression == 0
        decompressed .= buffer_compressed
    else
        if compression == 1
            codec = zlib
        elseif compression == 2
            codec = nothing
        else
            @error "invalid compression"
        end

        if compression == 1
            input = Buffer(buffer_compressed)
            output = Buffer(decompressed)
            error = Error()
            initialize(codec)
            _, _, e = process(codec, buffermem(input), buffermem(output), error)
            if e === :error
                throw(error[])
            end
            finalize(codec)
        elseif compression == 2
            zstd_uncompress!(buffer_compressed, decompressed)
        end
    end
    return decompressed
end

"""
    parse_ploidy(ploidy, d, idx, n_samples)
Parse ploidy part of the preamble.
"""
function parse_ploidy(ploidy::Union{UInt8,AbstractVector{UInt8}}, d::AbstractVector{UInt8},
    n_samples::Integer)
    missings = Int[]
    mask = 0x3f # 63 in UInt8
    mask_8 = 0x8080808080808080 # UInt64, mask for missingness
    idx1 = 9
    if typeof(ploidy) == UInt8 # if constant ploidy, just scan for missingness
        # check eight samples at a time

        if n_samples >= 8
            @inbounds for i in 0:8:(n_samples - (n_samples % 8) - 1)
                if mask_8 & unsafe_load_UInt64(d, idx1 + i) != 0
                    for j in (i+1):(i+8)
                        if d[idx1 + j - 1] & 0x80 != 0
                            push!(missings, j)
                        end
                    end
                end

            end
        end
        # remainder not in multiple of 8
        @inbounds for j in (n_samples - (n_samples % 8) + 1):n_samples
            if d[idx1 + j - 1] & 0x80 != 0
                push!(missings, j)
            end
        end
    else
        @inbounds for j in 1:n_samples
            ploidy[j] = mask & d[idx1 + j - 1]
            if d[idx1 + j - 1] & 0x80 != 0
                push!(missings, j)
            end
        end
    end
    return missings
end

@inline function get_max_probs(max_ploidy, n_alleles, phased)
    phased == 1 ?
        n_alleles : binomial(max_ploidy + n_alleles - 1, n_alleles - 1)
end

"""
    parse_preamble(d, idx, h, v)
Parse preamble of genotypes.
"""
function parse_preamble(d::AbstractVector{UInt8}, h::Header, v::Variant)
    startidx = 1
    if h.layout == 1
        n_samples = h.n_samples
        n_alleles = 2
        phased = false
        min_ploidy = 0x02
        max_ploidy = 0x02
        bit_depth = 16
    elseif h.layout == 2
        n_samples = reinterpret(UInt32, @view(d[startidx:startidx+3]))[1]
        startidx += 4
        @assert n_samples == h.n_samples "invalid number of samples"
        n_alleles = reinterpret(UInt16, @view(d[startidx:startidx+1]))[1]
        startidx += 2
        @assert n_alleles == v.n_alleles "invalid number of alleles"
        min_ploidy = d[startidx]
        startidx += 1
        max_ploidy = d[startidx]
        startidx += 1
    else
        @error "invalid layout"
    end
    constant_ploidy = (min_ploidy == max_ploidy)

    if constant_ploidy
        ploidy = max_ploidy
    else
        ploidy = Vector{UInt8}(undef, n_samples)
        fill!(ploidy, 0)
    end
    missings = []
    if h.layout == 2
        # this function also parses missingness.
        missings = parse_ploidy(ploidy, d, n_samples)
        startidx += n_samples
        phased = d[startidx]
        startidx += 1
        bit_depth = d[startidx]
        startidx += 1
    end
    max_probs = get_max_probs(max_ploidy, n_alleles, phased)
    Preamble(n_samples, n_alleles, phased, min_ploidy, max_ploidy, ploidy,
        bit_depth, max_probs, missings)
end

"""
    parse_layout1!(data, p, d, startidx)
Parse probabilities from layout 1.
"""
function parse_layout1!(data::AbstractArray{<:AbstractFloat},
    p::Preamble, d::AbstractArray{UInt8}, startidx::Integer
    )
    @assert length(data) == p.n_samples * p.max_probs
    factor = 1.0 / 32768
    idx = startidx
    @fastmath @inbounds for i in 1:p.max_probs:(p.n_samples * p.max_probs)
        j = idx
        data[i]     = reinterpret(UInt16, @view(d[j     : j + 1]))[1] * factor
        data[i + 1] = reinterpret(UInt16, @view(d[j + 2 : j + 3]))[1] * factor
        data[i + 2] = reinterpret(UInt16, @view(d[j + 4 : j + 5]))[1] * factor
        idx += 6

        # triple zero denotes missing for layout1
        if data[i] == 0.0 && data[i+1] == 0.0 && data[i+2] == 0.0
            data[i:i+2] .= NaN
            push!(p.missings, (i-1) ÷ p.max_probs + 1)
        end
    end
    return data
end

"""
   parse_layout2!(data, p, d, startidx)
Parse probabilities from layout 2.
"""
function parse_layout2!(data::AbstractArray{<:AbstractFloat},
    p::Preamble, d::AbstractArray{UInt8}, startidx::Integer
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
        idx1 = startidx
        @inbounds for offset in 1:3:(3 * nrows)
            idx2 = 2 * ((offset-1) ÷ 3)
            first = d[idx1 + idx2]
            second = d[idx1 + idx2 + 1]
            data[offset] = lookup[first + 1]
            data[offset + 1] = lookup[second + 1]
            data[offset + 2] = lookup[256 - first - second]
        end
    else
        idx1 = startidx
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
const one_255th = 1.0f0 / 255.0f0
const mask_odd = reinterpret(Vec{16, UInt16}, Vec{32, UInt8}(
    tuple(repeat([0xff, 0x00], 16)...)))
const mask_even = reinterpret(Vec{16, UInt16}, Vec{32, UInt8}(
    tuple(repeat([0x00, 0xff], 16)...)))
function ref_dosage_fast!(data::Vector{T}, p::Preamble,
    d::Vector{UInt8}, startidx::Integer, layout::UInt8
    ) where {T <:AbstractFloat}
    @assert length(data) == p.n_samples
    @assert layout == 2
    @assert p.bit_depth == 8 && p.max_probs == 3 && p.max_ploidy == p.min_ploidy
    idx1 = startidx

    if p.n_samples >= 16
        @inbounds for n in 1:16:(p.n_samples - p.n_samples % 16)
            idx_base = idx1 + ((n-1) >> 1) << 2
            r = reinterpret(Vec{16, UInt16}, vload(Vec{32, UInt8}, d, idx_base))
            second = (r & mask_even) >> 8
            first  = (r & mask_odd) << 1
            dosage_level = first + second
            dosage_level_float = one_255th * convert(
                Vec{16, T}, dosage_level)
            vstore(dosage_level_float, data, n)
        end
    end
    rem = p.n_samples % 16
    if rem != 0
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
    d::Vector{UInt8}, startidx::Integer, layout::UInt8
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
        j = startidx + bit_idx ÷ 8
        hom = (unsafe_load_UInt64(d, j) >> (bit_idx % 8)) & probs_mask
        bit_idx += p.bit_depth
        j = startidx + bit_idx ÷ 8
        het = (unsafe_load_UInt64(d, j) >> (bit_idx % 8)) & probs_mask
        bit_idx += p.bit_depth
        data[n] = ((hom * ploidy) + (het * half_ploidy)) * factor
        if layout == 1
            # layout 1 also stores hom_alt probability, and it indicates missing by
            # triple zero.
            j = startidx + bit_idx ÷ 8
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
    ref_dosage_phased!(data, p, d, idx, layout)
Dosage computation for phased genotypes.
"""
function ref_dosage_phased!(data::Vector{<:AbstractFloat}, p::Preamble,
    d::Vector{UInt8}, startidx::Integer, layout::UInt8
    )
    @assert length(data) == p.n_samples
    @assert layout == 2 "Phased genotypes not supported for Layout 1"
    ploidy = p.max_ploidy
    half_ploidy = ploidy / 2

    maxval = 2 ^ p.bit_depth - 1
    factor = 1.0 / maxval
    probs_mask = 0xFFFFFFFFFFFFFFFF >> (64 - p.bit_depth)
    bit_idx = 0
    for n = 1:p.n_samples
        if p.max_ploidy != p.min_ploidy
            ploidy = ploidy[n]
            half_ploidy = ploidy ÷ 2
        end
        ref_level = 0
        for _ in 1:ploidy
            j = startidx + bit_idx ÷ 8
            ref_level += (unsafe_load_UInt64(d, j) >> (bit_idx % 8)) & probs_mask
            bit_idx += p.bit_depth
        end
        data[n] = ref_level * factor
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

function _get_prob_matrix(d::Vector{T}, p::Preamble) where {T <: AbstractFloat}
    reshaped = reshape(d, p.max_probs, :)
    if p.phased == 1
        current = 1
        ragged = Matrix{T}(undef, p.max_ploidy * p.max_probs, p.n_samples)
        fill!(ragged, NaN)
        if p.max_ploidy == p.min_ploidy
            for i in 1:p.n_samples
                for j in 1:p.max_ploidy
                    first = (j-1) * p.max_probs + 1
                    last = j * p.max_probs
                    ragged[first:last, i] = reshaped[:, current]
                    current += 1
                end
            end
        else
            for (i, v) in enumerate(p.ploidy)
                for j in 1:v
                    first = (j-1) * p.max_probs + 1
                    last = j * p.max_probs
                    ragged[first:last, i] = reshaped[:, current]
                    current += 1
                end
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
The result is stored inside `v.genotypes.probs`, which can be cleared using
`clear!(v)`.

- T: type for the resutls
- `clear_decompressed`: clears decompressed byte string after execution if set `true`
"""
function probabilities!(b::Bgen, v::Variant;
        T=Float32, clear_decompressed=false, data=nothing, decompressed=nothing, is_decompressed=false)
    io, h = b.io, b.header
    if (decompressed !== nothing && !is_decompressed) ||
        (decompressed === nothing && (v.genotypes === nothing ||
        v.genotypes.decompressed === nothing))
        decompressed = decompress(io, v, h; decompressed=decompressed)
    else
        decompressed = v.genotypes.decompressed
    end
    if v.genotypes === nothing
        p = parse_preamble(decompressed, h, v)
        v.genotypes = Genotypes{T}(p, decompressed)
    else
        p = v.genotypes.preamble
    end
    startidx = 1
    if h.layout == 2
        startidx += 10 + h.n_samples
    end

    genotypes = v.genotypes
    data_size = get_data_size(p, h.layout)

    # skip parsing if already parsed
    if genotypes.probs !== nothing && length(genotypes.probs) >= data_size
        _get_prob_matrix(genotypes.probs, p)
    end

    if data !== nothing
        @assert length(data) == data_size
        genotypes.probs = data
    else
        genotypes.probs = Vector{T}(undef, data_size)
    end
    if h.layout == 1
        parse_layout1!(genotypes.probs, p, decompressed, startidx)
    elseif h.layout == 2
        parse_layout2!(genotypes.probs, p, decompressed, startidx)
    end
    if clear_decompressed
        clear_decompressed!(genotypes)
    end
    return _get_prob_matrix(genotypes.probs, p)
end

function ref_allele_dosage!(b::Bgen, v::Variant;
        T=Float32, mean_impute=false, clear_decompressed=false,
        data=nothing, decompressed=nothing, is_decompressed=false)
    io, h = b.io, b.header
    # just return it if already computed
    if v.genotypes !== nothing && v.genotypes.dose !== nothing
        genotypes = v.genotypes
        p = genotypes.preamble
        if genotypes.dose_mean_imputed && !mean_impute
            genotypes.dose[p.missings] .= NaN
        elseif !genotypes.dose_mean_imputed && mean_impute
            genotypes.dose[p.missings] .= mean(filter(!isnan, genotypes.dose))
        end
        if genotypes.minor_allele_dosage && genotypes.minor_idx != 1
            alt_dosage!(genotypes.dose, p)
        end
        return v.genotypes.dose
    end
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

    @assert p.n_alleles == 2 "allele dosages are available for biallelic variants"
    #@assert p.phased == 0

    genotypes = v.genotypes
    if data !== nothing
        @assert length(data) == h.n_samples
    else
        genotypes.dose = Vector{T}(undef, h.n_samples)
        data = genotypes.dose
    end

    if p.phased == 0
        if p.max_ploidy == p.min_ploidy && p.max_probs == 3 && p.bit_depth == 8 &&
                b.header.layout == 2
            ref_dosage_fast!(data, p, decompressed, startidx, h.layout)
        else
            ref_dosage_slow!(data, p, decompressed, startidx, h.layout)
        end
    else # phased
        ref_dosage_phased!(data, p, decompressed, startidx, h.layout)
    end
    genotypes.minor_idx = find_minor_allele(data, p)
    data[p.missings] .= NaN
    if mean_impute
        data[p.missings] .= mean(filter(!isnan, data))
    end
    if clear_decompressed
        clear_decompressed!(genotypes)
    end
    return data
end

"""
    minor_allele_dosage!(b::Bgen, v::Variant; T=Float32,
    mean_impute=false, clear_decompressed=false)
Given a `Bgen` struct and a `Variant`, compute minor allele dosage.
The result is stored inside `v.genotypes.dose`, which can be cleared using
`clear!(v)`.

- `T`: type for the results
- `mean_impute`: impute missing values with the mean of nonmissing values
- `clear_decompressed`: clears decompressed byte string after execution if set `true`
"""
function minor_allele_dosage!(b::Bgen, v::Variant;
        T=Float32, mean_impute=false, clear_decompressed=false,
        data=nothing, decompressed=nothing, is_decompressed=false)
    # just return it if already computed
    io, h = b.io, b.header
    if v.genotypes !== nothing && v.genotypes.dose !== nothing
        genotypes = v.genotypes
        p = genotypes.preamble
        if genotypes.dose_mean_imputed && !mean_impute
            genotypes.dose[p.missings] .= NaN
        elseif !genotypes.dose_mean_imputed && mean_impute
            genotypes.dose[p.missings] .= mean(filter(!isnan, genotypes.dose))
        end
        if !genotypes.minor_allele_dosage && genotypes.minor_idx != 1
            alt_dosage!(genotypes.dose, p)
        end
        return v.genotypes.dose
    end
    ref_allele_dosage!(b, v; T=T, mean_impute=mean_impute,
        clear_decompressed=clear_decompressed, data=data, decompressed=decompressed,
        is_decompressed=is_decompressed)
    genotypes = v.genotypes
    if data === nothing
        data = genotypes.dose
    end
    if genotypes.minor_idx != 1
        alt_dosage!(data, genotypes.preamble)
    end
    return data
end

"""
    clear!(g::Genotypes)
    clear!(v::Variant)
Clears cached decompressed byte representation, probabilities, and dose.
If `Variant` is given, it removes the corresponding `.genotypes` altogether.
"""
function clear!(g::Genotypes)
    g.decompressed = nothing
    g.probs = nothing
    g.dose = nothing
    return
end

"""
    clear_decompressed!(g::Genotypes)
Clears cached decompressed byte representation.
"""
function clear_decompressed!(g::Genotypes)
    g.decompressed = nothing
    return
end
