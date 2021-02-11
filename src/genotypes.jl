const lookup = [i / 255 for i in 0:510]

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
        decompressed = transcode(GzipDecompressor,
            (buffer_compressed[1:compressed_length]))
    elseif compression == 2
        decompressed = transcode(ZstdDecompressor,
            (buffer_compressed[1:compressed_length]))
    else
        @error "invalid compression"
    end
    return decompressed
end

struct Preamble
    n_samples::UInt32
    n_alleles::UInt16
    phased::UInt8
    min_ploidy::UInt8
    max_ploidy::UInt8
    ploidy::Vector{UInt8}
    bit_depth::UInt8
    max_probs::Int
    missings::Vector{Int}
end

function parse_ploidy!(ploidy::AbstractVector{UInt8}, d::AbstractVector{UInt8},
    idx::Vector{<:Integer})
    n_samples = length(ploidy)
    missings = Int[]
    mask = 0x3f # 63 in UInt8
    mask_8 = 0x8080808080808080 # UInt64, mask for missingness
    if ploidy[1] != 0 # if constant ploidy, just scan for missingness
        # check eight samples at a time
        @inbounds for i in 0:8:(n_samples - (n_samples % 8) - 1)
            if mask_8 & reinterpret(UInt64,
                @view(d[idx[1]+i:(idx[1]+i+7)]))[1] != 0
                for j in (i+1):(i+8)
                    if d[idx[1] + j - 1] & 0x80 != 0
                        push!(missings, j)
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
        missings = parse_ploidy!(ploidy, d, idx)
        phased = d[idx[1]]
        idx[1] += 1
        bit_depth = d[idx[1]]
        idx[1] += 1
    end
    max_probs = get_max_probs(max_ploidy, n_alleles, phased)
    Preamble(n_samples, n_alleles, phased, min_ploidy, max_ploidy, ploidy,
        bit_depth, max_probs, missings)
end

function parse_layout1!(data::AbstractArray{<:AbstractFloat},
    p::Preamble, d::AbstractArray{UInt8}, idx::Vector{<:Integer}
    )
    @assert length(data) == p.n_samples * p.max_probs
    factor = 1.0 / 32768
    @fastmath @inbounds @simd for i in 1:p.max_probs:(p.n_samples * p.max_probs)
        j = idx[1]
        data[i]     = reinterpret(UInt16, @view(d[j     : j + 1]))[1] * factor
        data[i + 1] = reinterpret(UInt16, @view(d[j + 2 : j + 3]))[1] * factor
        data[i + 2] = reinterpret(UInt16, @view(d[j + 4 : j + 5]))[1] * factor
        idx[1] += 6

        # triple zero denotes missing for layout1
        if data[i] == 0.0 && data[i+1] == 0.0 && data[i+2] == 0.0
            data[i:i+2] .= NaN
        end
    end
    return data
end

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
        idx2 = 0
        # fast path for unphased, ploidy==2, 8 bits per prob.
        @inbounds for offset in 1:3:(3 * nrows)
            first = d[idx[1] + idx2]
            second = d[idx[1] + idx2 + 1]
            data[offset] = lookup[first + 1]
            data[offset + 1] = lookup[second + 1]
            data[offset + 2] = lookup[256 - first - second]
            idx2 += 2
        end
    else
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
            @fastmath @inbounds for i = 1:n_probs
                j = idx[1] + bit_idx ÷ 8
                prob = ((reinterpret(UInt64, d[j:j+7])[1] >> (bit_idx % 8)) &
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

function parse_genotype!(data::AbstractArray{<:AbstractFloat},
    io::IOStream, h::Header, v::Variant;
    decompressed::Union{Nothing, AbstractVector{UInt8}}=nothing
    )
    if decompressed === nothing
        decompressed = decompress(io, v, h)
    end
    idx = [1]
    p = parse_preamble!(decompressed, idx, h, v)
    data_size = get_data_size(p, h.layout)
    @assert length(data) == data_size "incorrect length of data"
    if h.layout == 1
        parse_layout1!(data, p, decompressed, idx)
    elseif h.layout == 2
        parse_layout2!(data, p, decompressed, idx)
    end
    return p, data
end

function parse_genotype!(data::AbstractArray{<:AbstractFloat},
    b::Bgen, v::Variant;
    decompressed::Union{Nothing, AbstractVector{UInt8}}=nothing
    )
    parse_genotype!(data, b.io, b.header, v; decompressed=decompressed)
end

function parse_genotype(io::IOStream, h::Header, v::Variant;
    decompressed::Union{Nothing, AbstractVector{UInt8}}=nothing
    )
    if decompressed === nothing
        decompressed = decompress(io, v, h)
    end
    idx = [1]
    p = parse_preamble!(decompressed, idx, h, v)
    println(p)
    data_size = get_data_size(p, h.layout)
    data = Vector{Float64}(undef, data_size)
    if h.layout == 1
        parse_layout1!(data, p, decompressed, idx)
    elseif h.layout == 2
        parse_layout2!(data, p, decompressed, idx)
    end
    return p, data
end

function parse_genotype(b::Bgen, v::Variant;
    decompressed::Union{Nothing, AbstractVector{UInt8}}=nothing)
    parse_genotype(b.io, b.header, v; decompressed=decompressed)
end
