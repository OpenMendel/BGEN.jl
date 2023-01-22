"""
    filter(dest::AbstractString, b::Bgen, variant_mask::BitVector, sample_mask::BitVector;
    dest_sample=dest[1:end-5] * ".sample",
    sample_path=nothing, sample_names=b.samples,
    offsets=nothing, from_bgen_start=false)

Filters the input Bgen instance `b` based on `variant_mask` and `sample_mask`. The result
is saved in the new bgen file `dest`. Sample information is stored in `dest_sample`.
`sample_path` is the path of the `.sample` file for the input BGEN file, and 
`sample_names` stores the sample names in the BGEN file. 
`offsets` and `from_bgen_start` are arguments for the `iterator` function of `b`.

Only supports layout 2 and probibility bit depths should always be a multiple of 8. 
The output is always compressed in ZSTD. The sample names are stored in a separate .sample file, 
but not in the output .bgen file. 
"""
function filter(dest::AbstractString, b::Bgen, variant_mask::BitVector, 
    sample_mask::BitVector=trues(length(b.samples));
    dest_sample = dest[1:end-5] * ".sample",
    sample_path=nothing, sample_names=b.samples,
    offsets=nothing, from_bgen_start=false, use_zlib=false)

    @assert endswith(dest, ".bgen") "must use .bgen file"
    @assert b.header.layout == 2 "only layout 2 is supported."
    @assert length(variant_mask) == b.header.n_variants
    @assert length(sample_mask) == b.header.n_samples
    filter_samples(dest_sample, sample_mask; sample_path=sample_path, sample_names=sample_names)
    open(dest, "w") do io
        write(io, UInt32(20)) # offset, defaults to 20. 
        write(io, UInt32(20)) # length of header in bytes, defaults to 20.
        write(io, UInt32(sum(variant_mask))) # number of variants
        write(io, UInt32(sum(sample_mask))) # number of samples
        write(io, Vector{UInt8}("bgen")) # magic number
        if !use_zlib
            flag = 0x0000000a # zstd, layout 2, do not store sample info
        else
            flag = 0x00000009 # zlib, layout 2, do not store sample info
        end
        write(io, flag)
        v_it = iterator(b; offsets=offsets, from_bgen_start=from_bgen_start)
        for (i, v) in enumerate(v_it)
            if variant_mask[i]
                write_variant(io, b, v, sample_mask; use_zlib=use_zlib)
            end
        end

    end
end

function write_variant(io::IOStream, b::Bgen, v::Variant, sample_mask::BitVector; use_zlib=false)
    write(io, UInt16(length(v.varid)))  # length of varid
    write(io, v.varid)                  # varid
    write(io, UInt16(length(v.rsid)))   # length of rsid
    write(io, v.rsid)                   # rsid
    write(io, UInt16(length(v.chrom)))  # length of chrom
    write(io, v.chrom)                  # chrom
    write(io, v.pos)                    # position
    write(io, v.n_alleles)              # n of alleles
    for a_idx in 1:v.n_alleles
        write(io, UInt32(length(v.alleles[a_idx])))
        write(io, v.alleles[a_idx])
    end
    decompressed = decompress(b.io, v, b.header)
    if !all(sample_mask) # parse decompressed genotype, skip if all samples are chosen
        n_samples_new = sum(sample_mask)
        p = parse_preamble(decompressed, b.header, v)
        @assert p.bit_depth % 8 == 0 # probabilities should be byte aligned

        decompressed_new = Vector{UInt8}(undef, get_decompressed_length(p, decompressed, sample_mask))
        decompressed_new[1:4] = reinterpret(UInt8, UInt32[n_samples_new]) # update number of samples
        decompressed_new[5:6] .= decompressed[5:6] # number of alleles

        ploidy_old = @view decompressed[9 : 8 + p.n_samples]
        decompressed_new[9 : 8 + n_samples_new] = ploidy_old[sample_mask] # extract ploidies 
        ploidy_new = @view(decompressed_new[9 : 8 + n_samples_new]) .& 0x3f
        decompressed_new[7] = minimum(ploidy_new) # min ploidy
        decompressed_new[8] = maximum(ploidy_new) # max ploidy
        # phased flag
        decompressed_new[8 + n_samples_new + 1] = decompressed[8 + p.n_samples + 1] 
        # bit depth
        decompressed_new[8 + n_samples_new + 2] = decompressed[8 + p.n_samples + 2]
        offset = 10 + p.n_samples
        offset_new = 10 + n_samples_new
        base_bytes = p.bit_depth ÷ 8
        # write each genotype
        for (i, m) in enumerate(sample_mask)
            current_block_length = if p.phased == 1
                base_bytes * (ploidy_old[i] & 0x3f) * (p.n_alleles - 1)
            else
                z = ploidy_old[i] & 0x3f
                k = p.n_alleles
                base_bytes * (binomial(z + k - 1, k - 1) - 1)
            end
            if m
                decompressed_new[offset_new + 1 : offset_new + current_block_length] .=
                    decompressed[offset + 1 : offset + current_block_length]
                offset_new += current_block_length
            end
            offset += current_block_length
        end
        decompressed = decompressed_new
        @assert length(decompressed_new) == offset_new
    end
    if !use_zlib
        compressed = transcode(ZstdCompressor(), decompressed)
    else
        compressed = transcode(ZlibCompressor(), decompressed)
    end
    write(io, UInt32(length(compressed) + 4))
    write(io, UInt32(length(decompressed)))
    write(io, compressed)
end

function get_decompressed_length(p::Preamble, d::Vector{UInt8}, sample_mask::BitVector)
    cum_count = 10
    n_samples_new = sum(sample_mask)
    cum_count += n_samples_new
    base_bytes = p.bit_depth ÷ 8
    ploidy = d[9: 8 + p.n_samples]
    if p.phased == 1
        cum_count += base_bytes * sum(ploidy[sample_mask] .& 0x3f) * (p.n_alleles - 1)
    else
        for (i, m) in enumerate(sample_mask)
            if m
                z = ploidy[i] & 0x3f
                k = p.n_alleles
                cum_count += base_bytes * (binomial(z + k - 1, k - 1) - 1)
            end
        end
    end
    cum_count
end

function filter_samples(dest::AbstractString, sample_mask::BitVector; 
    sample_path=nothing, sample_names=nothing)
    io = open(dest, "w")
    if sample_path !== nothing
        sample_io = open(sample_path)
        println(io, readline(sample_io))
        println(io, readline(sample_io))
        for (i, sm) in enumerate(sample_mask)
            l = readline(sample_io)
            if sm
                println(io, l)
            end
        end
        close(sample_io)
    else
        @assert sample_names !== nothing "either `sample_path` or `sample_names` must be provided"
        println(io, "ID_1")
        println(io, "0")
        for (i, sm) in enumerate(sample_mask)
            if sm
                println(io, sample_names[i])
            end
        end
    end
    close(io)
end
