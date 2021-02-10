struct Header
    offset::UInt32
    header_length::UInt32
    n_variants::UInt32
    n_samples::UInt32
    compression::UInt8
    layout::UInt8
    has_sample_ids::Bool
end

function Header(io::IOStream)
    seek(io, 0)
    offset = read(io, UInt32)
    header_length = read(io, UInt32)
    n_variants = read(io, UInt32)
    n_samples = read(io, UInt32)
    magic = read(io, 4)
    # check magic number
    @assert String(magic) == "bgen" || all(magic .== 0) "Magic number mismatch"
    seek(io, header_length)
    flags = read(io, UInt32)
    compression = flags & 0x03
    layout = (flags & 0x3c) >> 2
    has_sample_ids = convert(Bool, flags & 0x80000000 >> 31)
    Header(offset, header_length, n_variants, n_samples, compression, layout,
        has_sample_ids)
end

function Header(filename::String)
    io = open(filename)
    Header(io)
end

@inline function offset_first_variant(h::Header)
    return h.offset + 4
end
