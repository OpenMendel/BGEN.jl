"""
    Variant(b::Bgen, offset::Integer)
    Variant(io, offset, compression, layout, expected_n)
Parse information of a single variant beginning from `offset`.
"""
function Variant(io::IOStream, offset::Integer,
        compression::Integer, layout::Integer,
        expected_n::Integer)
    seek(io, offset)
    if eof(io)
        @error "reached end of file"
    end
    if layout == 1
        n_samples = read(io, UInt32)
    else
        n_samples = expected_n
    end
    if n_samples != expected_n
        @error "number of samples does not match"
    end
    varid_len = read(io, UInt16)
    varid = String(read(io, varid_len))
    rsid_len = read(io, UInt16)
    rsid = String(read(io, rsid_len))
    chrom_len = read(io, UInt16)
    chrom = String(read(io, chrom_len))
    pos = read(io, UInt32)

    if layout == 1
        n_alleles = 2
    else
        n_alleles = read(io, UInt16)
    end

    alleles = Array{String}(undef, n_alleles)
    for i in 1:n_alleles
        allele_len = read(io, UInt32)
        alleles_bytes = read(io, allele_len)
        alleles[i] = String(alleles_bytes)
    end

    if compression == 0 && layout == 1
        geno_block_size = 6 * n_samples
    else
        geno_block_size = read(io, UInt32)
    end
    geno_offset = position(io)
    next_var_offset = geno_offset + geno_block_size

    Variant(offset, geno_offset, next_var_offset, geno_block_size, n_samples,
        varid, rsid, chrom, pos, n_alleles, alleles, [])
end

function Variant(b::Bgen, offset::Integer)
    h = b.header
    if offset >= b.fsize
        @error "reached end of file"
    end
    Variant(b.io, offset, h.compression, h.layout, h.n_samples)
end

"""
    destroy_genotypes!(v::Variant)
Destroy any parsed genotype information.
"""
function clear!(v::Variant)
    resize!(v.genotypes, 0)
    return
end
