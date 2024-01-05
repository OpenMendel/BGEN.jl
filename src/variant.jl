"""
    BgenVariant(b::Bgen, offset::Integer)
    BgenVariant(io, offset, compression, layout, expected_n)
Parse information of a single variant beginning from `offset`.
"""
function BgenVariant(io::IOStream, offset::Integer,
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

    BgenVariant(offset, geno_offset, next_var_offset, geno_block_size, n_samples,
        varid, rsid, chrom, pos, n_alleles, alleles, nothing)
end

function BgenVariant(b::Bgen, offset::Integer)
    h = b.header
    if offset >= b.fsize
        @error "reached end of file"
    end
    BgenVariant(b.io, offset, h.compression, h.layout, h.n_samples)
end

@inline n_samples(v::BgenVariant)::Int = v.n_samples
@inline varid(v::BgenVariant) = v.varid
@inline rsid(v::BgenVariant) = v.rsid
@inline chrom(v::BgenVariant) = v.chrom
@inline pos(v::BgenVariant)::Int = v.pos
@inline n_alleles(v::BgenVariant)::Int = v.n_alleles
@inline alleles(v::BgenVariant) = v.alleles

# The following functions are valid only after calling `probabilities!()`
# or `minor_allele_dosage!()`
@inline phased(v::BgenVariant) = v.genotypes.preamble.phased
@inline min_ploidy(v::BgenVariant) = v.genotypes.preamble.min_ploidy
@inline max_ploidy(v::BgenVariant) = v.genotypes.preamble.max_ploidy
@inline ploidy(v::BgenVariant) = v.genotypes.preamble.ploidy
@inline bit_depth(v::BgenVariant) = v.genotypes.preamble.bit_depth
@inline missings(v::BgenVariant) = v.genotypes.preamble.missings

# The below are valid after calling `minor_allele_dosage!()`
@inline function minor_allele(v::BgenVariant)
    midx = v.genotypes.minor_idx
    if midx == 0
        @error "`minor_allele_dosage!()` must be called before `minor_allele()`"
    else
        v.alleles[v.genotypes.minor_idx]
    end
end

@inline function major_allele(v::BgenVariant)
    midx = v.genotypes.minor_idx
    if midx == 0
        @error "`minor_allele_dosage!()` must be called before `minor_allele()`"
    else
        v.alleles[3 - midx]
    end
end

"""
    destroy_genotypes!(v::BgenVariant)
Destroy any parsed genotype information.
"""
function clear!(v::BgenVariant)
    v.genotypes = nothing
    return
end
