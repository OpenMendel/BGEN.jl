struct GenVariant
    chrom::String
    varid::String
    rsid::String
    pos::UInt
    alleles::Vector{String}
    probs::Matrix{Float64}
    ploidy::Vector{UInt8}
end

"""
    load_gen_data()
load data from "example.gen"
"""
function load_gen_data()
    variants = GenVariant[]
    open(BGEN.datadir("example.gen")) do io
        for l in readlines(io)
            tokens = split(l)
            chrom, varid, rsid, pos, ref, alt = tokens[1:6]
            pos = parse(UInt, pos)
            probs = reshape(parse.(Float64, tokens[7:end]), 3, :)
            nan_col_idx = transpose(sum(probs; dims=1)) .== 0
            for (i, v) in enumerate(nan_col_idx)
                if v
                    probs[:, i] .= NaN
                end
            end
            push!(variants, GenVariant(chrom, varid, rsid, pos, [ref, alt], probs, []))
        end
    end
    variants
end

"""
    parse_vcf_samples(format, samples)
Parse sample data from VCF.
"""
function parse_vcf_samples(format, samples)
    samples = [Dict(zip(split(format, ":"), split(x, ":"))) for x in samples]
    ks = Dict([("GT", r"[/|]"), ("GP", ","), ("HP", ",")])
    samples2 = [Dict{String, Vector{SubString{String}}}() for x in samples]
    for (x, y) in zip(samples, samples2)
        for k in keys(x)
            y[k] = split(x[k], ks[k])
        end
    end
    probs = [occursin("GP", format) ? x["GP"] : x["HP"] for x in samples2]
    probs = map(x -> parse.(Float64, x), probs)
    max_len = maximum(length(x) for x in probs)
    probs_out = Matrix{Float64}(undef, max_len, length(samples))
    fill!(probs_out, NaN)
    for i in 1:length(probs)
        probs_out[1:length(probs[i]), i] .= probs[i]
    end
    ploidy = [length(y["GT"]) for y in samples2]
    probs_out, ploidy
end

"""
    load_vcf_data()
Load data from "complex.vcf" for comparison
"""
function load_vcf_data()
    variants = GenVariant[]
    open(BGEN.datadir("complex.vcf")) do io
        for l in readlines(io)
            if startswith(l, "#")
                continue
            end
            tokens = split(l)
            chrom, pos, varid, ref, alts = tokens[1:5]
            pos = parse(UInt, pos)
            format = tokens[9]
            samples = tokens[10:end]
            varid = split(varid, ",")
            if length(varid) > 1
                rsid, varid = varid
            else
                rsid, varid = varid[1], ""
            end
            probs, ploidy = parse_vcf_samples(format, samples)
            var = GenVariant(string(chrom), varid, string(rsid), pos, vcat(String[ref],
                                                split(alts, ",")), probs, ploidy)
            push!(variants, var)
        end
    end
    variants
end

"""
    load_haps_data()
Load data from "haplotypes.haps" for comparison.
"""
function load_haps_data()
    variants = GenVariant[]
    open(BGEN.datadir("haplotypes.haps")) do io
        for l in readlines(io)
            tokens = split(l)
            chrom, varid, rsid, pos, ref, alt = tokens[1:6]
            pos = parse(UInt, pos)
            probs = tokens[7:end]
            probs = [x == "0" ? [1.0, 0.0] : [0.0, 1.0] for x in probs]
            probs = [probs[pos:pos+1] for pos in 1:2:length(probs)]
            probs = [hcat(x[1], x[2]) for x in probs]
            probs = hcat(probs...)
            probs = reshape(probs, 4, :)
            var = GenVariant(chrom, varid, rsid, pos, [ref, alt], probs, [])
            push!(variants, var)
        end
    end
    variants
end

"""
    epsilon(bit_depth)
Max difference expected for the bit depth
"""
@inline function epsilon(bit_depth)
    return 2 / (2 ^ (bit_depth - 1))
end

"""
    array_equal(truth, parsed, bit_depth)
check if the two arrays are sufficiently equal
"""
@inline function array_equal(truth, parsed, bit_depth)
    eps_abs = 3.2e-5
    all(isapprox.(truth, parsed; atol=max(eps_abs, epsilon(bit_depth)),
        nans=true))
end
