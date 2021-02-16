@testset "example files" begin

@testset "bit depths" begin
    for i in 1:32
        b = Bgen(BGEN.datadir("example.$(i)bits.bgen"); idx_path = nothing)
        v = parse_variants(b)
        for (j, (v1, v2)) in enumerate(zip(v, gen_data))
            @test v1.chrom == v2.chrom
            @test v1.varid == v2.varid
            @test v1.rsid == v2.rsid
            @test v1.pos == v2.pos
            @test all(v1.alleles .== v2.alleles)
            @test array_equal(v2.probs, probabilities!(b, v1), i)
        end
    end
end

@testset "zstd" begin
    b = Bgen(BGEN.datadir("example.16bits.zstd.bgen"); idx_path = nothing)
    v = parse_variants(b)
    for (j, (v1, v2)) in enumerate(zip(v, gen_data))
        @test v1.chrom == v2.chrom
        @test v1.varid == v2.varid
        @test v1.rsid == v2.rsid
        @test v1.pos == v2.pos
        @test all(v1.alleles .== v2.alleles)
        @test array_equal(v2.probs, probabilities!(b, v1), 16)
    end
end

@testset "v11" begin
    b = Bgen(BGEN.datadir("example.v11.bgen"); idx_path = nothing)
    v = parse_variants(b)
    for (j, (v1, v2)) in enumerate(zip(v, gen_data))
        @test v1.chrom == v2.chrom
        @test v1.varid == v2.varid
        @test v1.rsid == v2.rsid
        @test v1.pos == v2.pos
        @test all(v1.alleles .== v2.alleles)
        @test array_equal(v2.probs, probabilities!(b, v1), 16)
    end
end

@testset "haplotypes" begin
    b = Bgen(BGEN.datadir("haplotypes.bgen"); idx_path = nothing)
    v = parse_variants(b)
    for (j, (v1, v2)) in enumerate(zip(v, haps_data))
        @test v1.chrom == v2.chrom
        @test v1.varid == v2.varid
        @test v1.rsid == v2.rsid
        @test v1.pos == v2.pos
        @test all(v1.alleles .== v2.alleles)
        @test array_equal(v2.probs, probabilities!(b, v1), 16)
    end
end

@testset "complex" begin
    b = Bgen(BGEN.datadir("complex.bgen"); idx_path = nothing)
    v = parse_variants(b)
    for (j, (v1, v2)) in enumerate(zip(v, vcf_data))
        @test v1.chrom == v2.chrom
        @test v1.varid == v2.varid
        @test v1.rsid == v2.rsid
        @test v1.pos == v2.pos
        @test all(v1.alleles .== v2.alleles)
        @test array_equal(v2.probs, probabilities!(b, v1), 16)
    end
end

@testset "complex bit depths" begin
    for i in 1:32
        b = Bgen(BGEN.datadir("complex.$(i)bits.bgen"); idx_path = nothing)
        v = parse_variants(b)
        for (j, (v1, v2)) in enumerate(zip(v, vcf_data))
            @test v1.chrom == v2.chrom
            @test v1.varid == v2.varid
            @test v1.rsid == v2.rsid
            @test v1.pos == v2.pos
            @test all(v1.alleles .== v2.alleles)
            @test array_equal(v2.probs, probabilities!(b, v1), i)
        end
    end
end

@testset "null" begin
    @test_throws SystemError Bgen(BGEN.datadir("Hello_World.bgen"))
end
end
