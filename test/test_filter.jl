@testset "filter" begin
    b = Bgen(BGEN.datadir("example.8bits.bgen"))
    vidx = falses(b.header.n_variants)
    vidx[1:10] .= true
    BGEN.filter("test.bgen", b, vidx)
    b2 = Bgen("test.bgen"; sample_path="test.sample")
    @test all(b.samples .== b2.samples)
    for (v1, v2) in zip(iterator(b), iterator(b2)) # length of two iterators are different.
        # it stops when the shorter one (b2) ends.
        @test v1.varid == v2.varid
        @test v1.rsid == v2.rsid
        @test v1.chrom == v2.chrom
        @test v1.pos == v2.pos
        @test v1.n_alleles == v2.n_alleles
        @test all(v1.alleles .== v2.alleles)
        decompressed1 = BGEN.decompress(b.io, v1, b.header)
        decompressed2 = BGEN.decompress(b2.io, v2, b2.header)
        @test all(decompressed1 .== decompressed2)
    end

    sidx = falses(b.header.n_samples)
    sidx[1:10] .= true
    BGEN.filter("test2.bgen", b, trues(b.header.n_variants), sidx)
    b3 = Bgen("test2.bgen"; sample_path="test2.sample")
    for (v1, v3) in zip(iterator(b), iterator(b3))
        @test v1.varid == v3.varid
        @test v1.rsid == v3.rsid
        @test v1.chrom == v3.chrom
        @test v1.pos == v3.pos
        @test v1.n_alleles == v3.n_alleles
        @test all(v1.alleles .== v3.alleles)
        @test isapprox(probabilities!(b, v1)[:, 1:10], probabilities!(b3, v3); nans=true)
    end

    b4 = Bgen(BGEN.datadir("complex.24bits.bgen"))
    BGEN.filter("test3.bgen", b4, trues(b4.header.n_variants), BitVector([false, false, true, true]))
    b5 = Bgen("test3.bgen"; sample_path = "test3.sample")
    for (v4, v5) in zip(iterator(b4), iterator(b5))
        @test v4.varid == v5.varid
        @test v4.rsid == v5.rsid
        @test v4.chrom == v5.chrom
        @test v4.pos == v5.pos
        @test v4.n_alleles == v5.n_alleles
        @test all(v4.alleles .== v5.alleles)
        @test isapprox(probabilities!(b4, v4)[:, 3:4], probabilities!(b5, v5); nans=true)
    end
    rm("test.bgen", force=true)
    rm("test.sample", force=true)
    rm("test2.bgen", force=true)
    rm("test2.sample", force=true)
    rm("test3.bgen", force=true)
    rm("test3.sample", force=true)
end