@testset "select_region" begin

@testset "select_region_null" begin
    chrom, start, stop = "01", 5000, 50000
    b = Bgen(example_16bits)
    @test b.idx !== nothing
    @test length(select_region(b, "02")) == 0
end

@testset "select_whole_chrom" begin
    chrom, start, stop = "01", 5000, 50000
    b = Bgen(example_16bits)
    lt = (x, y) -> isless((x.chrom, x.pos), (y.chrom, y.pos))
    variants = collect(select_region(b, chrom))
    for (x, y) in zip(sort(variants; lt=lt), sort(gen_data; lt=lt))
        @test (x.rsid, x.chrom, x.pos) == (y.rsid, y.chrom, y.pos)
    end
end

@testset "select_after_position" begin
    chrom, start, stop = "01", 5000, 50000
    b = Bgen(example_16bits)
    lt = (x, y) -> isless((x.chrom, x.pos), (y.chrom, y.pos))
    variants = collect(select_region(b, chrom; start=start))
    gen_data_f = filter(x -> x.pos >= start, gen_data)
    for (x, y) in zip(sort(variants; lt=lt), sort(gen_data_f; lt=lt))
        @test (x.rsid, x.chrom, x.pos) == (y.rsid, y.chrom, y.pos)
    end
end

@testset "select_in_region" begin
    chrom, start, stop = "01", 5000, 50000
    b = Bgen(example_16bits)
    lt = (x, y) -> isless((x.chrom, x.pos), (y.chrom, y.pos))
    variants = collect(select_region(b, chrom; start=start, stop=stop))
    gen_data_f = filter(x -> start <= x.pos <= stop, gen_data)
    for (x, y) in zip(sort(variants; lt=lt), sort(gen_data_f; lt=lt))
        @test (x.rsid, x.chrom, x.pos) == (y.rsid, y.chrom, y.pos)
    end
end
end
