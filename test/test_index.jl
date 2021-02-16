@testset "index" begin
@testset "index_opens" begin
    @test Bgen(BGEN.datadir("example.15bits.bgen")).idx === nothing
    @test Bgen(BGEN.datadir("example.16bits.bgen")).idx !== nothing
end
@testset "index_region" begin
    chrom = "01"
    start = 5000
    stop = 50000

    b = Bgen(BGEN.datadir("example.16bits.bgen"))
    idx = b.idx
    @test length(select_region(b, chrom)) == length(gen_data)
    @test length(select_region(b, "02")) == 0
    @test length(select_region(b, chrom;
        start=start * 100, stop=stop * 100)) == 0

    after_pos_offsets = select_region(b, chrom; start=start)
    @test length(after_pos_offsets) ==
        length(filter(x -> x.pos >= start, gen_data))

    in_region_offsets = select_region(b, chrom; start=start, stop=stop)
    @test length(in_region_offsets) == length(filter(x -> start <= x.pos <=
        stop, gen_data))
end

end
