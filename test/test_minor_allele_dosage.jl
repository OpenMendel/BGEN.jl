@testset "minor allele dosage" begin
@testset "slow" begin
    path = BGEN.datadir("example.16bits.zstd.bgen")
    b = Bgen(path)
    for v in iterator(b)
        probs = probabilities!(b, v)
        dose = minor_allele_dosage!(b, v)

        # dosages for each allele
        a1 = probs[1, :] .* 2 .+ probs[2, :]
        a2 = probs[3, :] .* 2 .+ probs[2, :]

        dose_correct = sum(a1[.!isnan.(a1)]) < sum(a2[.!isnan.(a2)]) ? a1 : a2
        minor_allele_index = sum(a1[.!isnan.(a1)]) < sum(a2[.!isnan.(a2)]) ? 1 : 2
        @test v.genotypes[1].minor_idx[1] == minor_allele_index
        @test all(isapprox.(dose, dose_correct; atol=2e-7, nans=true))
    end
end
@testset "fast" begin
    path = BGEN.datadir("example.8bits.bgen")
    b = Bgen(path)
    for v in iterator(b)
        probs = probabilities!(b, v)
        dose = minor_allele_dosage!(b, v)

        # dosages for each allele
        a1 = probs[1, :] .* 2 .+ probs[2, :]
        a2 = probs[3, :] .* 2 .+ probs[2, :]

        dose_correct = sum(a1[.!isnan.(a1)]) < sum(a2[.!isnan.(a2)]) ? a1 : a2
        minor_allele_index = sum(a1[.!isnan.(a1)]) < sum(a2[.!isnan.(a2)]) ? 1 : 2
        @test v.genotypes[1].minor_idx[1] == minor_allele_index
        @test all(isapprox.(dose, dose_correct; atol=2e-7, nans=true))
    end
end
@testset "v11" begin
    path = BGEN.datadir("example.v11.bgen")
    b = Bgen(path)
    for v in iterator(b)
        probs = probabilities!(b, v)
        dose = minor_allele_dosage!(b, v)

        # dosages for each allele
        a1 = probs[1, :] .* 2 .+ probs[2, :]
        a2 = probs[3, :] .* 2 .+ probs[2, :]

        dose_correct = sum(a1[.!isnan.(a1)]) < sum(a2[.!isnan.(a2)]) ? a1 : a2
        minor_allele_index = sum(a1[.!isnan.(a1)]) < sum(a2[.!isnan.(a2)]) ? 1 : 2
        @test v.genotypes[1].minor_idx[1] == minor_allele_index
        @test all(isapprox.(dose, dose_correct; atol=7e-5, nans=true))
    end
end
end
