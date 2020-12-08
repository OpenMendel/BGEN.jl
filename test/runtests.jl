using BGEN
using Test, Printf

example_10bits = open(BGEN.datadir("example.10bits.bgen"))

@testset "header" begin
header_test = BGEN.Header(example_10bits)
@test header_test.offset == 0x0000178c # 6028
@test header_test.header_length == 0x00000014 # 20
@test header_test.n_variants == 0x000000c7 # 199
@test header_test.n_samples == 0x000001f4 # 500
@test header_test.compression == 1
@test header_test.layout == 2
@test header_test.has_sample_ids == 1

n_samples = header_test.n_samples
samples_test = BGEN.Samples(example_10bits, n_samples)
samples_correct = [(@sprintf "sample_%03d" i) for i in 1:n_samples]
@test all(samples_correct .== samples_test.samples)
samples_test2 = BGEN.Samples(BGEN.datadir("example.sample"), n_samples)
@test all(samples_correct .== samples_test2.samples)
samples_test3 = BGEN.Samples(n_samples)
@test all([string(i) for i in 1:n_samples] .== samples_test3.samples)
end
