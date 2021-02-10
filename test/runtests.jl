using BGEN
using Test, Printf

example_10bits = BGEN.datadir("example.10bits.bgen")
example_sample = BGEN.datadir("example.sample")

header_ref = BGEN.Header(example_10bits)

@testset "header" begin
header_test = BGEN.Header(example_10bits)
@test header_test.offset == 0x0000178c # 6028
@test header_test.header_length == 0x00000014 # 20
@test header_test.n_variants == 0x000000c7 # 199
@test header_test.n_samples == 0x000001f4 # 500
@test header_test.compression == 1
@test header_test.layout == 2
@test header_test.has_sample_ids == 1
end

@testset "samples_separate" begin
n_samples = 500
samples_test2 = BGEN.Samples(example_sample, n_samples)
samples_correct = [(@sprintf "sample_%03d" i) for i in 1:n_samples]
@test all(samples_correct .== samples_test2.samples)
samples_test3 = BGEN.Samples(n_samples)
@test all([string(i) for i in 1:n_samples] .== samples_test3.samples)
end

@testset "bgen" begin
bgen = BGEN.Bgen(example_10bits)
@test bgen.fsize == 223646
@test bgen.header == header_ref
n_samples = bgen.header.n_samples
samples_correct = [(@sprintf "sample_%03d" i) for i in 1:n_samples]
@test all(samples_correct .== bgen.samples.samples)
var = bgen.variants[4]
@test length(bgen.variants) == 199
@test var.offset == 0x0000000000002488
@test var.geno_offset == 0x00000000000024b1
@test var.next_var_offset == 0x0000000000002902
@test var.geno_block_size == 0x00000451
@test var.n_samples == 0x000001f4
@test var.varid == "SNPID_5"
@test var.rsid == "RSID_5"
@test var.chrom == "01"
@test var.pos == 0x00001388
@test var.n_alleles == 2
@test all(var.alleles.== ["A", "G"])
end
