@testset "getters" begin
b = BGEN.Bgen(example_8bits)
@test fsize(b) == 128746
@test all(samples(b) .== b.samples)
@test n_samples(b) == 500
@test n_variants(b) == 199
@test compression(b) == "Zlib"

v = first(iterator(b; from_bgen_start=true))
@test n_samples(v::Variant) == 500
@test varid(v::Variant) == "SNPID_2"
@test rsid(v::Variant) == "RSID_2"
@test chrom(v::Variant) == "01"
@test pos(v::Variant) == 2000
@test n_alleles(v::Variant) == 2
@test length(alleles(v::Variant)) == 2
@test all(alleles(v) .== ["A", "G"])

minor_allele_dosage!(b, v)
@test phased(v) == 0
@test min_ploidy(v) == 2
@test max_ploidy(v) == 2
@test all(ploidy(v) .== 2)
@test bit_depth(v) == 8
@test length(missings(v)) == 1
@test missings(v)[1] == 1
@test minor_allele(v) == "A"
@test major_allele(v) == "G"
end
