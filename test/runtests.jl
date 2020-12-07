using Bgen
using Test

example_10bits = open(Bgen.datadir("example.10bits.bgen"))
@testset "header" begin
header_test = Bgen.Header(example_10bits)
@test header_test.offset == 0x0000178c # 6028
@test header_test.header_length == 0x00000014 # 20
@test header_test.n_variants == 0x000000c7 # 199
@test header_test.n_samples == 0x000001f4 # 500
@test header_test.compression == 1
@test header_test.layout == 2
@test header_test.has_sample_ids == 1

end
