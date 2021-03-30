@testset "utils" begin
    @testset "counts" begin
        path = BGEN.datadir("example.8bits.bgen")
        b = Bgen(path)
        for v in iterator(b)
            cnt = counts!(b, v)
            dose = ref_allele_dosage!(b, v)
            correct_cnt = zeros(Int, 512)
            for v in dose
                if !isnan(v)
                    ind = convert(Int, round(v * 255)) + 1
                    correct_cnt[ind] += 1
                end
            end
            correct_cnt[512] = count(isnan.(dose))
            @test all(cnt .== correct_cnt)
        end
    end
end