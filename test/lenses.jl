@testset "Test ideal lens" begin
    
    arr = randn((10, 10))
    @test all(abs.(lens(10.0, 12.0, size(arr), radius=1000)) .≈ 1)

    arr2 = zeros(size(arr))
    arr2[6, 6] = 1
    @test arr2 ≈ abs.(lens(10.0, 12.0, size(arr), radius=0))
end


@testset "Test microlens" begin
    mla = micro_lens_array(1000e-3, 10e-3, 10e-3, (64, 64))
    l = lens(1000e-3, 10e-3, (64, 64), radius=10e-3/2)
    @test mla ≈ l
    
    mla = micro_lens_array(1000e-3, 10e-3, 10e-3, (64, 64), aperture=false)
    l = lens(1000e-3, 10e-3, (64, 64), radius=Inf)
    @test mla ≈ l

    
    mla = micro_lens_array(1000e-3, 20e-3, 20e-3 / 2, (128, 128), aperture=false, centered=false)
    l = lens(1000e-3, 10e-3, (64, 64), radius=Inf)
    @test mla[1:64, 1:64] ≈ l
    
    mla = micro_lens_array(1000e-3, 20e-3, 20e-3 / 2, (128, 128), aperture=false, centered=true)
    l = lens(1000e-3, 10e-3, (64, 64), radius=Inf)
    @test center_extract(mla, (64, 64)) ≈ l


end
