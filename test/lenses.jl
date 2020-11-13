@testset "Test ideal lens" begin
    
    arr = randn((10, 10))
    @test all(abs.(lens(10.0, 12.0, size(arr), radius=1000)) .≈ 1)

    arr2 = zeros(size(arr))
    arr2[6, 6] = 1
    @test arr2 ≈ abs.(lens(10.0, 12.0, size(arr), radius=0))


end
