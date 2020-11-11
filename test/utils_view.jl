
@testset "Test fft_view_tuple"  begin
    N = 100
    x = randn((N, N))

    yi, yr, ya = fft_view_tuple(x)
    @test ya ≈ abs.(yi * 1im + yr)
    
    yi, yr, ya = fft_view_tuple(x, false)
    @test ya ≈ abs.(yi * 1im + yr)
    
    yi, yr, ya = fft_view_tuple(x, false, [1])
    @test ya ≈ abs.(yi * 1im + yr)
    
    yi, yr, ya = fft_view_tuple(x, false, logarithmic=true)
    x_fft = fft(x)
    @test yi ≈ log.(1 .+ abs.(imag(x_fft)))

end
