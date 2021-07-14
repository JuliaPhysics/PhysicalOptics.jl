@testset "Generate PSF method" begin
    # large aperture is delta peak
    out = zeros((5, 5))
    out[1,1] = 1

    @test out ≈ simple_psf((5, 5), 100)
    out[1,1] = 0
    out[3,3] = 1
    @test out ≈ simple_psf((5, 5), 100, shift=true)

    # pinhole aperture
    out = ones((10, 10))
    out ./= sum(out)
    @test out ≈ simple_psf((10, 10), 0.01)

    # normalized
    @test 1 ≈ sum(simple_psf((100, 100), 10))
end


@testset "Generate jinc PSF method" begin
    # large aperture is delta peak
    out = zeros((5, 5))
    out[1,1] = 1

    # large aperture is roughly a delta
    out[1,1] = 0
    out[3,3] = 1
    @test out ≈ jinc_psf((5, 5), 1, 1, shift=true)

    # pinhole aperture
    out = ones((10, 10))
    out ./= sum(out)
    @test out ≈ jinc_psf((10, 10), 1, 1e-12)

    # normalized
    @test 1 ≈ sum(jinc_psf((100, 100), 10f0, 2f0, 1; λ=1f0, shift=true))
    @test 1 ≈ sum(jinc_psf((100, 100), 10f0, 2f0, 1e-20; λ=1f0, shift=true))

    
    # check consistency with Rayleigh resolution
    d = 1
    f = 1
    r = rayleigh_criterion(f, d)

    psf = jinc_psf((12, 12), 2 * r, d, f)
    @test psf[1, 6] > psf[1, 7] && psf[1, 7] < psf[1, 8] 
    
    psf = jinc_psf((1024, 1024), 2 * r, d, f)
    @test 1 .+ psf[513, 513] ≈ 1
end
