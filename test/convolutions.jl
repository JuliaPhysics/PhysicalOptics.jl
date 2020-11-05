@testset "Convolution methods" begin
    function conv_test(psf, img, img_out, dims, s)
        otf = fft(psf, dims)
        otf_r = rfft(psf, dims)
        otf_p, conv_p = plan_conv_r(psf, dims)
        @testset "$s" begin
            @test img_out ≈ conv_psf(img, psf, dims)
            @test img_out ≈ conv_otf(img, otf, dims)
            @test img_out ≈ conv_otf_r(img, otf_r, dims)
            @test img_out ≈ conv_p(img, otf_p)
        end
    end
    

    N = 5
    psf = zeros((N, N))
    psf[1, 1] = 1
    img = randn((N, N))
    conv_test(psf, img, img, [1,2], "Convolution random image with delta peak")

    N = 5
    psf = abs.(randn((N, N, 2)))
    img = randn((N, N, 2))
    dims = [1, 2]
    img_out = conv_psf(img, psf, dims)
    conv_test(psf, img, img_out, dims, "Convolution with random 3D PSF and random 3D image over 2D dimensions")

    N = 5
    psf = abs.(randn((N, N, N, N, N)))
    img = randn((N, N, N, N, N))
    dims = [1, 2, 3, 4]
    img_out = conv_psf(img, psf, dims)
    conv_test(psf, img, img_out, dims, "Convolution with random 5D PSF and random 5D image over 4 Dimensions")

    N = 5
    psf = abs.(zeros((N, N, N, N, N)))
    for i = 1:N
        psf[1,1,1,1, i] = 1
    end
    img = randn((N, N, N, N, N))
    dims = [1, 2, 3, 4]
    img_out = conv_psf(img, psf, dims)
    conv_test(psf, img, img, dims, "Convolution with 5D delta peak and random 5D image over 4 Dimensions")



    function conv_test_complex(psf, img, img_out, dims, s)
        otf = fft(psf, dims)
        @testset "$s" begin
            @test img_out ≈ conv_psf(img, psf, dims, real_res=false)
            @test img_out ≈ conv_otf(img, otf, dims, real_res=false)
        end
    end

    N = 5
    psf = abs.(zeros((N, N, N, N, N)))
    for i = 1:N
        psf[1,1,1,1, i] = 1
    end
    img = randn((N, N, N, N, N))
    img = img .+ exp.(1im .* img)
    dims = [1, 2, 3, 4]
    img_out = conv_psf(img, psf, dims, real_res=false)
    conv_test_complex(psf, img, img, dims, "Complex Convolution with 5D delta peak and random 5D image over 4 Dimensions")

    N = 5
    psf = zeros((N, N))
    psf[1, 1] = 1
    img = randn((N, N))
    img = img .+ exp.(1im .* img)
    conv_test_complex(psf, img, img, [1,2], "Complex Convolution random image with delta peak")

    N = 5
    psf = abs.(randn((N, N, 2)))
    img = randn((N, N, 2))
    img = img .+ exp.(1im .* img)
    dims = [1, 2]
    img_out = conv_psf(img, psf, dims, real_res=false)
    conv_test_complex(psf, img, img_out, dims, "Complex Convolution with random 3D PSF and random 3D image over 2D dimensions")

    N = 5
    psf = abs.(randn((N, N, N, N, N)))
    img = randn((N, N, N, N, N))
    img = img .+ exp.(1im .* img)
    dims = [1, 2, 3, 4]
    img_out = conv_psf(img, psf, dims, real_res=false)
    conv_test_complex(psf, img, img_out, dims, "Complex Convolution with random 5D PSF and random 5D image over 4 Dimensions")



end

