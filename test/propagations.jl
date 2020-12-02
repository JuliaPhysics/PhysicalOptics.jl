@testset "1D Test Propagations with single slit example" begin
    

    function test_prop_1D(M, d, z, λ, L_slit)
	    slit = zeros(Float64, (1, M))
	    xs = LinRange(-L_slit / 2, L_slit / 2, M)
	    slit[1, (abs.(xs) .< d / 2)] .= 1


        slit_prop = abs2.(propagate(slit, L_slit, z, kernel=fresnel_kernel, λ=λ))
        slit_proprs = abs2.(propagate(slit, L_slit, z, kernel=rs_kernel, λ=λ))
	    # normalize for a better plot
	    slit_prop ./= maximum(slit_prop)
	    slit_proprs ./= maximum(slit_proprs)


	    analytical(x) = begin
	    	x_0 = π .* d .* x ./ λ ./ z
	    	return (sin.(x_0) ./ x_0).^2
	    end
        
        @test   ≈(analytical.(xs)[200:800], slit_prop[200:800], rtol=0.01)
        @test   ≈(analytical.(xs)[200:800], slit_proprs[200:800], rtol=0.01)

        return slit_prop
    end

    x_1 = test_prop_1D(1000, 100e-6, 0.2, 550e-9, 10e-3)
    x_2 = test_prop_1D(1000, 100e-6, 0.2, 450e-9, 10e-3)
    test_prop_1D(1000, 120e-6, 0.2, 450e-9, 10e-3)

    @test ~ ≈(x_1, x_2, rtol=0.01)
end


@testset "Confirm that lens_propagate, four_f_propagate matches analytical solution with 0.5 percent relative error" begin
    function lens_test(L, N, f, radius, shift)
	    E0 = zeros((N, N))
        # place point source at center
	    E0[center_pos(N), center_pos(N)] = 1
	    E1, L1 = lens_propagate(E0, L, f)
	    E1_ = circ(E1, radius, L1)
	    E2, L2 = lens_propagate(E1_, L1, f)
	    psf_ = abs2.(E2)
	    psf = psf_ ./ sum(psf_)

        psf_analytical = jinc_psf((N, N), L, radius, f; shift=shift)
        

        psf_4f, L2_4f = four_f_propagate(E0, L, f, f, radius/f) 
        if ~ shift
            psf_4f = ifftshift(psf_4f)
        end
        @test L2 == L2_4f

        # exclude borders, which sometimes cause trouble
	    psf1 = psf ./ maximum(psf)
	    psf2 = psf_analytical ./ maximum(psf_analytical)
	    psf_4f = normabs2(psf_4f) 
        if shift == false
            psf1 = ifftshift(psf1)
            slice = 1 : round(Int, 0.5 * N)
        else
            slice = round(Int, 0.25 * N) : round(Int, 0.75 * N)
        end
	    @test ≈(psf1[slice, slice], psf2[slice, slice], rtol=0.005)
	    @test ≈(psf_4f[slice, slice], psf2[slice, slice], rtol=0.005)

        

    end

    lens_test(10e-3, 2048, 100e-3, 1e-3, true)
    lens_test(1e-3, 2049, 100e-3, 10e-3, true)
    lens_test(100e-6, 513, 1e-3, 0.2e-3, false)
end
