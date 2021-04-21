
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


@testset "Confirm that lens_propagate, four_f_propagate, Debye integral match analytical solution and are consistent to each other" begin
    function lens_test(L, N, f, radius, shift)

        NA = radius/f
        
	    E0 = zeros((N, N))
        # place point source at center
	    E0[center_pos(N), center_pos(N)] = 1
	    E1, L1 = lens_propagate(E0, L, f)
	    E1_ = circ(E1, L1, radius)
	    E2, L2 = lens_propagate(E1_, L1, f)
	    psf_ = abs2.(E2)
	    psf = psf_ ./ sum(psf_)
	    psf1 = psf ./ maximum(psf)

        psf_analytical = jinc_psf((N, N), L, radius, f; shift=shift)
	    psf2 = psf_analytical ./ maximum(psf_analytical)
        
        
        E_deb = debye_focus(L, N, NA, 0.0)
        psf_deb = abs2.(E_deb)
        psf_deb ./= sum(psf_deb)
        psf_deb ./= maximum(psf_deb)


        psf_4f, L2_4f = four_f_propagate(E0, L, f, f, NA) 
	    psf_4f = normabs2(psf_4f) 
        if ~ shift
            psf_deb = ifftshift(psf_deb)
            psf1 = ifftshift(psf1)
            psf_4f = ifftshift(psf_4f)
            # exlude borders of the fields because they can differ
            slice = 1 : round(Int, 0.5 * N)
        else
            slice = round(Int, 0.25 * N) : round(Int, 0.75 * N)
        end


        @test L2 == L2_4f
	    @test ≈(psf1[slice, slice], psf2[slice, slice], rtol=0.005)
	    @test ≈(psf_4f[slice, slice], psf2[slice, slice], rtol=0.005)
        @test ≈(psf_deb[slice, slice], psf_4f[slice, slice], rtol=0.01)
        
    end

    lens_test(10e-3, 2048, 100e-3, 1e-3, true)
    lens_test(1e-3, 2049, 100e-3, 10e-3, true)
    lens_test(100e-6, 513, 1e-3, 0.2e-3, false)
end


@testset "Point source propagation tests" begin
    N = 100
    L = 100e-6
    ps = point_source_propagate(L, (N, N), Point(0.0, 0.0, 0.0))
    ref = zeros(ComplexF64, (N, N))
    ref[center_pos(N), center_pos(N)] = 1

    @test ps == ref

    ps = point_source_propagate(L, (N, N), Point(10e-6, 0.0, 10.0e-6))
    ref = zeros(ComplexF64, (N, N))

    k = 2π/550e-9
    for (j, x) in enumerate(fftpos(L, N))
        for (i, y) in enumerate(fftpos(L, N))
            r = sqrt((x-10e-6)^2 + (y)^2 + (10.0e-6)^2)
            ref[i, j] = 1/r .* exp(-1im * k * r);
        end
    end
    ref ./= ref[argmax(abs2.(ref))]
    @test ps ≈ ref


    ps = point_source_propagate(L, (N, N), Point(-50.0e-6, -50.0e-6, 0.0))
    @test ps[1, 1] == 1 

end

@testset "Compare Debye with 4F for different NAs" begin

    function compare(N, L, NA)
	    out_s = (N, N)
	    #f = 50-3
	    out = debye_focus(L, N, NA, 0e-6)
	    h1 = abs2.(out)
	    h1 ./= sum(h1)
        
        E0 = zeros(ComplexF64, (N, N))
	    E0[center_pos(N), center_pos(N)] = 1
	    E1, L2 = four_f_propagate(E0, L, 100e-3, 100e-3, NA)
	    h2 = abs2.(E1)
	    h2 ./= sum(h2)
        
        slice = round(Int, 0.3*N):round(Int, 0.7 * N)
        @test ≈(h1[slice, slice], h2[slice, slice], rtol=0.05) 
    end

    compare(256, 20e-6, 0.2)
    compare(256, 50e-6, 0.1)
    compare(256, 10e-6, 0.7)
    compare(512, 10e-6, 0.4)
end




