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
