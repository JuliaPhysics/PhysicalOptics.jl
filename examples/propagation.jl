### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 79d92320-2412-11eb-14a7-b1b5246bbc3c
begin
	using Revise, ImageView, Images, FFTW
	using Plots, PhysicalOptics
end

# ╔═╡ 645e2eac-2433-11eb-0fed-6982080e4951
md" ## Single slit diffraction

The position of the diffraction peaks is:

$$x_k = (2k+1)\frac{z \lambda}{2d}$$

where $z$ is the distance from the slit to the screen, $\lambda$ the wavelength and $d$ the slit width. $\lambda$ is the wavelength.

The analytical solution of the field behind a single slit is in 1D:

$$I \sim \left( \frac{\sin\left(\frac{\pi\,d\,x}{\lambda \, z}\right)}{\left(\frac{\pi\,d\,x}{\lambda \, z}\right)} \right)^2$$

"

# ╔═╡ b0c8db42-2439-11eb-3d9b-b7af3c007f4c
begin
	# define the parameters and define slit
	M = 1024
	slit = zeros(Float64, (1, M))
	d = 100e-6
	L_slit = 10e-3
	z = 0.2
	λ = 550e-9
	xs = LinRange(-L_slit / 2, L_slit / 2, M)
	slit[1, (abs.(xs) .< d / 2)] .= 1
	plot(xs, slit[1, :])
end


# ╔═╡ ff9e9bf6-2445-11eb-0ca7-0f84a22d3c02
begin
	# propagate the initial field and take abs2 for intensity
	slit_prop = abs2.(propagate(slit, L_slit, z, kernel=fresnel_kernel))
	# normalize for a better plot
	slit_prop ./= maximum(slit_prop)
	

	analytical = begin
		x_0 = π .* d .* xs ./ λ ./ z
		(sin.(x_0) ./ x_0).^2
	end
	
	# plots results
	plot(xs, slit_prop[1, :])
	plot!(xs, analytical)
end

# ╔═╡ b0989266-2439-11eb-16da-e75ada56f4c0
collect((fftfreq(100, 1)))

# ╔═╡ 343774e2-2440-11eb-0545-374afacad874


# ╔═╡ 3425c3a0-2440-11eb-08af-570a0aa4a373
md"## Diffraction pattern of a square aperture"

# ╔═╡ 340e8258-2440-11eb-2390-ffb2a54e889c
begin
	N = 256
	rect = zeros(Float64, (N, N))
	L_rect = 0.5
	w_rect = 0.051*2
	x_rect = LinRange(-L_rect/2, L_rect/2, N)
	for (j, y) in enumerate(x_rect)
		for (i, x) in enumerate(x_rect)
			if abs(x) <= w_rect /2 && abs(y) <= w_rect / 2
				rect[i, j] = 1.0
			end
		end
	end
end

# ╔═╡ c798a162-2444-11eb-1a76-ffe0d355671a
begin
	z_rect = 2e3
	# fresnel kernel
	out_rec = propagate(rect, L_rect, z_rect, kernel=fresnel_kernel, λ=500e-9)
	# full Rayleigh Sommerfeld Kernel
	out_recrs = propagate(rect, L_rect, z_rect, kernel=rs_kernel, λ=500e-9)

	imshow([abs2.(out_rec) abs2.(out_recrs) abs2.(rect)])
end

# ╔═╡ e34619f8-2444-11eb-3ce6-41e9921d4b0a
begin
	plot(x_rect, abs2.(out_rec[129, :]))
	plot!(x_rect, abs2.(rect[129, :]))
	plot!(x_rect, abs2.(out_recrs[129, :]))
end

# ╔═╡ a4494ece-2448-11eb-2443-714c62428f9b
fresnel_kernel([1,2], [1,2],1,1)

# ╔═╡ Cell order:
# ╠═79d92320-2412-11eb-14a7-b1b5246bbc3c
# ╟─645e2eac-2433-11eb-0fed-6982080e4951
# ╠═b0c8db42-2439-11eb-3d9b-b7af3c007f4c
# ╠═ff9e9bf6-2445-11eb-0ca7-0f84a22d3c02
# ╠═b0989266-2439-11eb-16da-e75ada56f4c0
# ╠═343774e2-2440-11eb-0545-374afacad874
# ╟─3425c3a0-2440-11eb-08af-570a0aa4a373
# ╠═340e8258-2440-11eb-2390-ffb2a54e889c
# ╠═c798a162-2444-11eb-1a76-ffe0d355671a
# ╠═e34619f8-2444-11eb-3ce6-41e9921d4b0a
# ╠═a4494ece-2448-11eb-2443-714c62428f9b
