### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 1ae29a38-24fb-11eb-33b0-9d171e0cbcd0
begin
	using Revise, Images, ImageView, PhysicalOptics, FFTW, Plots
end

# ╔═╡ 2e643bb6-24fb-11eb-20e7-cf0175edcdb2
md"## Point Spread Function 4f system
Here we simulate the point spread function of a 4f-system. We use the following qualitative procedure:
* Propagate over distance f to the lens
* Multiply with lens transmission
* Propagate another f
* cut out the pupil
* Propagate another f to the second lens
* Multiply with lens transmission
* Propagate to the image plane 

Luckily there is a numerical procedure to calculate the following steps with one Fourier transform:
* Propagate over distance f to the lens
* Multiply with lens transmission
* Propagate another f

Therefore we do only:
* Fourier Transform
* cut out the pupil
* Fourier Transform
"

# ╔═╡ 2e48f63a-24fb-11eb-2e00-47522330d613
begin
	L = 10e-3
	N = 2048
	f = 100e-3
	radius = 0.5e-3
	E0 = zeros((N, N))
	E0[1025, 1025] = 1
	E1, L1 = lens_propagate(E0, L, f)
	@show L1
	E1_ = circ(E1, radius, L1)
	E2, L2 = lens_propagate(E1_, L1, f)
	psf_ = abs2.(E2)
	psf = psf_ ./ sum(psf_)
end

# ╔═╡ 23584178-259a-11eb-13c5-9db519daef75
md"## Compare with analytical solution"

# ╔═╡ 2ca65ee2-2531-11eb-07c0-275d6a39d8cf
begin
	psf_analytical = jinc_psf((N, N), L, radius, f; shift=true)
end

# ╔═╡ 13283706-2534-11eb-19c9-e74e50186bae
begin
	psf1 = psf ./ maximum(psf)
	psf2 = psf_analytical ./ maximum(psf_analytical)
	imshow(cat(psf1, psf2, psf1 .- psf2, dims=3))	
end

# ╔═╡ 8051b950-25a0-11eb-2e2c-f7986397d24e
begin
	xs = fftpos(L, N)
	x1 = 2 * N ÷5 
	x2 = 3 * N ÷5
	slice = (N ÷2 + 1, x1:x2) 
	plot(xs[x1:x2], psf1[slice...])
	plot!(xs[x1:x2], psf2[slice...])
end

# ╔═╡ Cell order:
# ╠═1ae29a38-24fb-11eb-33b0-9d171e0cbcd0
# ╟─2e643bb6-24fb-11eb-20e7-cf0175edcdb2
# ╠═2e48f63a-24fb-11eb-2e00-47522330d613
# ╟─23584178-259a-11eb-13c5-9db519daef75
# ╠═2ca65ee2-2531-11eb-07c0-275d6a39d8cf
# ╠═13283706-2534-11eb-19c9-e74e50186bae
# ╠═8051b950-25a0-11eb-2e2c-f7986397d24e
