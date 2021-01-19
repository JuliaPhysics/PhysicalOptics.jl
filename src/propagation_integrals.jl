export debye_focus


"""
    debye_focus(L, N, NA, z; λ=550e-9, n=1)

Evaluates the Debye integral to get the electrical field at position `z` with 
respect to the focus of the lens. 
`z>0` means farther away from the lens than the focus
and vice versa.
Lens has a numerical aperture `NA`. 
`λ` is the wavelength and `n` the refractive index.

Returns a quadratic array with size `(N, N)`.
The side length is `L`, specifically `fftpos(L, N)` are the coordinates.
"""
function debye_focus(L, N, NA, z; λ=550e-9, n=1)
    # half opening angle
    α = asin(NA / n)
    # wavevector
    k = calc_k(λ, n)

    # positions we evaluate the integral at
    # the sqrt(2) because we evaluate the integral over the whole 
    # diagonal of the quadratic output array
    # later we rotate the values to get the full solution over a
    # quadratic
    # we only evaluate at positive values, because of symmetry
    
    # we choose more data points  (we actually would only need sqrt(2)* N/2
    # but since we do an interpolation later, we go for more data points
    xs = range(0, sqrt(2) * L / 2, length=10 * N) 

    # 1D array storing the diagonal integral values
    out = zeros(ComplexF64, length(xs))

    # integrand
    fi(θ, r) = besselj0(k*r*sin(θ))*exp(-1im*k*z*cos(θ))*sin(θ)
    # loop over the xs positions we want to evaluate at
    Threads.@threads for j in eachindex(xs)
        integral, err = quadgk(θ -> fi(θ, xs[j]), zero(eltype(α)), α)
        out[j] = integral
    end
    # some prefactor
    pre_c = 2*π*im / λ
    out .*= pre_c

    # output coordinates where we finally want to have the values
    xs_out = fftpos(L, N)
    ys_out = copy(xs_out) 
    
    out_2D = apply_rot_symmetry(out, xs, xs_out, ys_out)

    # normalize to 1
    return out_2D
end

