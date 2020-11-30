export propagate
export fresnel_kernel, rs_kernel
export point_source_propagate
export lens_propagate


"""
    fresnel_kernel(fx, fy, z, λ, n=1.0)

Calculates the Fresnel propagation kernel in Fourier space.
`fx`, `f_y` are the spatial frequency. `z` is the propagation distance.
`λ` is the wavelength and `n` the refractive index.
"""
function fresnel_kernel(fx, fy, z, λ, n=1)
    κ = calc_κ(λ, n) 
    return @. exp(1im * 2 * π * κ * z) * exp(-1im * π * z / κ * (fx^2 + fy^2))
end

"""
    rs_kernel(fx, fy, z, λ, n=1.0)

Calculates the Rayleigh-Sommerfeld propagation kernel in Fourier space.
`fx`, `f_y` are the spatial frequencies. `z` is the propagation distance.
`λ` is the wavelength and `n` the refractive index.
"""
function rs_kernel(fx, fy, z, λ, n=1.0)
    κ = calc_κ(λ, n) 
    c = @. complex(κ^2 - (fx^2 + fy^2))
    return @. exp(1im * z * 2 * π * sqrt(c))
end

"""
    propagate(arr, L, z; kernel=rs_kernel, λ=550e-9, n=1)

Propagate an electric field `arr` with a array size (in physical dimensions like meter etc) of `L` over the distance `z`.
Per default `kernel=rs_kernel` (Rayleigh-Sommerfeld) is the propagation kernel. 
`λ` is the wavelength and `n` the refractive index.
"""
function propagate(arr, L, z; kernel=rs_kernel, λ=550e-9, n=1)
    # array with the frequencies in Fourier space
    freq_x = fftfreq(size(arr)[2], size(arr)[2] / L)
    freq_y = fftfreq(size(arr)[1], size(arr)[1] / L)

    # go to Fourier space
    arr_ft = fft(arr)

    out_ft = copy(arr_ft)
    # multiply the kernel for each entry with the field in Fourier space
    for (j, fx) in enumerate(freq_x)
        for (i, fy) in enumerate(freq_y)
            out_ft[i, j] = arr_ft[i, j] * kernel(fx, fy, z, λ, n)
        end
    end
    # go back to real space
    out = ifft(out_ft)

    return out
end

"""
    point_source_propagate(L, size, x0, y0, d; λ=550e-9, n=1)

Propagate a point source in a field of width `L` with array size `size` 
over a distance `d` in a medium with refractive index `n`.
Point source is located at position `(x0, y0)` and then propagated over distance
`d`.

This is based on the analytical solution in real space and not on Fourier
space propagation. The latter one suffers from artifacts while microscopic large `L`.
This function should be always preferred for point sources.
"""
function point_source_propagate(L, size, x0, y0, d; λ=550e-9, n=1)
    out = zeros(ComplexF64, size)
    k = calc_k(λ, n)

    # if we are at the plane of the point source
    # return numerical delta
    if iszero(d)
        out[center_pos(size[1]), center_pos(size[2])] = 1
        return out
    end
    # calculate values
    for (j, x) in enumerate(fftpos(L, size[2]))
        for (i, y) in enumerate(fftpos(L, size[1]))
            r = sqrt((x-x0)^2 + (y-y0)^2 + d^2)
            out[i, j] = 1/r .* exp(1im * k * r);
        end
    end
    out ./= out[argmax(abs2.(out))]
    return out
end



"""
    lens_propagate(arr, L, f; λ=550e-9, n=1, d=nothing)

Propagate an electrical field of width `L` at a distance `d` in
front of a lens with focal length `f` in medium with refractive index `n`
to the back focal plane of the lens.
Per default `d=nothing` meaning that we set `d=f`.
Assume same sampling and size in both dimensions.
As return parameter we get the resulting field and the new field size.
Essentially the lens (assuminig infinite large of that lens) performs
a scaled Fourier transform. Therefore we get a new field size.


 # Reference
Based on:
* "Computational Fourier Optics. A MATLAB Tutorial", D. Voelz, (2011).
* Goodman, Joseph W. Introduction to Fourier optics
"""
function lens_propagate(arr, L, f; λ=550e-9, n=1, d=nothing)
    if isnothing(d)
        d = f
    end
    κ = calc_κ(λ, n)
    # shift to center again
    out = fftshift(fft(ifftshift(arr)))
    dx = L / size(arr)[2]
    dy = L / size(arr)[1]

    # lens performs scaled fourier transform. Therefore new field size
    L_new = λ * f / dx

    # complex pre factor of final results
    for (j, x) in enumerate(fftpos(L_new, size(arr)[2]))
        for (i, y) in enumerate(fftpos(L_new, size(arr)[1]))
            c = 1 / (1im * λ * f) * exp(1im * π * κ / f * (1 - d / f) * (x^2 + y^2))
            out[i, j] *= c * dx * dy
        end
    end
    return out, L_new
end
