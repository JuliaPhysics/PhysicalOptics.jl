export propagate
export fresnel_kernel, rs_kernel

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
