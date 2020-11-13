export simple_psf
export jinc_psf


function shift_and_norm(psf, shift, norm)
    if norm
        psf ./= sum(psf)
    end

    if shift
        psf = fftshift(psf)
    end
    return psf
end


"""
    simple_psf(psf_size, radius; shift=false)

Generation of an approximate 2D PSF.
`psf_size` is the output size of the PSF. The PSF will be centered
around the point [1, 1],
`radius` indicates the pupil diameter in pixel from which the PSF is generated.

# Examples
```julia-repl
julia> simple_psf([5, 5], 2)
5×5 Array{Float64,2}:
 0.36       0.104721    0.0152786    0.0152786    0.104721
 0.104721   0.0304627   0.00444444   0.00444444   0.0304627
 0.0152786  0.00444444  0.000648436  0.000648436  0.00444444
 0.0152786  0.00444444  0.000648436  0.000648436  0.00444444
 0.104721   0.0304627   0.00444444   0.00444444   0.0304627
```
"""
function simple_psf(psf_size, radius; shift=false)
    mask = rr(psf_size) .<= radius
    mask_ft = fft(mask)
    psf = abs2.(mask_ft)
    return shift_and_norm(psf, shift, true) 
end


"""
    jinc_psf(psf_size, L, radius[, f]; λ=550e-9, shift=false)

Generate the normalized, incoherent 2D jinc PSF of a circular aperture.
`psf_size` is output array shape. `L` is the width of the array
expressed in meter. `radius` is the aperture radius in meter.
Keyword arguments `λ` and `f=100e-3` represent wavelength and focal
length of the lens respectively.

Reference: Mertz, J. (2019). Introduction to Optical Microscopy (2nd ed.).
"""
function jinc_psf(psf_size, L, radius, f=100e-3; λ=550e-9, shift=false)
    κ = calc_κ(λ) 
    Δk⊥ = 2 * κ * radius / f
    # create real output psf
    
    psf = zeros(Float64, psf_size)
    # calculate each point
    for (j, x) in enumerate(fftpos(L, psf_size[2]))
        for (i, y) in enumerate(fftpos(L, psf_size[1]))
            r = sqrt(x^2 + y^2)
            psf[i, j] = jinc(π * Δk⊥ * r).^2
        end
    end
    # shift center to top left corner
    psf = ifftshift(psf)
    return shift_and_norm(psf, shift, true) 
end
