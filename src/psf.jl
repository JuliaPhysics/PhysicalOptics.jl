export generate_simple_psf
export generate_jinc_psf


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
    generate_simple_psf(psf_size, radius; shift=false)

Generation of an approximate 2D PSF.
`psf_size` is the output size of the PSF. The PSF will be centered
around the point [1, 1],
`radius` indicates the pupil diameter in pixel from which the PSF is generated.

# Examples
```julia-repl
julia> generate_psf([5, 5], 2)
5×5 Array{Float64,2}:
 0.36       0.104721    0.0152786    0.0152786    0.104721
 0.104721   0.0304627   0.00444444   0.00444444   0.0304627
 0.0152786  0.00444444  0.000648436  0.000648436  0.00444444
 0.0152786  0.00444444  0.000648436  0.000648436  0.00444444
 0.104721   0.0304627   0.00444444   0.00444444   0.0304627
```
"""
function generate_simple_psf(psf_size, radius; shift=false)
    mask = rr_2D(psf_size) .<= radius
    mask_ft = fft(mask)
    psf = abs2.(mask_ft)
    return shift_and_norm(psf, shift, true) 
end


"""
    generate_jinc_psf(psf_size, L, radius[, f]; λ=550e-9)

Generate the normalized, incoherent 2D jinc PSF of a circular aperture.
`psf_size` is output array shape. `L` is the width of the array
expressed in meter. `radius` is the aperture radius in meter.
Keyword arguments `λ` and `f=100e-3` represent wavelength and focal
length of the lens respectively.

Reference: Mertz, J. (2019). Introduction to Optical Microscopy (2nd ed.).
"""
function generate_jinc_psf(psf_size, L, radius, f=100e-3; λ=550e-9, shift=false)
    # create grid with radial distances
    rr = rr_2D(psf_size, norm=true) .* L ./ 2
    κ = 1 / λ
    Δk⊥ = 2 * κ * radius / f
    psf = jinc.(π .* Δk⊥ .* rr).^2
    psf = ifftshift(psf)
    return shift_and_norm(psf, shift, true) 
end
