export propagate
export fresnel_kernel, rs_kernel
export point_source_propagate
export lens_propagate, four_f_propagate, mla_propagate

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
function rs_kernel(fx, fy, z, λ=λ0, n=1)
    κ = calc_κ(λ, n) 
    c = complex(κ^2 - (fx^2 + fy^2))
    return exp(1im * z * 2 * π * sqrt(c))
end

"""
    propagate!(arr, L, z; kernel=rs_kernel, λ=550e-9, n=1)

Propagate an electric field `arr` with a array size (in physical dimensions like meter etc) of `L` over the distance `z`.
Per default `kernel=rs_kernel` (Rayleigh-Sommerfeld) is the propagation kernel. 
`λ` is the wavelength and `n` the refractive index.
"""
function propagate!(arr, L, z; kernel=rs_kernel, λ=λ0, n=1)
    # trivial propagation
    if iszero(z)
        # return a copy instead of arr
        # therefore we can rely on, that we don't modify the initial array
        return copy(arr)
    end

    # array with the frequencies in Fourier space
    freq_x = to_gpu_or_cpu(arr, Zygote.@ignore eltype(arr).(fftfreq(size(arr)[2], size(arr)[2] / L)'))
    freq_y = to_gpu_or_cpu(arr, Zygote.@ignore eltype(arr).(fftfreq(size(arr)[1], size(arr)[1] / L)))
    arr_ft = fft!(arr)
    
    κ = eltype(arr)(calc_κ(λ, n))
    c_sqrt = Zygote.@ignore (κ^2 .- (freq_x.^2 .+ freq_y.^2))
    c = Zygote.@ignore 1im * eltype(arr)(z * 2π)
    c_exp = Zygote.@ignore (exp.(c .* sqrt.(c_sqrt)))
    out_ft = arr_ft .* c_exp 
    #out_ft = arr_ft .* kernel.(freq_x, freq_y, Ref(z), Ref(λ), Ref(n))
    out = ifft!(out_ft)
    return out
end

function propagate(arr, L, z; kernel=rs_kernel, λ=λ0, n=1)
    propagate!(complex.(arr), L, z, kernel=kernel, λ=λ, n=n)
end



function propagate_tilted_plane(arr, L, z, tilting_ϕ, dim=1; kernel=rs_kernel, λ=λ0, n=1)
    Zygote.@ignore GC.gc(true)
    Δzs = z .+ fftpos(L, size(arr, 1)) .* tan(deg2rad(tilting_ϕ))
    
    arr2 = propagate!(arr, L, Δzs[1])
    if dim == 1
        out = arr2[1:1, :]
    else
        out = arr2[:, 1:1] 
    end
    for (i, Δz) in enumerate(Δzs)
        Δz = Δzs.step
        if isone(i)
            continue
        end
        arr2 = propagate!(arr2, L, Δz, kernel=rs_kernel, λ=λ, n=n)
        if dim == 1
            out = vcat(out, arr2[i:i, :])
        else
            out = hcat(out, arr2[:, i:i])
        end
    end
    return out
end



"""
    point_source_propagate(L, size, point; λ=550e-9, n=1, dtype=ComplexF64)

Propagate a point source in a field of width `L` with array size `size` 
over a distance `z` in a medium with refractive index `n`.
Point source is located at position `point(x, y, z)` and then propagated to 
a location z=0. `z<0` in `point` means forward propagation.

This is based on the analytical solution in real space and not on Fourier
space propagation. The latter one suffers from artifacts while microscopic large `L`.
This function should be always preferred for point sources.
"""
function point_source_propagate(L, size, point::Point; λ=λ0, n=1, dtype=ComplexF64)
    x0, y0, z = point.x, point.y, point.z
    out = zeros(dtype, size)
    k = calc_k(λ, n)
    # if we are at the plane of the point source
    # return numerical delta
    if iszero(z)
        j = argmin(map(x -> abs(x - x0), fftpos(L, size[2])))
        i = argmin(map(y -> abs(y - y0), fftpos(L, size[1])))
        out[i, j] = 1
        return out
    end
    # calculate values
    for (j, x) in enumerate(fftpos(L, size[2]))
        for (i, y) in enumerate(fftpos(L, size[1]))
            # z>0 means that we back propagate the PS to z=0
            r = sign(-z) * sqrt((x-x0)^2 + (y-y0)^2 + z^2)
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
function lens_propagate(arr, L, f; λ=λ0, n=1, d=nothing)
    if isnothing(d)
        d = f
    end
    κ = calc_κ(λ, n)
    # shift to center again
    if typeof(arr) <: Array{<:Real}
        arr = complex.(arr)
    end
    out_f = fftshift(fft!(ifftshift(arr)))
    dx = L / size(arr)[2]
    dy = L / size(arr)[1]

    # lens performs scaled fourier transform. Therefore new field size
    L_new = λ * f / dx
    aux_x =  Zygote.@ignore fftpos(L_new, size(arr)[2])'
    aux_y = Zygote.@ignore fftpos(L_new, size(arr)[1])
    x = Zygote.@ignore to_gpu_or_cpu(arr, aux_x)
    y = Zygote.@ignore to_gpu_or_cpu(arr, aux_y)
    c_exp = 1im * eltype(arr)(π * κ / f * (1 - d / f))
    c = Zygote.@ignore  1 ./ (1im .* λ .* f) .* exp.(c_exp .* (x.^2 .+ y.^2))
    out = out_f .* c .* eltype(out_f)(dx * dy)
    
    return out, L_new
end




"""
    four_f_propagate(arr, L, f1, f2, NA)

Propagate the electrical field `arr` (field size `L`) from front
focal plane to the back focal plane of a 4f system.
The focal length of the first lens is `f1` and the second lens `f`.
Magnification of the system is then `f2/f1`.
"""
function four_f_propagate(arr, L, f1, f2, NA)
	radius = NA * 2 * f1 / 2 
	
    out = lens_propagate(arr, L, f1)
    E1 = out[1]
    L1 = out[2]
    # check that the field size is large enough that the pupil fits
    @assert radius < L1 / 2
	E1_ = circ(E1, L1, radius)

    # tuple of E2 and L2
    out2 = lens_propagate(E1_, L1, f2)
    return out2 
end



"""
    mla_propagate(E, L, fmla, p_mla, z1, z2)

Propagate an electrical field `E` placed a distance `z1` in front
of a microlens array, to a distance `z2` behind the microlens array.
`L` is the size of the field and `fmla` is the focal length of the lenslets.
`p_mla` is microlenses pitch.
"""
function mla_propagate(E, L, fmla, p_mla, z1, z2)

    E_front = propagate(E, L, z1)
	mla = micro_lens_array(fmla, L, p_mla, size(E))
	E_mla = E_front .* mla
	E_behind = propagate(E_mla, L, z2);
	return E_behind, real(mla)
end
