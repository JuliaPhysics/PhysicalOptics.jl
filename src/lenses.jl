export lens 
export micro_lens_array


"""
    lens(f, x, y; λ=550e-9, n=1)

Return complex transmission value at position `x` and `y` 
for a lens.
The lens has a focal length `f` and the surrounding medium has 
refractive index `n`. Wavelength is `λ`.
"""
function lens(f::Number, x::Number, y::Number; radius=Inf, λ=550e-9, n=1)
    κ = calc_κ(λ, n)
    r2 = @. x ^ 2 + y ^ 2

    out = @. circ(radius, x, y) * exp(-1im * π * κ / f * (x^2 + y^2))
    return out
end

"""
    lens(f, L, size=(512, 512); λ=550e-9, n=1)

Return complex lens transmission value for a field of `size` 
and a length of `L`
The lens has a focal length `f` and the surrounding medium has 
refractive index `n`. Wavelength is `λ`.

"""
function lens(f::Number, L::Number, size::Tuple=(512, 512); radius=Inf, λ=550e-9, n=1)
    # initialize output array with complex numbers
    out = zeros(ComplexF64, size)

    # iterate over values
    for (j, x) in enumerate(fftpos(L, size[2]))
        for (i, y) in enumerate(fftpos(L, size[1]))
            out[i, j] = circ(radius, x, y) * lens(f, x, y; radius=radius, λ=λ, n=n)
        end
    end

    return out
end


"""
    micro_lens_array(f, L, pitch, size=(512, 512); centered=true, λ=550e-9, n=1)

Return a micro lens array (MLA) with pitch size `pitch` in the field size of `L`.
All lenslets have the same focal length `f`.
As default `centered=true` means that the there is a micro lens centered around
the center of the array.
"""
function micro_lens_array(f, L, pitch, size=(512, 512); centered=true, λ=550e-9, n=1)
    # pitch size in indices
    mla_size = (size .* pitch) .÷  L
    # ÷ doesn't convert to Int tuple
    mla_size = convert(Tuple{Int, Int}, mla_size)

    # if centered, we embed a odd number of micro lenses
    # into a larger arrays
    # at the end, we do a center extract of this larger array
    # otherwise it is troublesome to calculate the offset and
    # cut off micro lenses
    if centered
        N_mla = ceil(Int, L / pitch)
        if mod(N_mla, 2) == 0
            N_mla = N_mla + 1 
        end

        size_out = size
        size = N_mla .* mla_size
    end
    

    out = zeros(ComplexF64, size)
    # put the micro lenses in the out field
    for j = 1:mla_size[2]:size[2]
        for i = 1:mla_size[1]:size[1]
            # min is required in the case when a full micro lens
            # doesn't fit anymore
            mla_i = min(mla_size[1], size[1] - (i - 1))
            mla_j = min(mla_size[2], size[2] - (j - 1))

            out_i = min(size[1], i - 1 + mla_size[1])
            out_j = min(size[2], j - 1 + mla_size[2])
           
            # add small microlens
            out[i:out_i, j:out_j] = lens(f, pitch, mla_size, λ=λ, n=n, radius=pitch/2)[1:mla_i, 1:mla_j]
        end
    end

    
    if centered
        out = center_extract(out, size_out)
    end

    return out
end
