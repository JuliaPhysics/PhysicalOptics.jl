export lens 

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
    for (j, x) in enumerate(fftpos(L, size(arr)[2]))
        for (i, y) in enumerate(range(L, size(arr)[1]))
            out[i, j] = circ(radius, x, y) * lens(f, x, y; radius=radius, λ=λ, n=n)
        end
    end

    return out
end
