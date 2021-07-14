export calc_k, calc_λ, calc_κ
export calc_NA
export fresnel_number

"""
    calc_NA(focal_length, diameter[, n])

Calculate the numerical aperture of a system with `focal_length` 
and `diameter` of the lens. Per default `n=1` is the refractive index.

 # Examples
```julia-repl
julia> calc_NA(100e-3, 200e-3)
1.0

julia> calc_NA(100e-3, 200e-3, 1.33)
1.33
```
"""
function calc_NA(focal_length, diameter, n=1)
    return n * diameter / 2 / focal_length
end


"""
    calc_k(λ, n=1)

Calculate the value of the (angular) wave number with vacuum wavelength `λ` in
medium with refractive index `n`.
It holds: \$ k = \\kappa \\cdot 2 \\pi \$
"""
function calc_k(λ, n=1)
    return 2 * oftype(λ, π) * n / λ
end

"""
    calc_κ(λ, n=1)

Calculate the value of the (non angular) wave number with vacuum wavelength `λ` in
medium with refractive index `n`.
It holds: \$ \\kappa = \\frac{k}{2 \\pi} \$
"""
function calc_κ(λ, n=1)
    return calc_k(λ, n) / 2 / oftype(λ, π)
end


"""
    calc_λ(k, n=1)

Calculate the vacuum wavelength `λ` from the angular wave number `k` in 
medium with refractive index `n`.
It holds: \$ \\lambda = \\frac{n \\cdot 2 \\pi}{k} \$
""" 
function calc_λ(k, n=1)
    return 2 * oftype(k, π) * n / k
end


"""
    fresnel_number(a, L, λ=λ0)

Calculate the Fresnel number where `a` is the characteristic size,
`L` the distance from screen to aperture and `λ` the wavelength.
"""
function fresnel_number(a, L, λ=λ0)
    return a^2 / L / λ
end
