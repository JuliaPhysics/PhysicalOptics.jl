export calc_NA
export rayleigh_criterion


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
    rayleigh_criterion(focal_length, diameter; λ=500e-9)

Calculates the resolution of a microscope according to Rayleigh.
`focal_length` is the focal length of the lens.
`λ` is the wavelength and `diameter` is the diameter of the aperture of the lens.
See [Wikipedia](https://en.wikipedia.org/wiki/Angular_resolution#The_Rayleigh_criterion).

 # Examples
```julia-repl
julia> rayleigh_criterion(100, 200)
3.05e-7

julia> rayleigh_criterion(100, 200, λ=100e-9)
6.099999999999999e-8
```
"""
function rayleigh_criterion(focal_length, diameter; λ=500e-9)
    return rayleigh_criterion(calc_NA(focal_length, diameter), λ=λ) 
end

"""
    rayleigh_criterion(NA; λ=500e-9)

Calculates the resolution of a microscope according to Rayleigh.
`λ` is the wavelength and `NA` the numerical aperture of the system.
See [Wikipedia](https://en.wikipedia.org/wiki/Angular_resolution#The_Rayleigh_criterion).

 # Examples
```julia-repl
julia> rayleigh_criterion(1.0)
3.05e-7

julia> rayleigh_criterion(1.0, λ=100e-9)
6.099999999999999e-8
```
"""
function rayleigh_criterion(NA; λ=500e-9)
    return 0.61 * λ / NA
end
