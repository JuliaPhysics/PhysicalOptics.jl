export circ, circ!
export quadratic!, quadratic

"""
    circ(radius_aperture, x, y)

Respresenting an aperture.
Being 0 if `x^2 + y^2 > radius_aperture^2` and otherwise 1.

# Examples
```julia-repl
julia> circ(1, 0, 1)
1

julia> circ(1, 0.5, 0.5)
1

julia> circ(1, 1, 1)
0
```

"""
function circ(radius_aperture, x, y)
    if x^2 + y ^2 > radius_aperture^2
        return 0 
    else
        return 1 
    end
end


"""
    circ(radius_aperture, r)

Respresenting an aperture.
Being 0 if `r^2 > radius_aperture^2` and otherwise 1.

# Examples
```julia-repl
julia> circ(1, 0)
1

julia> circ(1, 1)
1

julia> circ(1, 1.01)
0
```
"""
function circ(radius_aperture, r)
    if r^2 > radius_aperture ^2
        return 0
    else
        return 1 
    end
end


"""
    circ(arr, L, radius_aperture)

Apply circular aperture of radius `radius_aperture` to an
array which has field width of `L`.
`circ!` is also available.

```jldoctest
julia> circ(ones((5, 5)), 5, 2.5)
5×5 Array{Float64,2}:
 0.0  0.0  1.0  0.0  0.0
 0.0  1.0  1.0  1.0  0.0
 1.0  1.0  1.0  1.0  1.0
 0.0  1.0  1.0  1.0  0.0
 0.0  0.0  1.0  0.0  0.0
```
"""
function circ(arr::AbstractArray, L, radius_aperture)
    return circ!(copy(arr), L, radius_aperture)
end

function circ!(arr::AbstractArray, L, radius_aperture)
    radius_arr = PhysicalOptics.rr(size(arr), L=L)
    arr .*= (radius_aperture .≥ radius_arr)
    return arr
end


"""
    quadratic(arr, L, diameter, [, diameter=(0,0)])

Apply a quadratic aperture of `diameter` to an 
array which has a field width of `L`.
The aperture is shifted by `Δx, Δy` with respect to the center

There is also an in-place version `quadratic!`

```jldoctest
julia> quadratic(ones((4, 4)), 10, 5)
4×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  1.0  1.0  1.0
 0.0  1.0  1.0  1.0
 0.0  1.0  1.0  1.0

julia> quadratic(ones((5, 5)), 10, 5)
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> quadratic(ones((5, 5)), 5, 2.5, (0.0, 1.25)) 
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  0.0  0.0  0.0
```
"""
function quadratic(arr, L, diameter,  offset::NTuple{2, <: Number}=(0, 0))
    return quadratic!(copy(arr), L, diameter, offset)
end

function quadratic!(arr, L, diameter, offset::NTuple{2, <: Number}=(0,0))
    # need to take abs to compare with diameter
    x = abs.(-offset[2] .+ fftpos(L[2], size(arr, 2)))'
    y = abs.(-offset[1] .+ fftpos(L[1], size(arr, 1)))
    
    arr .*= (diameter[2] / 2 .≥ x)
    arr .*= (diameter[1] / 2 .≥ y)

    return arr
end

function quadratic!(arr, L::Number, diameter::Number, offset::NTuple{2, <:Number}=(0,0))
    return quadratic!(arr, (L, L), (diameter, diameter), offset)
end
