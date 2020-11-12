export circ, circ!

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
    circ(arr, radius_aperture, L)

Apply circular aperture of radius `radius_aperture` to an
array which has field width of `L`.
`circ!` is also available.

```jldoctest
julia> circ(ones((5, 5)), 2.5, 5)
5×5 Array{Float64,2}:
 0.0  0.0  1.0  0.0  0.0
 0.0  1.0  1.0  1.0  0.0
 1.0  1.0  1.0  1.0  1.0
 0.0  1.0  1.0  1.0  0.0
 0.0  0.0  1.0  0.0  0.0
```
"""
function circ(arr::AbstractArray, radius_aperture, L)
    return circ!(copy(arr), radius_aperture, L)
end

function circ!(arr::AbstractArray, radius_aperture, L)
    for (j, x) in enumerate(fftpos(L, size(arr)[2]))
        for (i, y) in enumerate(fftpos(L, size(arr)[1]))
            arr[i, j] *= circ(radius_aperture, x, y)
        end
    end
    return arr
end
