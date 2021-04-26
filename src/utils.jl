export center_extract, center_set!, get_indices_around_center, center_pos
export rr, rr_3D
export fftpos
export normabs2

"""
    normabs2(arr)

Takes abs2 and normalizes maximum to 1

 # Examples
```jldoctest
julia> normabs2([10.0, 12.1])
2-element Array{Float64,1}:
 0.6830134553650707
 1.0
```
"""
function normabs2(arr)
    arr2 = abs2.(arr)
    return arr2 ./ maximum(arr2)
end


"""
    fftpos(L, N)

Construct a range from -L/2 to L/2.
However, we ensure that everything is centered around the center
in a way that a FFT interpretes it correctly.
For odd sequences it is indeed in the real center.
For even sequences the center is at `N/2 + 1`.

 # Examples
```jldoctest
julia> collect(fftpos(1, 4))
4-element Array{Float64,1}:
 -0.5
 -0.25
  0.0
  0.25

julia> collect(fftpos(1, 5))
5-element Array{Float64,1}:
 -0.5
 -0.25
  0.0
  0.25
  0.5
```
"""
function fftpos(l, N)
    if iszero(l) || isone(N)
        return zero(l):zero(l)
    end
    if N % 2 == 0
        dx = l / N
        return range(-l/2, l/2-dx, length=N)
    else
        return range(-l/2, l/2, length=N) 
    end
end


function fftpos_old(l, N)
    if iszero(l) && isone(N)
        return l:l
    end
    if N % 2 == 0
        dx = l / N
        return -l/2:dx:l/2-dx
    else
        dx = l / (N-1)
        return -l/2:dx:l/2
    end
end




"""
    get_indices_around_center(i_in, i_out)
A function which provides two output indices `i1` and `i2`
where `i2 - i1 = i_out`
The indices are chosen in a way that the set `i1:i2`
cuts the interval `1:i_in` in a way that the center frequency
stays at the center position.
Works for both odd and even indices
"""
function get_indices_around_center(i_in, i_out)
    if (mod(i_in, 2) == 0 && mod(i_out, 2) == 0 
     || mod(i_in, 2) == 1 && mod(i_out, 2) == 1) 
        x = (i_in - i_out) ÷ 2
        return 1 + x, i_in - x
    elseif mod(i_in, 2) == 1 && mod(i_out, 2) == 0
        x = (i_in - 1 - i_out) ÷ 2
        return 1 + x, i_in - x - 1 
    elseif mod(i_in, 2) == 0 && mod(i_out, 2) == 1
        x = (i_in - (i_out - 1)) ÷ 2
        return 1 + x, i_in - (x - 1)
    end
end

"""
    center_extract(arr, new_size_array)
Extracts a center of an array. 
`new_size_array` must be list of sizes indicating the output
size of each dimension. Centered means that a center frequency
stays at the center position. Works for even and uneven.
If `length(new_size_array) < length(ndims(arr))` the remaining dimensions
are untouched and copied.
# Examples
```jldoctest
julia> center_extract([1 2; 3 4], [1]) 
1×2 Array{Int64,2}:
 3  4
julia> center_extract([1 2; 3 4], [1, 1])
1×1 Array{Int64,2}:
 4
julia> center_extract([1 2 3; 3 4 5; 6 7 8], [2 2])
2×2 Array{Int64,2}:
 1  2
 3  4
```
"""
function center_extract(arr::AbstractArray, new_size_array)
    new_size_array = collect(new_size_array)

    # we construct two lists
    # the reason is, that we don't change higher dimensions which are not 
    # specified in new_size_array
    out_indices1 = [get_indices_around_center(size(arr)[x], new_size_array[x]) 
                    for x = 1:length(new_size_array)]
    
    out_indices1 = [x[1]:x[2] for x = out_indices1]
    
    # out_indices2 contains just ranges covering the full size of each dimension
    out_indices2 = [1:size(arr)[i] for i = (1 + length(new_size_array)):ndims(arr)]
    return arr[out_indices1..., out_indices2...]
end


"""
    center_set!(arr_large, arr_small)
Puts the `arr_small` central into `arr_large`.
The convention, where the center is, is the same as the definition
as for FFT based centered.
Function works both for even and uneven arrays.
# Examples
```jldoctest
julia> center_set!([1, 1, 1, 1, 1, 1], [5, 5, 5])
6-element Array{Int64,1}:
 1
 1
 5
 5
 5
 1
```
"""
function center_set!(arr_large, arr_small)
    out_is = []
    for i = 1:ndims(arr_large)
        a, b = get_indices_around_center(size(arr_large)[i], size(arr_small)[i])
        push!(out_is, a:b)
    end

    #rest = ones(Int, ndims(arr_large) - 3)
    arr_large[out_is...] = arr_small
    
    return arr_large
end


"""
    center_pos(x)
Calculate the position of the center frequency.
Size of the array is `x`
# Examples
```jldoctest
julia> center_pos(3)
2
julia> center_pos(4)
3
```
"""
function center_pos(x::Integer)
    # integer division
    return div(x, 2) + 1
end


function rr_3D(s)
    rarr = zeros((s...))
    for k = 1:s[3]
        for j = 1:s[2]
            for i = 1:s[1]
                rarr[i, j, k] = sqrt( (i-center_pos(s[1]))^2 + (j-center_pos(s[2]))^2 + (k-center_pos(s[3]))^2)
            end
        end
    end
    return rarr
end

"""
    rr(s; L=nothing)

Generate a array with values being the distance to the center pixel.
`s` specifies the output size of the 2D array.
`L` specifies the length. Can be a single number indicating
the same length for both dimensions or a tuple specifying a length for each 
dimension separately.

# Examples
```julia-repl
julia> PhysicalOptics.rr((4,4))
4×4 Matrix{Float64}:
 2.82843  2.23607  2.0  2.23607
 2.23607  1.41421  1.0  1.41421
 2.0      1.0      0.0  1.0
 2.23607  1.41421  1.0  1.41421

julia> PhysicalOptics.rr((4,4), L=1)
4×4 Matrix{Float64}:
 0.707107  0.559017  0.5   0.559017
 0.559017  0.353553  0.25  0.353553
 0.5       0.25      0.0   0.25
 0.559017  0.353553  0.25  0.353553

julia> PhysicalOptics.rr((4,4), L=(5, 1))
4×4 Matrix{Float64}:
 2.54951  1.34629  0.5   1.34629
 2.51247  1.27475  0.25  1.27475
 2.5      1.25     0.0   1.25
 2.51247  1.27475  0.25  1.27475

julia> PhysicalOptics.rr((5,5), L=(1, 1))
5×5 Matrix{Float64}:
 0.707107  0.559017  0.5   0.559017  0.707107
 0.559017  0.353553  0.25  0.353553  0.559017
 0.5       0.25      0.0   0.25      0.5
 0.559017  0.353553  0.25  0.353553  0.559017
 0.707107  0.559017  0.5   0.559017  0.707107
```
"""
function rr_old(s; L=nothing)
    # in case of odd/even array we need to fix the size
    s2 = s[2] %2 == 0 ? s[2] : s[2] -1
    s1 = s[1] %2 == 0 ? s[1] : s[1] -1
    x = fftpos(s2, s[2])'
    y = fftpos(s1, s[1])

    if isnothing(L) == false
        if typeof(L) <: Number
            x *= L ./ 2 ./ abs(x[1])
            y *= L ./ 2 ./  abs(y[1])
        else
            x *= L[1] ./ 2 ./ abs(x[1])
            y *= L[2] ./ 2 ./ abs(y[1])
        end
    end
    

    r = sqrt.(x.^2 .+ y .^ 2)
    return r
end

function rr(s; L=nothing, cuda=false)
    if isnothing(L) == false
        if typeof(L) <: Number
            x = to_gpu_or_cpu(cuda, fftpos(L, s[2])')
            y = to_gpu_or_cpu(cuda, fftpos(L, s[1]))
        else
            x = to_gpu_or_cpu(cuda, fftpos(L[2], s[2])')
            y = to_gpu_or_cpu(cuda, fftpos(L[1], s[1]))
        end
    else
        s2 = s[2] %2 == 0 ? s[2] : s[2] -1
        s1 = s[1] %2 == 0 ? s[1] : s[1] -1
        x = to_gpu_or_cpu(cuda, fftpos(s2, s[2])')
        y = to_gpu_or_cpu(cuda, fftpos(s1, s[1]))
    end
    

    r = sqrt.(x.^2 .+ y .^ 2)
    return r
end



"""
    jinc(x)

Computes the jinc function which is \$\\text{jinc} = \\frac{J_1(x)}{x}\$
where \$J_1\$ being the first Bessel function.

In future, probably switch to: https://github.com/JuliaMath/SpecialFunctions.jl/issues/264

# Examples
```julia-repl
julia> jinc(0.0)
0.5

julia> jinc(0.0im)
0.5 + 0.0im

julia> jinc(0.0001im)
0.5000000006250005 - 3.0616170016954073e-17im

julia> jinc(1.0f0)
0.44005057f0

julia> jinc.([0.0 0.5 1.0])
1×3 Array{Float64,2}:
 0.5  0.484537  0.440051
```
"""
function jinc(x::T) where T
    if iszero(x)
		return convert(T, 0.5)
	else
		return besselj1(x) / x
	end
end



"""
    apply_rot_symmetry(arr, xs, xs_out, ys_out)

Takes a 1D array `arr` with rotational symmetry about the center point.
`arr` should have values evaluates at the positions `xs`.

`xs_out` and `ys_out` are the positions for the output array.
We interpolate `arr` and evaluate it at distances `sqrt(xs_out[j]^2+ys_out[i]^2)`.

"""
function apply_rot_symmetry(arr, xs, xs_out, ys_out)
    itp = Interpolations.scale(
            interpolate(arr, 
                        BSpline(Cubic(Line(Interpolations.OnCell())))),
            xs)

    out = zeros(eltype(arr), (length(ys_out), length(xs_out)))
    
    for j in axes(out, 2)
        for i in axes(out, 1)
            r = sqrt(xs_out[j]^2+ys_out[i]^2)
            out[i, j] = itp(r) 
        end
    end

    return out

end
