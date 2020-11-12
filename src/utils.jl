export center_extract, center_set!, get_indices_around_center, center_pos
export rr, rr_3D
export jinc
export fftpos


"""
    fftpos(l, N)

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
    if N % 2 == 0
        dx = l / N
        return range(-l/2, l/2-dx, length=N)
    else
        return range(-l/2, l/2, length=N) 
    end
end

"""
    get_indices_around_center(i_in, i_out)

A function which provides two output indices i1 and i2
where i2 - i1 = i_out
The indices are choosen in a way that the set i1:i2
cuts the intervall 1:i_in in a way that the center frequency
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
        return 2 + x, i_in - x 
    elseif mod(i_in, 2) == 0 && mod(i_out, 2) == 1
        x = (i_in - (i_out - 1)) ÷ 2
        return 1 + x, i_in - (x - 1)
    end
end


"""
    center_extract(arr, new_size)

Extracts a center of an array. 
`new_size` must be list of sizes indicating the output
size of each dimension. Centered means that a center frequency
stays at the center position. Works for even and uneven.
If `length(new_size) < length(size(arr))` the remaining dimensions
are untouched and copied.

# Examples
```julia-repl
julia> center_extract([[1,2] [3, 4]], [1])
1×2 Array{Int64,2}:
 2  4

julia> center_extract([[1,2] [3, 4]], [1, 1])
1×1 Array{Int64,2}:
4
```
"""
function center_extract(arr, index_arrays)
    index_arrays = collect(index_arrays)
    out_indices1 = [get_indices_around_center(size(arr)[x], index_arrays[x]) 
                    for x = 1:length(index_arrays)]
    
    out_indices1 = [x[1]:x[2] for x = out_indices1]


    out_indices2 = map(eval, [1:size(arr)[length(out_indices1) + i] for i = (1 + size(index_arrays)[1]):ndims(arr)])
    return view(arr, out_indices1..., out_indices2...)
end


"""
    center_set!(arr_large, arr_small)

Puts the `arr_small` central into `arr_large`.
The convention, where the center is, is the same as the definition
as for FFT based centered.
Function works both for even and uneven arrays.

# Examples
```julia-repl
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
```julia-repl
julia> center_pos(3)
2
julia> center_pos(4)
4
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
    rr(s; norm=false)

Generate a image with values being the distance to the center pixel.
`s` specifies the output size of the 2D array.
`norm` normalizes the values to the total size.

# Examples
```julia-repl
julia> rr((6, 6))
6×6 Array{Float64,2}:
 4.24264  3.60555  3.16228  3.0  3.16228  3.60555
 3.60555  2.82843  2.23607  2.0  2.23607  2.82843
 3.16228  2.23607  1.41421  1.0  1.41421  2.23607
 3.0      2.0      1.0      0.0  1.0      2.0
 3.16228  2.23607  1.41421  1.0  1.41421  2.23607
 3.60555  2.82843  2.23607  2.0  2.23607  2.82843

julia> rr((4,4), norm=true)
4×4 Array{Float64,2}:
 1.41421  1.11803   1.0  1.11803
 1.11803  0.707107  0.5  0.707107
 1.0      0.5       0.0  0.5
 1.11803  0.707107  0.5  0.707107
```
"""
function rr(s; norm=false)
    rarr = zeros((s...)) 
    for j = 1:s[2]
        for i = 1:s[1]
                if norm
                    rarr[i, j] = sqrt( 1/(s[1] ÷ 2)^2 * (i - center_pos(s[1]))^2 + 1 / (s[2] ÷ 2)^2 * (j - center_pos(s[2]))^2)
                else    
                    rarr[i, j] = sqrt( (i - center_pos(s[1]))^2 + (j - center_pos(s[2]))^2)
                end
        end
    end
    return rarr
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
