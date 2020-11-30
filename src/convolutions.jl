# this part includes useful convolution methods in PhysicalOptics 


export conv_v_ft, conv
export conv_psf, conv_otf, conv_otf_r, plan_conv_r



"""
    conv_v_ft(u, v_ft[, dims]; real_res=false)
Convolve `u` with `v_ft` over `dims` dimensions.
Based on FFT convolution.

# Arguments
* `u` is an array in real space.
* `v_ft` is the array to be convolved with in Fourier space. 
    Therefore you have to check yourself that `v` was shifted correctly
    in real space.
* Per default `dims=1:ndims(v_ft)` means that we perform the convolution 
    over all dimensions of `v_ft`. 
    If `dims` is an array with integers, we perform convolution 
    only over these dimensions. Eg. `dims=[1,3]` would perform the convolution
    over the first and third dimension. Second dimension is not convolved.
* Per default `real_res=false` means that the output will be real. Otherwise we cut
    off the imaginary part.


# Examples
Convolution with delta peak is an identity operation
```jldoctest
julia> u = [1 2 3 4 5]
1×5 Array{Int64,2}:
 1  2  3  4  5

julia> v = [1 0 0 0 0]
1×5 Array{Int64,2}:
 1  0  0  0  0

julia> v_ft = fft(v)
1×5 Array{Complex{Float64},2}:
 1.0+0.0im  1.0+0.0im  1.0+0.0im  1.0+0.0im  1.0+0.0im

julia> conv_v_ft(u, v_ft)
1×5 Array{Complex{Float64},2}:
 1.0+0.0im  2.0+0.0im  3.0+0.0im  4.0+0.0im  5.0+0.0im

julia> conv_v_ft(u, v_ft, real_res=true)
1×5 Array{Float64,2}:
 1.0  2.0  3.0  4.0  5.0
```
"""
function conv_v_ft(u, v_ft, dims=1:ndims(v_ft); real_res=false)
    res = ifft(fft(u, dims) .* v_ft, dims)
    if real_res 
        return real(res)
    else
        return res
    end
end


 # define custom adjoint for conv_v_ft
 # so far only defined for the derivative regarding the first component
function ChainRulesCore.rrule(::typeof(conv_v_ft), u, v_ft, dims=1:ndims(v_ft); real_res=false)
    Y = conv_v_ft(u, v_ft, dims, real_res=real_res)
    function conv_pullback(barx)
        z = zero(eltype(u))
        return z, conv_v_ft(barx, conj(v_ft), dims, real_res=real_res), z, z, z
    end 
    return Y, conv_pullback
end

"""
    conv(u, v[, dims]; shift=false, real_res=false)
Convolve `u` with `v` over `dims` dimensions.

# Arguments
 
* `u` is an array in real space.
* `v` is the array to be convolved.
* Per default `dims=1:ndims(v)` means that we perform the convolution 
    over all dimensions of `v`. 
    If `dims` is an array with integers, we perform convolution 
    only over these dimensions. Eg. `dims=[1,3]` would perform the convolution
    over the first and third dimension. Second dimension is not convolved.
* Per default `shift=false` therefore we assume that the center point of `v`
    is alreay ifftshifted to the first entry of the array.
* Per default `real_res=false` means that the output will be real. Otherwise we cut
    off the imaginary part.

 # Examples
1D with FFT over all dimensions. We choose `v` to be a delta peak.
Therefore convolution should act as identity.
```jldoctest
julia> u = [1 2 3 4 5]
1×5 Array{Int64,2}:
 1  2  3  4  5

julia> v = [0 0 1 0 0]
1×5 Array{Int64,2}:
 0  0  1  0  0

julia> conv(u, v)
1×5 Array{Complex{Float64},2}:
 4.0+0.0im  5.0+0.0im  1.0+0.0im  2.0+0.0im  3.0+0.0im

julia> conv(u, v, real_res=true) # since v is not ifftshifted with peak at the first entry, we see a wrong result.
1×5 Array{Float64,2}:
 4.0  5.0  1.0  2.0  3.0

julia> conv(u, v, shift=true, real_res=true)
1×5 Array{Float64,2}:
 1.0  2.0  3.0  4.0  5.0

julia> conv(u, ifftshift(v), real_res=true)
1×5 Array{Float64,2}:
 1.0  2.0  3.0  4.0  5.0
```
2D with FFT with different `dims` arguments.
```jldoctest
julia> u = [1 2 3; 4 5 6]
2×3 Array{Int64,2}:
 1  2  3
 4  5  6

julia> v = [1 0 0; 1 0 0]
2×3 Array{Int64,2}:
 1  0  0
 1  0  0

julia> conv(u, v, [2])
2×3 Array{Complex{Float64},2}:
 1.0+0.0im  2.0+0.0im  3.0+0.0im
 4.0+0.0im  5.0+0.0im  6.0+0.0im

julia> conv(u, v, [2], real_res=true)
2×3 Array{Float64,2}:
 1.0  2.0  3.0
 4.0  5.0  6.0

julia> conv(u, v, [1, 2], real_res=true) # now we do a 2D convolution which is not identity anymore
2×3 Array{Float64,2}:
 5.0  7.0  9.0
 5.0  7.0  9.0

julia> conv(u, v, real_res=true) # same statement as above
2×3 Array{Float64,2}:
 5.0  7.0  9.0
 5.0  7.0  9.0
```
"""
function conv(u, v, dims=1:ndims(v); shift=false, real_res=false)
    # this means that we must shift the center pixel of v in real space
    # to the first entry. Simply apply ifftshift on v 
    if shift
        v_ft = fft(ifftshift(v), dims)
    else
        v_ft = fft(v, dims)
    end
    
    # hand over to second function
    return conv_v_ft(u, v_ft, dims, real_res=real_res)
end





"""
    conv_psf(obj, psf [, dims]; shift=false)
Convolve `obj` with `psf` over `dims` dimensions.
Based on FFT convolution. 
This function calls `conv`, check the help of this method.
Wrapper for `conv(obj, psf, dims, shift=shift, real_res=true)`
"""
function conv_psf(obj, psf, dims=1:ndims(psf); shift=false)
    return conv(obj, psf, dims, shift=shift, real_res=true)
end



"""
    conv_otf(obj, otf [ , dims])
Performs a FFT-based convolution of an `obj` with `otf`.
Wrapper for `conv_v_ft(obj, otf, dims, real_res=true)`.
Check the help of `conv_v_ft` for more details and examples.
"""
function conv_otf(obj, otf, dims=1:ndims(otf))
    return conv_v_ft(obj, otf, dims, real_res=true)
end


"""
    conv_otf_r(obj, otf [, dims])
Performs a FFT-based convolution of an `obj`
with an `otf`.
Same arguments as `conv_otf` but with `obj` being real and `otf=rfft(psf)`.
All FFTs are computed with `rfft` and `irfft`.
"""
function conv_otf_r(obj, otf, dims=1:ndims(otf); real_res=true)
    res_fft = rfft(obj, dims) .* otf
    # the output size is not unique 
    # therefore we must specify to irfft which size it should return
    # rfft always shrinkens the array over the first FFT dimension
    # therefore we use dims[1] to find out the output size
    out_size = size(obj)[dims[1]]
    res = irfft(res_fft, out_size, dims)
    return real(res)
end


"""
    plan_conv_r(psf [, dims])
Pre-plan an optimized convolution for array shaped like `psf` (based on pre-plan FFT)
along the given dimenions `dims`.
`dims = 1:ndims(psf)` per default.
The 0 frequency of `psf` must be located at the first entry.
We return first the `otf` (obtained by `rfft(psf))`.
The second return is the convolution function `pconv`.
`pconv` itself has two arguments. `pconv(obj, otf)` where `obj` is the object and `otf` the otf.
This function achieves faster convolution than `conv_psf(obj, psf)`.


 # Examples

```jldoctest
julia> u = [1 2 3 4 5]
1×5 Array{Int64,2}:
 1  2  3  4  5

julia> v = [1 0 0 0 0]
1×5 Array{Int64,2}:
 1  0  0  0  0

julia> otf, pconv = plan_conv_r(v)
(Complex{Float64}[1.0 + 0.0im 1.0 + 0.0im … 1.0 + 0.0im 1.0 + 0.0im], PhysicalOptics.var"#conv#49"{FFTW.rFFTWPlan{Float64,-1,false,2,UnitRange{Int64}},AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,2,UnitRange{Int64}},Float64}}(FFTW real-to-complex plan for 1×5 array of Float64
(rdft2-rank>=2/1
  (rdft2-r2hc-rank0-x5)
  (dft-direct-5 "n1fv_5_avx2_128")), 0.2 * FFTW complex-to-real plan for 1×5 array of Complex{Float64}
(rdft2-rank>=2/1
  (rdft2-hc2r-rank0
    (rdft-rank0-iter-ci/1-x5))
  (dft-direct-5 "n1bv_5_avx2_128"))))

julia> pconv(u, otf)
1×5 Array{Float64,2}:
 1.0  2.0  3.0  4.0  5.0
```
"""
function plan_conv_r(psf, dims=1:ndims(psf))
    # do the preplanning step
    P = plan_rfft(psf, dims)
    otf = P * psf 
    P_inv = plan_irfft(otf, size(psf)[dims[1]], dims)
    
    # construct the efficient conv function
    # P and P_inv can be understood like matrices
    # but their computation is fast
    conv(obj, otf) = real.(P_inv * ((P * obj) .* otf))
    return otf, conv
end

