export conv_psf, conv_otf, conv_otf_r, plan_conv_r

"""
    conv_psf(obj, psf [, dims]; mode=true)
Convolve `obj` with `psf` over `dims` dimensions.
Based on FFT convolution.
This function calls `conv_otf`, check the help of this method.

"""
function conv_psf(obj, psf, dims=[1, 2]; real_res=true)
    return conv_otf(obj, fft(psf, dims), dims, real_res=real_res)
end


"""
    conv_otf(obj, otf [ , dims]; real_res=true)
Performs a FFT-based convolution of an `obj`
with an `otf`. `otf = fft(psf)`. The 0 frequency of the `otf` must be located
at position [1, 1, 1].
The `obj` can be of arbitrary dimension but `ndims(obj) ≥ ndims(otf)`.
The convolution happens over the `dims` array. Any further dimensions 
are broadcasted.
Per default `dims = [1, 2]`.
Default `real_res=true` means that the output will be real.
That is useful once `obj` and `psf` are real and therefore mathematically the result 
would be as well.
"""
function conv_otf(obj, otf, dims=[1, 2]; real_res=true)
    res = ifft(fft(obj, dims) .* otf, dims)
    if real_res
        return real(res)
    end
    return res
end


"""
    conv_otf_r(obj, otf [, dims; real_res=true])
Performs a FFT-based convolution of an `obj`
with an `otf`.
Same arguments as `conv_otf` but with `obj` being real and `otf=rfft(psf)`.
All FFTs are computed with `rfft` and `irfft`.
"""
function conv_otf_r(obj, otf, dims=[1, 2]; real_res=true)
    res_fft = rfft(obj, dims) .* otf
    # the output size is not unique 
    # therefore we must specify to irfft which size it should return
    # rfft always shrinkens the array over the first FFT dimension
    # therefore we use dims[1] to find out the output size
    out_size = size(obj)[dims[1]]
    res = irfft(res_fft, out_size, dims)
    if real_res
        return real(res)
    end
    return res
end


"""
    plan_conv_r(psf [, dims]; real_res=true)
Pre-plan an optimized convolution for array shaped like `psf` (based on pre-plan FFT)
along the given dimenions `dims`.
`dims = [1, 2]` per default.
The 0 frequency of `psf` must be located at the first entry.
We return first the `otf` (obtained by `rfft(psf))`.
The second return is the convolution function `conv`.
`conv` itself has two arguments. `conv(obj, otf)` where `obj` is the object and `otf` the otf.
This function achieves faster convolution than `conv_psf(obj, psf)`.
"""
function plan_conv_r(psf, dims=[1, 2]; real_res=true)
    # do the preplanning step
    P = plan_rfft(psf, dims)
    otf = P * psf 
    P_inv = plan_irfft(otf, size(psf)[dims[1]], dims)

    # construct the efficient conv function
    # P and P_inv can be understood like matrices
    # but their computation is fast
    if real_res
        conv(obj, otf) = real.(P_inv * ((P * obj) .* otf))
    else
        conv(obj, otf) = P_inv * ((P * obj) .* otf)
    end
    return otf, conv
end

