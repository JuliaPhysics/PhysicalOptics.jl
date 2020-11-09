export fft_view_tuple

"""
    fft_view_tuple(x[, shift, dims])


"""
function fft_view_tuple(x, shift=true, dims=nothing; logarithmic=false)
    if isnothing(dims)
        dims = size(x)
    end
    if shift
        x_fft = fftshift(fft(x))
    else
        x_fft = fft(x)
    end

    x_fft_i = imag(x_fft)
    x_fft_r = real(x_fft)
    x_fft_abs = abs.(x_fft)

    if logarithmic
        l(x) = log(1 + abs.(x))
    else
        l = identity
    end

    return l.(x_fft_i), l.(x_fft_r), l.(x_fft_abs)
end
