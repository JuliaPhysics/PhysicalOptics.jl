export gaussian_beam!

function gaussian_beam!(arr, L, z, w₀; λ=λ0, n=1)
    L_new = let 
        if typeof(L) <: Number
            (L, L)
        else
            L
        end
    end
    # array with the frequencies in Fourier space
    x = to_gpu_or_cpu(arr, Zygote.@ignore eltype(arr).(fftpos(L_new[2], size(arr, 2))'))
    y = to_gpu_or_cpu(arr, Zygote.@ignore eltype(arr).(fftpos(L_new[1], size(arr, 1))))

    k = calc_k(λ, n)
    r2 = x .^ 2 .+  y .^ 2 
   
    zᵣ = gaussian_zᵣ(z, w₀, λ=λ, n=n)

    wz = gaussian_w(z, zᵣ, w₀)
    Rz = gaussian_R(z, zᵣ)
    ψz = gouy_phase(z, zᵣ)

    arr_out = arr .* w₀ ./ wz .* exp.(-r2/wz.^2) .* exp.(-1im .* (k .* z .+ k .* r2 ./ (2 .* Rz) .- ψz)) 
    return arr_out
end

# paraxial approx
function gaussian_θ(w₀; λ=λ0, n=1)
    return λ / oftype(w₀, π) / n / w₀
    #return atan(gaussian_w(z, gaussian_zᵣ(z, w₀, λ=λ, n=n), w₀) / z)
end


function gaussian_zᵣ(z, w₀; λ=λ0, n=1)
    return oftype(z, π) * w₀^2 * n / λ
end

function gouy_phase(z, zᵣ)
    return atan(z / zᵣ)
end

function gaussian_R(z, zᵣ)
    if iszero(z)
        return oftype(z, Inf)
    else
        return z * (1 + (zᵣ / z)^2)
    end
end

function gaussian_w(z, zᵣ, w₀)
    return w₀ * √(1 + (z / zᵣ)^2)
end
