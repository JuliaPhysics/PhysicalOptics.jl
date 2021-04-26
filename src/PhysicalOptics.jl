module PhysicalOptics

using ChainRulesCore
using FFTW
using SpecialFunctions
using QuadGK
using Interpolations
using IndexFunArrays
using FourierTools
using Zygote
using CUDA


to_gpu_or_cpu(arr::AbstractArray, x::AbstractArray) = x
is_cuda(arr::AbstractArray) = false
if CUDA.functional()
    is_cuda(arr::CuArray) = true
    to_gpu_or_cpu(b::Bool, x::AbstractArray) = b ? CuArray(x) : x
    to_gpu_or_cpu(x::AbstractArray) = CuArray(x)
    to_gpu_or_cpu(arr::CuArray, x::AbstractArray) = CuArray(x)
else
    to_gpu_or_cpu(b::Bool, x::AbstractArray) = b ? throw("CUDA.functional() = false") : x
    to_gpu_or_cpu(arr::AbstractArray) = arr 
end

export Point



 # define point struct
struct Point{T}
    x::T
    y::T
    z::T
end



 # some useful functions
include("utils.jl")
 # some optics related conversions
include("physical_conversions.jl")
 
 # simple equations to calculate some resolution criterias
include("resolution_equations.jl")

 # fast FFT based convolutions
include("convolutions.jl")


 # apertures 
include("apertures.jl")

 # utils to get useful objects to view 
include("utils_view.jl")

 # light propagation
include("propagations.jl")

 # # light propagation based on analytical integrals
include("propagation_integrals.jl")

 # some lens functions
include("lenses.jl")

 # point spread functions
include("psf.jl")



end # module
