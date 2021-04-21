module PhysicalOptics

using ChainRulesCore
using FFTW
using SpecialFunctions
using QuadGK
using Interpolations
using IndexFunArrays


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
