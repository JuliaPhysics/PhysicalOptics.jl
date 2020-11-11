using PhysicalOptics
using Test
using FFTW
using Random

 # fix seed for reproducibility
Random.seed!(42)

include("convolutions.jl")


include("utils.jl")

include("psf.jl")

include("physical_conversions.jl")

include("resolution_equations.jl")

include("propagations.jl")

include("utils_view.jl")
