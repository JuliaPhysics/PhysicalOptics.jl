module PhysicalOptics

using FFTW
using SpecialFunctions

 # convolutions
include("convolutions.jl")


include("utils.jl")

 # utils to view conveniently images
include("utils_view.jl")


include("psf.jl")


 # simple equations to calculate some resolutions
include("resolution_equations.jl")


end # module
