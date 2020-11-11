module PhysicalOptics

using FFTW
using SpecialFunctions


 # some useful functions
include("utils.jl")
 # some optics related conversions
include("physical_conversions.jl")
 

 # simple equations to calculate some resolution criterias
include("resolution_equations.jl")


 # fast FFT based convolutions
include("convolutions.jl")

 # utils to get useful objects to view 
include("utils_view.jl")



 # light propagation
include("propagations.jl")


 # point spread functions
include("psf.jl")



end # module
