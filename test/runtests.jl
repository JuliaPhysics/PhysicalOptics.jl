using PhysicalOptics
using Test
using FFTW
using Random

 # fix seed for reproducibility
Random.seed!(42)

tests = [
    "apertures",
    "convolutions",
    "lenses",
    "physical_conversions",
    "propagations",
    "psf",
    "resolution_equations",
    "utils",
    "utils_view",
]

for t in tests
    @testset "$(t)" begin
        include("$t.jl")
    end
end
