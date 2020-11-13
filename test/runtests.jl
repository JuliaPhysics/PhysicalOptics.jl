using PhysicalOptics
using Test
using FFTW
using Random

 # fix seed for reproducibility
Random.seed!(42)

tests = [
    "convolutions",
    "utils",
    "psf",
    "apertures",
    "physical_conversions",
    "resolution_equations",
    "propagations",
    "utils_view"
]

for t in tests
    @testset "$(t)" begin
        include("$t.jl")
    end
end
