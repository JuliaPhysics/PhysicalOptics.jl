using PhysicalOptics
using Test
using FFTW
using Random

 # fix seed for reproducibility
Random.seed!(42)

include("convolutions.jl")
