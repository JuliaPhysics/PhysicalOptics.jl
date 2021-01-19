using PhysicalOptics
using Test
using FFTW
using Random
using SpecialFunctions


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


@testset "Point test" begin
    p = Point(1.0, 2.0, 3.0)

    @test typeof(p) === Point{Float64}
    @test p.x == 1.0
    @test p.y == 2.0
    @test p.z == 3.0
    

    p = Point(1, 2, 3)
    @test typeof(p) === Point{Int64}
    @test p.x == 1
    @test p.y == 2
    @test p.z == 3

end


for t in tests
    @testset "$(t)" begin
        include("$t.jl")
    end
end
