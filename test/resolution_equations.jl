@testset "Test equations for Rayleigh resolution" begin
    f = 100e-3
    D = 50e-3

    @test rayleigh_criterion(calc_NA(f, D)) ≈ rayleigh_criterion(f, D)
    @test rayleigh_criterion(calc_NA(f, D), λ=500e-9) ≈ rayleigh_criterion(f, D)
    @test ~(rayleigh_criterion(calc_NA(f, D)) ≈ rayleigh_criterion(f, D, λ=200e-9))

end


@testset "Calc NA" begin
    @test calc_NA(100, 100) ≈ 0.5
    @test calc_NA(10e-9, 40e-9, 2) ≈ 4
end
