@testset "Check physical conversion between different quantities" begin

    @test 12 * 20 / 2 / 100 ≈ calc_NA(100, 20, 12)
    @test 1.33 ≈ calc_NA(100e-3, 200e-3, 1.33)

    λ =  500e-9
    n = 1.123
    k = 2 * π / λ * n
    @test calc_k(λ, n) ≈ k
    @test calc_κ(λ, n) ≈ k / 2 / π
    @test calc_λ(k, n) ≈ λ

    n = 1
    λ =  500e-9
    k = 2 * π / λ * n
    @test calc_k(λ) ≈ k
    @test calc_κ(λ) ≈ k / 2 / π
    @test calc_λ(k) ≈ λ

end
