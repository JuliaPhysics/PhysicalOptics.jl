

@testset "Apertures" begin
    @testset "Circular apertures" begin
        x = 0
        y = 0

        @test circ(1, x, y) == 1
        @test circ(0, x, y) == 1
        @test circ(0, sqrt(x^2 + y^2)) == 1


        x = 1
        y = 1
        @test circ(1, x, 0) == 1
        @test circ(1, 0, y) == 1
        @test circ(1, sqrt(x^2 + y^2)) == 0
        @test circ(1, x, y) == 0
    
        x = 10.0
        y = 10.0
        r = 11.0
        @test circ(r, x, 0) ≈ 1
        @test circ(r, 0, y) ≈ 1
        @test circ(r, x, y) ≈ 0
        @test circ(r, sqrt(x^2 + y^2)) ≈ 0

    
        arr = ones((10, 10))
        @test all(circ(arr, 10, 30) .≈ 1)
        
        arr = ones((3,3))
        @test all([0 1 0; 1 1 1; 0 1 0] .≈ circ(arr, 2, 1))
        
        arr = ones((4,4))
        @test all([0 0 1 0; 0 1 1 1; 1 1 1 1; 0 1 1 1] .≈ circ(arr, 2, 1))
        @test all(arr .≈ circ!(arr, 2, 1))
    end

    @testset "Quadratic aperture" begin

        @test quadratic(ones((5, 5)), 5, 2.5, (0.0, 1.25)) == [0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0]
        @test quadratic(ones((5, 5)), 10, 5) == [0.0 0.0 0.0 0.0 0.0; 0.0 1.0 1.0 1.0 0.0; 0.0 1.0 1.0 1.0 0.0; 0.0 1.0 1.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0]
        @test quadratic(ones((4, 4)), 10, 5) == [0.0 0.0 0.0 0.0; 0.0 1.0 1.0 1.0; 0.0 1.0 1.0 1.0; 0.0 1.0 1.0 1.0]
        
        arr = ones(4,4)
        quadratic!(arr, 10, 5) 
        @test arr == [0.0 0.0 0.0 0.0; 0.0 1.0 1.0 1.0; 0.0 1.0 1.0 1.0; 0.0 1.0 1.0 1.0]

        @test quadratic(ones((5, 6)), (5, 6), (4.9, 5)) == [0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 1.0 1.0 1.0 1.0; 0.0 1.0 1.0 1.0 1.0 1.0; 0.0 1.0 1.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0]

    end


end




