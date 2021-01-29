


function center_test(x1, x2, x3, y1, y2, y3)
    arr1 = randn((x1, x2, x3))
    arr2 = zeros((y1, y2, y3))

    center_set!(arr2, arr1)
    arr3 = center_extract(arr2, (x1, x2, x3))
    @test arr1 ≈ arr3
end

 # test center set and center extract methods
@testset "center methods" begin
    center_test(4, 4, 4, 6,7,4)
    center_test(5, 4, 4, 7, 8, 4)
    center_test(5, 4, 4, 8, 8, 8)
    center_test(6, 4, 4, 7, 8, 8)


    @test 1 == center_pos(1)
    @test 2 == center_pos(2)
    @test 2 == center_pos(3)
    @test 3 == center_pos(4)
    @test 3 == center_pos(5)
    @test 513 == center_pos(1024)

    @test PhysicalOptics.get_indices_around_center((5), (2)) == (2, 3)
    @test PhysicalOptics.get_indices_around_center((5), (3)) == (2, 4)
    @test PhysicalOptics.get_indices_around_center((4), (3)) == (2, 4)
    @test PhysicalOptics.get_indices_around_center((4), (2)) == (2, 3)
end



@testset "rr methods" begin
    out = [1.4142135623730951 1.0 1.4142135623730951; 1.0 0.0 1.0; 1.4142135623730951 1.0 1.4142135623730951]
    @test out ≈ rr((3,3))
    @test [0] ≈ rr((1, 1))
    @test [2,1,0,1,2] ≈ rr((5, 1))
    @test [3,2,1,0,1,2] ≈ rr((6, 1))

    out = [ sqrt(2)  1.0; 1.0      0.0]
    @test out == rr((2,2))
    
    out = [1.4142135623730951 1.118033988749895 1.0 1.118033988749895; 1.118033988749895 0.7071067811865476 0.5 0.7071067811865476; 1.0 0.5 0.0 0.5; 1.118033988749895 0.7071067811865476 0.5 0.7071067811865476]
    @test out == rr((4,4), norm=true)

    @test rr((5, 5), norm = true) == [1.4142135623730951 1.118033988749895 1.0 1.118033988749895 1.4142135623730951; 1.118033988749895 0.7071067811865476 0.5 0.7071067811865476 1.118033988749895; 1.0 0.5 0.0 0.5 1.0; 1.118033988749895 0.7071067811865476 0.5 0.7071067811865476 1.118033988749895; 1.4142135623730951 1.118033988749895 1.0 1.118033988749895 1.4142135623730951]

    @test rr((5, 5), norm = false) == [2.8284271247461903 2.23606797749979 2.0 2.23606797749979 2.8284271247461903; 2.23606797749979 1.4142135623730951 1.0 1.4142135623730951 2.23606797749979; 2.0 1.0 0.0 1.0 2.0; 2.23606797749979 1.4142135623730951 1.0 1.4142135623730951 2.23606797749979; 2.8284271247461903 2.23606797749979 2.0 2.23606797749979 2.8284271247461903]

    out = [1.7320508075688772 1.4142135623730951; 1.4142135623730951 1.0; 1.4142135623730951 1.0; 1.0 0.0]
    out = reshape(out, (2,2,2))
    @test out ≈ rr_3D((2,2,2))

end


@testset "Test fftpos" begin
    @test fftpos(1, 10) == -0.5:0.1:0.4

    @test collect(fftpos(0.01, 3)) ≈ [-0.005, 0.0, 0.005]

    @test fftpos(0.1, 20) ≈ -0.05:0.005:0.045

    @test fftpos(0.1, 21) == -0.05:0.005:0.05

    for N = 2:30
        @test length(fftpos(123, N)) == N
    end
end



@testset "Test normabs2" begin
    x = [1.0, 1.0]
    @test normabs2(x) == x
    x = [1.0 3.0; 4.0 123123123.2]
    @test normabs2(x) ≈ abs2.(x) / maximum(abs2.(x))
end


@testset "Test apply_rot_symmetry" begin
    
    function test_asr(f, xpos, xpos_out, ypos_out, rtol)
        out_ref = [f(sqrt(x^2 + y^2)) for x in xpos_out, y in ypos_out]


        ev_1D = [f(x) for x in xpos]

        out = PhysicalOptics.apply_rot_symmetry(ev_1D, xpos, xpos_out, ypos_out)
        
        @test ≈(out_ref, out, rtol=rtol)
    end

    x = range(0, 2 * √2, length=250)
    x_out = range(0, 2, length=20)
    test_asr(sin, x, x_out, x_out, 1e-9)
    
    x = range(0, 10 * √2, length=1000)
    x_out = range(0, 10, length=20)
    test_asr(SpecialFunctions.jinc, x, x_out, x_out, 1e-7)

end
