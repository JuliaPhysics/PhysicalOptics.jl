


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
end



@testset "rr methods" begin
    out = [1.4142135623730951 1.0 1.4142135623730951; 1.0 0.0 1.0; 1.4142135623730951 1.0 1.4142135623730951]
    out ≈ rr((3,3))
    @test [0] ≈ rr((1, 1))
    @test [2,1,0,1,2] ≈ rr((5, 1))
    @test [3,2,1,0,1,2] ≈ rr((6, 1))

    out = [1.7320508075688772 1.4142135623730951; 1.4142135623730951 1.0; 1.4142135623730951 1.0; 1.0 0.0]
    out = reshape(out, (2,2,2))
    @test out ≈ rr_3D((2,2,2))

end
