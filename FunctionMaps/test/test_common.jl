function test_dimension()
    @test euclideandimension(Int) == 1
    @test euclideandimension(Float64) == 1
    @test euclideandimension(ComplexF64) == 1
    @test euclideandimension(SVector{2,Float64}) == 2
    @test euclideandimension(MVector{2,Float64}) == 2
    @test euclideandimension(Tuple{Int,Int}) == 2
    @test euclideandimension(Tuple{Int,Float64}) == 2
    @test_throws ArgumentError euclideandimension(Vector{Float64})
end

function test_prectype()
    @test prectype(1.0) == Float64
    @test prectype(big(1.0)) == BigFloat
    @test prectype(1) == typeof(float(1))
    @test prectype(SVector(1,2)) == typeof(float(1))
    @test prectype(SVector(1,big(2))) == typeof(float(big(2)))
    @test prectype(1.0+2.0im) == Float64
    @test prectype([1.0+2.0im, 3.0]) == Float64
    @test prectype(NTuple{2,Int}) == Float64
    @test prectype(typeof((1.0,))) == Float64
    @test prectype(typeof((1.0, 2.0))) == Float64
    @test prectype(typeof((1.0, 2.0, 3.0))) == Float64
    @test prectype(typeof((1.0, big(2.0), 3.0+im))) == BigFloat
    @test prectype(NTuple{4,Int}) == Float64
    @test @inferred(prectype(1, 2.0)) == Float64
    @test @inferred(prectype(typeof((1, 2.0, 3, 40+im)))) == Float64

    @test convert_prectype(Float64, 2) == 2
    @test convert_prectype(Float64, 2) isa Float64
    @test convert_prectype(BigFloat, 1.0+im) == 1+im
    @test convert_prectype(BigFloat, 1.0+im) isa Complex{BigFloat}
    @test convert_prectype(Float64, SA[1,2]) == SA[1.0,2.0]
    @test convert_prectype(Float64, SA[1,2]) isa SVector{2,Float64}
    @test convert_prectype(Float64, MVector(1,2)) == MVector(1.0,2.0)
    @test convert_prectype(Float64, MVector(1,2)) isa MVector
    @test convert_prectype(Float64, SA[1 2; 3 4]) == SA[1.0 2.0; 3.0 4.0]
    @test convert_prectype(Float64, SA[1 2; 3 4]) isa SMatrix{2,2,Float64}
    @test convert_prectype(Float64, MMatrix{2,2}(1, 2, 3, 4)) == SA[1.0 3.0; 2.0 4.0]
    @test convert_prectype(Float64, MMatrix{2,2}(1, 2, 3, 4)) isa MMatrix{2,2,Float64}
    @test convert_prectype(BigFloat, SA[1,2+im]) isa SVector{2,Complex{BigFloat}}
    @test_throws ArgumentError convert_prectype(BigFloat, "a")

    @test promote_prectype(2) == 2
    @test promote_prectype(2, 3.0) isa Tuple{Float64,Float64}
    @test promote_prectype(2, 3.0+im, big(4)) isa Tuple{BigFloat,Complex{BigFloat},BigFloat}
end


function test_numtype()
    @test numtype(1.0) == Float64
    @test numtype(big(1.0)) == BigFloat
    @test numtype(1) == Int
    @test numtype([1,2,3], [4,5,6]) == Int
    @test numtype(SVector(1,2)) == Int
    @test numtype(SVector(1,big(2))) == BigInt
    @test numtype(Array{Float64,2}) == Float64
    @test numtype(1.0+2.0im) == Complex{Float64}
    @test numtype([1.0+2.0im, 3.0]) == Complex{Float64}
    @test numtype(NTuple{2,Complex{Int}}) == Complex{Int}
    @test numtype(typeof((1.0,))) == Float64
    @test numtype(typeof((1.0, 2.0))) == Float64
    @test numtype((1.0, 2.0, 3.0)) == Float64
    @test numtype((1.0, 2.0, 3.0, 4.0)) == Float64
    @test numtype(1.0, big(2.0), 3.0+im) == Complex{BigFloat}
    @test numtype(typeof((1.0, big(2.0), 3.0+im))) == Complex{BigFloat}
    @test @inferred(numtype(1, 2.0)) == Float64
    @test @inferred(numtype(typeof((1, 2.0, 3, 40+im)))) == Complex{Float64}

    @test convert_numtype(Float64, 2) == 2
    @test convert_numtype(Float64, 2) isa Float64
    @test convert_numtype(Float64, SA[1,2]) == SA[1.0,2.0]
    @test convert_numtype(Float64, SA[1,2]) isa SVector{2,Float64}
    @test convert_numtype(Float64, MVector(1,2)) == MVector(1.0,2.0)
    @test convert_numtype(Float64, MVector(1,2)) isa MVector
    @test convert_numtype(Float64, SA[1 2; 3 4]) == SA[1.0 2.0; 3.0 4.0]
    @test convert_numtype(Float64, SA[1 2; 3 4]) isa SMatrix{2,2,Float64}
    @test convert_numtype(Float64, MMatrix{2,2}(1, 2, 3, 4)) == SA[1.0 3.0; 2.0 4.0]
    @test convert_numtype(Float64, MMatrix{2,2}(1, 2, 3, 4)) isa MMatrix{2,2,Float64}
    @test_throws ArgumentError convert_numtype(BigFloat, "a")

    @test promote_numtype(2) == 2
    @test promote_numtype(2, 3.0) isa Tuple{Float64,Float64}
    @test promote_numtype(2, 3.0+im, big(4)) isa Tuple{Complex{BigFloat},Complex{BigFloat},Complex{BigFloat}}
end

@testset "common functionality" begin
    @testset "dimension" begin
        test_dimension()
    end
    @testset "prectype" begin
        test_prectype()
    end
    @testset "numtype" begin
        test_numtype()
    end
end
