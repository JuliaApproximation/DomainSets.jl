
using DomainSets: convert_numtype, convert_prectype,
    promote_numtype, promote_prectype

function test_dimension()
    @test DomainSets.euclideandimension(Int) == 1
    @test DomainSets.euclideandimension(Float64) == 1
    @test DomainSets.euclideandimension(ComplexF64) == 1
    @test DomainSets.euclideandimension(SVector{2,Float64}) == 2
    @test DomainSets.euclideandimension(MVector{2,Float64}) == 2
    @test DomainSets.euclideandimension(Tuple{Int,Int}) == 2
    @test DomainSets.euclideandimension(Tuple{Int,Float64}) == 2
    @test_throws MethodError DomainSets.euclideandimension(Vector{Float64})
end

function test_components()
    @test components([1,2,3]) == ()
    d = UnionDomain(Point(1),Point(2))
    @test iscomposite(d)
    @test ncomponents(d) == length(components(d))
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
    @test prectype((1.0,)) == Float64
    @test prectype((1.0, 2.0)) == Float64
    @test prectype((1.0, 2.0, 3.0)) == Float64
    @test prectype((1.0, big(2.0), 3.0+im)) == BigFloat
    @test prectype(NTuple{4,Int}) == Float64
    @test @inferred(prectype(1, 2.0)) == Float64
    @test @inferred(prectype((1, 2.0, 3, 40+im))) == Float64

    @test convert_prectype(2, Float64) == 2
    @test convert_prectype(2, Float64) isa Float64
    @test convert_prectype(1.0+im, BigFloat) == 1+im
    @test convert_prectype(1.0+im, BigFloat) isa Complex{BigFloat}
    @test convert_prectype(SA[1,2], Float64) == SA[1.0,2.0]
    @test convert_prectype(SA[1,2], Float64) isa SVector{2,Float64}
    @test convert_prectype(SA[1,2+im], BigFloat) isa SVector{2,Complex{BigFloat}}
    @test_throws ErrorException convert_prectype("a", BigFloat)

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
    @test numtype((1.0,)) == Float64
    @test numtype((1.0, 2.0)) == Float64
    @test numtype((1.0, 2.0, 3.0)) == Float64
    @test numtype((1.0, 2.0, 3.0, 4.0)) == Float64
    @test numtype(1.0, big(2.0), 3.0+im) == Complex{BigFloat}
    @test numtype((1.0, big(2.0), 3.0+im)) == Complex{BigFloat}
    @test @inferred(numtype(1, 2.0)) == Float64
    @test @inferred(numtype((1, 2.0, 3, 40+im))) == Complex{Float64}

    @test convert_numtype(2, Float64) == 2
    @test convert_numtype(2, Float64) isa Float64
    @test convert_numtype(SA[1,2], Float64) == SA[1,2]
    @test convert_numtype(SA[1,2], Float64) isa SVector{2,Float64}
    @test_throws ErrorException convert_numtype("a", BigFloat)

    @test promote_numtype(2) == 2
    @test promote_numtype(2, 3.0) isa Tuple{Float64,Float64}
    @test promote_numtype(2, 3.0+im, big(4)) isa Tuple{Complex{BigFloat},Complex{BigFloat},Complex{BigFloat}}
end

@testset "common functionality" begin
    @testset "dimension" begin
        test_dimension()
    end
    @testset "elements" begin
        test_components()
    end
    @testset "prectype" begin
        test_prectype()
    end
    @testset "numtype" begin
        test_numtype()
    end
end
