
function test_subeltype()
    @test subeltype([1,2,3]) == Int
    @test subeltype(([1,2,3], [4,5,6])) == Int
    @test subeltype(Array{Float64,2}) == Float64
    @test subeltype(Array{Array{BigFloat,1},2}) == BigFloat
end

function test_dimension()
    @test dimension(Int) == 1
    @test dimension(SVector(1,2)) == 2
    @test dimension(SVector{2,Float64}) == 2
    @test dimension((1,2)) == 2
    @test dimension((1,2,3,4,5)) == 5
    @test dimension(Tuple{Int,Int,Float64,Float64,BigFloat}) == 5
    @test dimension((1,2,3,4,5)) == 5
    @test dimension(CartesianIndex(1,2)) == 2
end

const SystemFloat = typeof(1.0)

function test_prectype()
    @test prectype(1.0) == SystemFloat
    @test prectype(big(1.0)) == BigFloat
    @test prectype(1) == typeof(float(1))
    @test prectype(SVector(1,2)) == typeof(float(1))
    @test prectype(SVector(1,big(2))) == typeof(float(big(2)))
    @test prectype(1.0+2.0im) == SystemFloat
    @test prectype([1.0+2.0im, 3.0]) == SystemFloat
    @test prectype((1.0, 2.0, 3.0)) == SystemFloat
    @test prectype((1.0, big(2.0), 3.0+im)) == BigFloat
end

function test_numtype()
    @test numtype(1.0) == SystemFloat
    @test numtype(big(1.0)) == BigFloat
    @test numtype(1) == Int
    @test numtype(SVector(1,2)) == Int
    @test numtype(SVector(1,big(2))) == BigInt
    @test numtype(1.0+2.0im) == Complex{SystemFloat}
    @test numtype([1.0+2.0im, 3.0]) == Complex{SystemFloat}
    @test numtype((1.0, 2.0, 3.0)) == SystemFloat
    @test numtype(1.0, big(2.0), 3.0+im) == Complex{BigFloat}
    @test numtype((1.0, big(2.0), 3.0+im)) == Complex{BigFloat}
end

test_subeltype()
test_dimension()
test_prectype()
test_numtype()
