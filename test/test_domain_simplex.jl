function test_simplex()
    @test UnitSimplex(2) isa VectorUnitSimplex{Float64}
    @test UnitSimplex(Val(2)) isa EuclideanUnitSimplex{2,Float64}
    @test UnitSimplex{Float64}() isa StaticUnitSimplex{Float64}
    @test UnitSimplex{Float64}(1) isa StaticUnitSimplex{Float64}
    @test_throws AssertionError UnitSimplex{Float64}(2)
    @test UnitSimplex{SVector{2,Float64}}(Val(2)) isa EuclideanUnitSimplex{2,Float64}
    @test_throws AssertionError UnitSimplex{SVector{2,Float64}}(Val(3))
    @test UnitSimplex{SVector{2,Float64}}(2) isa EuclideanUnitSimplex{2,Float64}
    @test_throws AssertionError UnitSimplex{SVector{2,Float64}}(3)
    @test UnitSimplex{SVector{2,Float64}}() isa EuclideanUnitSimplex{2,Float64}
    @test UnitSimplex{Vector{Float64}}(2) isa VectorUnitSimplex{Float64}
    @test_throws MethodError UnitSimplex{Vector{Float64}}()

    @test UnitSimplex{Float64,:open}() isa StaticUnitSimplex{Float64,:open}
    @test UnitSimplex{Float64,:closed}(1) isa StaticUnitSimplex{Float64,:closed}
    @test_throws AssertionError UnitSimplex{Float64,:closed}(2)
    @test UnitSimplex{SVector{2,Float64},:open}(Val(2)) isa EuclideanUnitSimplex{2,Float64,:open}
    @test_throws AssertionError UnitSimplex{SVector{2,Float64},:closed}(Val(3))
    @test UnitSimplex{SVector{2,Float64},:open}(2) isa EuclideanUnitSimplex{2,Float64,:open}
    @test_throws AssertionError UnitSimplex{SVector{2,Float64},:closed}(3)
    @test UnitSimplex{SVector{2,Float64},:open}() isa EuclideanUnitSimplex{2,Float64,:open}
    @test UnitSimplex{Vector{Float64},:closed}(2) isa VectorUnitSimplex{Float64,:closed}
    @test_throws MethodError UnitSimplex{Vector{Float64},:open}()

    @test StaticUnitSimplex(Val(3)) isa StaticUnitSimplex{SVector{3,Float64}}
    @test DynamicUnitSimplex{Float64}(1) isa DynamicUnitSimplex{Float64}
    @test_throws AssertionError DynamicUnitSimplex{Float64}(2)

    @test repr(UnitSimplex(Val(2))) == "UnitSimplex(Val(2))"
    @test repr(UnitSimplex(3)) == "UnitSimplex(3)"

    d = UnitSimplex(Val(2))
    # We test a point in the interior, a point on each of the boundaries and
    # all corners.
    @test SA[0.2,0.2] ∈ d
    @test SA[0.0,0.2] ∈ d
    @test SA[0.2,0.0] ∈ d
    @test SA[0.5,0.5] ∈ d
    @test SA[0.0,0.0] ∈ d
    @test SA[1.0,0.0] ∈ d
    @test SA[0.0,1.0] ∈ d
    # And then some points outside
    @test SA[0.6,0.5] ∉ d
    @test SA[0.5,0.6] ∉ d
    @test SA[-0.2,0.2] ∉ d
    @test SA[0.2,-0.2] ∉ d
    @test boundingbox(d) == UnitCube{SVector{2,Float64}}()

    # issue #102
    @test !([0.3,0.4,0.2] ∈ UnitSimplex(2))
    z1 = @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain 'UnitSimplex(Val(2))' with dimension 2. Returning false.") !([0.3,0.4,0.2] ∈ UnitSimplex(Val(2)))
    @test z1
    z2 = @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain 'UnitSimplex(Val(2))' with dimension 2. Returning false.") !(SVector(0.3,0.4,0.2) ∈ UnitSimplex(Val(2)))
    @test z2

    @test approx_in(SA[-0.1,-0.1], d, 0.1)
    @test !approx_in(SA[-0.1,-0.1], d, 0.09)

    @test corners(d) == [ SA[0.0,0.0], SA[1.0,0.0], SA[0.0,1.0]]

    @test convert(Domain{SVector{2,BigFloat}}, d) == EuclideanUnitSimplex{2,BigFloat}()

    @test isclosedset(d)
    @test !isopenset(d)
    @test isopenset(interior(d))
    @test closure(d) == d
    @test point_in_domain(d) ∈ d

    # open/closed
    d2 = EuclideanUnitSimplex{2,Float64,:open}()
    @test !isclosedset(d2)
    @test isopenset(d2)
    @test SA[0.3,0.1] ∈ d2
    @test SA[0.0,0.1] ∉ d2
    @test SA[0.3,0.0] ∉ d2
    @test approx_in(SA[-0.01,0.0], d2, 0.1)
    @test !approx_in(SA[-0.01,0.0], d2, 0.001)

    d3 = EuclideanUnitSimplex{3,BigFloat}()
    @test point_in_domain(d3) ∈ d3
    x0 = big(0.0)
    x1 = big(1.0)
    x2 = big(0.3)
    @test SA[x0,x0,x0] ∈ d3
    @test SA[x1,x0,x0] ∈ d3
    @test SA[x0,x1,x0] ∈ d3
    @test SA[x0,x0,x1] ∈ d3
    @test SA[x2,x0,x0] ∈ d3
    @test SA[x0,x2,x0] ∈ d3
    @test SA[x0,x0,x2] ∈ d3
    @test SA[x2,x2,x2] ∈ d3
    @test SA[-x2,x2,x2] ∉ d3
    @test SA[x2,-x2,x2] ∉ d3
    @test SA[x2,x2,-x2] ∉ d3
    @test SA[x1,x1,x1] ∉ d3

    D = VectorUnitSimplex(2)
    @test isopenset(interior(D))
    @test closure(D) == D
    @test SA[0.2,0.2] ∈ D
    @test SA[0.0,0.2] ∈ D
    @test SA[0.2,0.0] ∈ D
    @test SA[0.5,0.5] ∈ D
    @test SA[0.0,0.0] ∈ D
    @test SA[1.0,0.0] ∈ D
    @test SA[0.0,1.0] ∈ D
    # And then some points outside
    @test SA[0.6,0.5] ∉ D
    @test SA[0.5,0.6] ∉ D
    @test SA[-0.2,0.2] ∉ D
    @test SA[0.2,-0.2] ∉ D
    @test convert(Domain{Vector{BigFloat}}, D) == VectorUnitSimplex{BigFloat}(2)
    @test corners(D) == [ [0.0,0.0], [1.0,0.0], [0.0,1.0]]
    @test boundingbox(D) == UnitCube(4)
end
