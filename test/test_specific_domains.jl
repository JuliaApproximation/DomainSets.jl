# test_specific_domains.jl

const v = TypeFactory{SVector}()

function test_specific_domains()
    @testset "$(rpad("Specific domains",80))" begin
        test_emptyspace()
        test_fullspace()
        test_interval()
        test_unitball()
        test_cube()
        test_simplex()
        test_sphere()
        test_arithmetics()
        test_tensorproduct_domain()
    end
end

function test_emptyspace()
    println("- an empty space")
    d1 = EmptySpace()
    @test eltype(d1) == Float64
    @test 0.5 ∉ d1
    @test d1 ∩ d1 == d1
    @test d1 ∪ d1 == d1
    d2 = interval()
    @test d1 ∩ d2 == d1
    @test d1 ∪ d2 == d2

    d2 = EmptySpace(SVector{2,Float64})
    @test v[0.1,0.2] ∉ d2
end

function test_fullspace()
    println("- a full Euclidean space")
    d1 = FullSpace()
    @test 0.5 ∈ d1
    @test d1 ∪ d1 == d1
    @test d1 ∩ d1 == d1
    d2 = interval()
    @test d1 ∪ d2 == d1
    @test d1 ∩ d2 == d2

    d2 = FullSpace(SVector{2,Float64})
    @test v[0.1,0.2] ∈ d2
end

function test_interval(T = Float64)
    println("- intervals")

    d = interval(zero(T), one(T))
    @test T(0.5) ∈ d
    @test T(1.1) ∉ d
    @test leftendpoint(d) == zero(T)
    @test rightendpoint(d) == one(T)
    @test isclosed(d)
    @test !isopen(d)
    @test iscompact(d)
    @test typeof(similar_interval(d, one(T), 2*one(T))) == typeof(d)

    d = UnitInterval{T}()
    @test leftendpoint(d) == zero(T)
    @test rightendpoint(d) == one(T)
    @test isclosed(d)
    @test !isopen(d)
    @test iscompact(d)

    d = ChebyshevInterval{T}()
    @test leftendpoint(d) == -one(T)
    @test rightendpoint(d) == one(T)
    @test isclosed(d)
    @test !isopen(d)
    @test iscompact(d)

    d = halfline(T)
    @test leftendpoint(d) == zero(T)
    @test rightendpoint(d) == T(Inf)
    @test !isclosed(d)
    @test !isopen(d)
    @test !iscompact(d)

    d = negative_halfline(T)
    @test leftendpoint(d) == -T(Inf)
    @test rightendpoint(d) == zero(T)
    @test !isclosed(d)
    @test isopen(d)
    @test !iscompact(d)

    ## Some mappings preserve the interval structure
    # Translation
    d = interval(zero(T), one(T))
    d2 = d + one(T)
    @test typeof(d2) == typeof(d)
    @test leftendpoint(d2) == one(T)
    @test rightendpoint(d2) == 2*one(T)

    # Does a shifted unit interval return an interval?
    d = UnitInterval{T}()
    d2 = d + one(T)
    @test typeof(d2) <: AbstractInterval

    # Scaling
    d = interval(zero(T), one(T))
    d3 = T(2) * d
    @test typeof(d3) == typeof(d)
    @test leftendpoint(d3) == zero(T)
    @test rightendpoint(d3) == T(2)

    d = interval(zero(T), one(T))
    d4 = d / T(2)
    @test typeof(d4) == typeof(d)
    @test leftendpoint(d4) == zero(T)
    @test rightendpoint(d4) == T(1)/T(2)

    # Union and intersection of intervals
    i1 = interval(zero(T), one(T))
    i2 = interval(one(T)/3, one(T)/2)
    i3 = interval(one(T)/2, 2*one(T))
    i4 = interval(T(2), T(3))
    # - union of completely overlapping intervals
    du1 = i1 ∪ i2
    @test typeof(du1) <: AbstractInterval
    @test leftendpoint(du1) == leftendpoint(i1)
    @test rightendpoint(du1) == rightendpoint(i1)

    # - intersection of completely overlapping intervals
    du2 = i1 ∩ i2
    @test typeof(du2) <: AbstractInterval
    @test leftendpoint(du2) == leftendpoint(i2)
    @test rightendpoint(du2) == rightendpoint(i2)

    # - union of partially overlapping intervals
    du3 = i1 ∪ i3
    @test typeof(du3) <: AbstractInterval
    @test leftendpoint(du3) == leftendpoint(i1)
    @test rightendpoint(du3) == rightendpoint(i3)

    # - intersection of partially overlapping intervals
    du4 = i1 ∩ i3
    @test typeof(du4) <: AbstractInterval
    @test leftendpoint(du4) == leftendpoint(i3)
    @test rightendpoint(du4) == rightendpoint(i1)

    # - union of non-overlapping intervals
    du5 = i1 ∪ i4
    @test typeof(du5) <: UnionDomain

    # - intersection of non-overlapping intervals
    du6 = i1 ∩ i4
    @test typeof(du6) == EmptySpace{T}
end

function test_unitball()
    C = disk(2.0)
    @test in(v[1.4, 1.4], C)
    @test !in(v[1.5, 1.5], C)
    @test typeof(1.2*C)==typeof(C*1.2)
    @test in(v[1.5,1.5],1.2*C)
    @test in(v[1.5,1.5],C*1.2)

    S = ball(2.0)
    @test v[1.9,0.0,0.0] ∈ S
    @test in(v[0,-1.9,0.0],S)
    @test in(v[0.0,0.0,-1.9],S)
    @test !in(v[1.9,1.9,0.0],S)
end

function test_cube()
    #Square
    D = cube(Val{2})
    @test v[0.9, 0.9] ∈ D
    @test v[1.1, 1.1] ∉ D

    #Cube
    D = cube(-1.5, 2.2, 0.5, 0.7, -3.0, -1.0)
    @test v[0.9, 0.6, -2.5] ∈ D
    @test v[0.0, 0.6, 0.0] ∉ D
end

function test_simplex()
    d = simplex(Val{2})
    # We test a point in the interior, a point on each of the boundaries and
    # all corners.
    @test v[0.2,0.2] ∈ d
    @test v[0.0,0.2] ∈ d
    @test v[0.2,0.0] ∈ d
    @test v[0.5,0.5] ∈ d
    @test v[0.0,0.0] ∈ d
    @test v[1.0,0.0] ∈ d
    @test v[0.0,1.0] ∈ d
    # And then some points outside
    @test v[0.6,0.5] ∉ d
    @test v[0.5,0.6] ∉ d
    @test v[-0.2,0.2] ∉ d
    @test v[0.2,-0.2] ∉ d

    d3 = simplex(Val{3}, BigFloat)
    x0 = big(0.0)
    x1 = big(1.0)
    x2 = big(0.3)
    @test v[x0,x0,x0] ∈ d3
    @test v[x1,x0,x0] ∈ d3
    @test v[x0,x1,x0] ∈ d3
    @test v[x0,x0,x1] ∈ d3
    @test v[x2,x0,x0] ∈ d3
    @test v[x0,x2,x0] ∈ d3
    @test v[x0,x0,x2] ∈ d3
    @test v[x2,x2,x2] ∈ d3
    @test v[-x2,x2,x2] ∉ d3
    @test v[x2,-x2,x2] ∉ d3
    @test v[x2,x2,-x2] ∉ d3
    @test v[x1,x1,x1] ∉ d3
end

function test_sphere()
end

function test_arithmetics()
    # joint domain
    D = cube(Val{3})
    S = ball(2.0)
    DS = D ∪ S
    @test v[0.0, 0.6, 0.0] ∈ DS
    @test v[0.9, 0.6,-2.5] ∉ DS

    # domain intersection
    DS = D ∩ S
    @test v[0.1, 0.1, 0.1] ∈ DS
    @test v[0.1, -0.1, 0.1] ∉ DS

    # domain difference
    DS = D-S
    @test v[0.1, 0.1, 0.1] ∉ DS
end

function test_tensorproduct_domain()
    # ProductDomain 1
    T = tensorproduct(interval(-1.0, 1.0), 2)
    @test v[0.5,0.5] ∈ T
    @test v[-1.1,0.3] ∉ T

    # ProductDomain 2
    T = ProductDomain(disk(1.05), interval(-1.0,1.0))
    @test (v[0.5,0.5,0.8]) ∈ T
    @test (v[-1.1,0.3,0.1]) ∉ T
end
