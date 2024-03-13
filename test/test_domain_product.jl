using DomainSets: unionbox, intersectbox

function test_product_domains()
    @test productdomain() == ()
    @test productdomain(2) == 2

    @testset "VcatDomain" begin
        d1 = VcatDomain(-1.0..1.0, -1.0..1.0)
        @test d1 isa DomainSets.VcatDomain
        @test d1.domains isa Tuple
        @test eltype(d1) == SVector{2,typeof(1.0)}
        @test SA[0.5,0.5] ∈ d1
        @test SA[-1.1,0.3] ∉ d1
        @test @inferred(VcatDomain(-1.0..1, -1.0..1)) === d1
        # Test promotion
        @test convert(Domain{SVector{2,BigFloat}}, d1) isa ProductDomain{SVector{2,BigFloat}}
        d1w = convert(Domain{SVector{2,BigFloat}}, d1)
        @test eltype(d1w) == SVector{2,BigFloat}
        @test eltype(component(d1w, 1)) == BigFloat
        @test eltype(component(d1w, 2)) == BigFloat

        @test VcatDomain( (-1..1, -2..2)) isa VcatDomain{2,Int,(1,1),Tuple{ClosedInterval{Int64}, ClosedInterval{Int64}}}

        show(io,d1)
        @test String(take!(io)) == "($(-1.0..1.0)) × ($(-1.0..1.0))"

        bnd = boundary(d1)
        @test bnd isa EuclideanDomain
        @test bnd isa UnionDomain
        @test [-1.0, 0.5] ∈ bnd
        @test [1.0, 0.5] ∈ bnd
        @test [0.5, -1.0] ∈ bnd
        @test [0.5, 1.0] ∈ bnd
        @test [0.5, 0.2] ∉ bnd

        # A single domain
        @test VcatDomain(UnitCircle()) isa VcatDomain{2}

        # Test vectors of wrong length
        @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain '($(-1.0..1.0)) × ($(-1.0..1.0))' with dimension 2. Returning false.") SA[0.0,0.0,0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 1 and domain '($(-1.0..1.0)) × ($(-1.0..1.0))' with dimension 2. Returning false.") SA[0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 3 and domain '($(-1.0..1.0)) × ($(-1.0..1.0))' with dimension 2. Returning false.") [0.0,0.0,0.0] ∉ d1
        @test_logs (:warn, "`in`: incompatible combination of vector with length 1 and domain '($(-1.0..1.0)) × ($(-1.0..1.0))' with dimension 2. Returning false.") [0.0] ∉ d1

        d2 = VcatDomain(-1.0 .. 1.0, -1.5 .. 2.5)
        @test SA[0.5,0.5] ∈ d2
        @test SA[-1.1,0.3] ∉ d2

        d3 = VcatDomain(1.05 * UnitDisk(), -1.0 .. 1.0)
        @inferred(cross(1.05 * UnitDisk(), -1.0 .. 1.0)) === d3
        @test d3 isa VcatDomain
        @test eltype(d3) == SVector{3,Float64}
        @test SA[0.5,0.5,0.8] ∈ d3
        @test SA[-1.1,0.3,0.1] ∉ d3
        @test choice(d3) ∈ d3

        v = SA[0.1,0.2,0.3]
        @test map_domain(Translation(v), productdomain(UnitDisk(), UnitInterval())) isa VcatDomain
    end
    @testset "mixed intervals" begin
        d = (0..1) × (0.0..1)
        @test SA[0.1,0.2] ∈ d
        @test SA[0.1,1.2] ∉ d
        @test SA[1.1,1.3] ∉ d
        @test d isa EuclideanDomain{2}
        # Make sure promotion of domains happened
        @test eltype(component(d,1)) == Float64
        @test eltype(component(d,2)) == Float64
        @test choice(d) ∈ d
    end
    @testset "vector domains" begin
        d1 = VectorProductDomain([0..1.0, 0..2.0])
        @test d1 isa VectorDomain{Float64}
        @test d1.domains isa Vector
        @test dimension(d1) == 2
        @test [0.1,0.2] ∈ d1
        @test SA[0.1,0.2] ∈ d1
        @test choice(d1) ∈ d1
        @test convert(Domain{Vector{BigFloat}}, d1) == d1
        d1big = convert(Domain{Vector{BigFloat}}, d1)
        @test eltype(d1big) == Vector{BigFloat}

        # Test an integer type as well
        d2 = VectorProductDomain([0..1, 0..3])
        @test dimension(d2) == 2
        @test [0.1,0.2] ∈ d2
        @test choice(d2) ∈ d2

        # other constructor calls
        @test VectorProductDomain(0..1, 1..2) isa VectorProductDomain{Vector{Int}}
        @test VectorProductDomain((0..1, 1..2)) isa VectorProductDomain{Vector{Int}}
        @test VectorProductDomain{Vector{Int}}(0..1, 1..2) isa VectorProductDomain{Vector{Int}}
        @test VectorProductDomain{Vector{Int}}((0..1, 1..2)) isa VectorProductDomain{Vector{Int}}

        bnd = boundary(d1)
        @test bnd isa VectorDomain
        @test bnd isa UnionDomain
        @test dimension(bnd) == 2
        @test [0.0, 0.5] ∈ bnd
        @test [1.0, 0.5] ∈ bnd
        @test [0.2, 0.0] ∈ bnd
        @test [0.2, 2.0] ∈ bnd
        @test [0.2, 0.5] ∉ bnd

        d3 = VectorProductDomain([0..i for i in 1:10])
        @test d3 isa VectorProductDomain
        @test eltype(d3) == Vector{Int}
        @test rand(10) ∈ d3
        @test 2 .+ rand(10) ∉ d3
        d4 = VectorProductDomain{Vector{Float64}}([0..i for i in 1:10])
        @test d4 isa VectorProductDomain
        @test eltype(d4) == Vector{Float64}
        @test rand(10) ∈ d4
        @test 2 .+ rand(10) ∉ d4

        @test eltype(VectorProductDomain{SVector{2,Float64}}(SVector(0..1, 0..2)).domains[1]) == Float64
    end
    @testset "Tuple product domain" begin
        # Use the constructor ProductDomain{T} directly
        d1 = TupleProductDomain{Tuple{Float64,Float64}}(0..0.5, 0..0.7)
        @test d1 isa TupleProductDomain{Tuple{Float64,Float64}}
        @test d1.domains isa Tuple
        @test eltype(d1) == Tuple{Float64,Float64}
        @test dimension(d1) == 2
        @test (0.2,0.6) ∈ d1
        @test (0.2,0.8) ∉ d1
        @test (true,0.6) ∉ d1
        if VERSION < v"1.6-"
            @test_logs (:warn, "`in`: incompatible combination of point: SArray{Tuple{2},Float64,1,2} and domain eltype: Tuple{Float64,Float64}. Returning false.") SA[0.2,0.6] ∉ d1
            @test_logs (:warn, "`in`: incompatible combination of point: Array{Float64,1} and domain eltype: Tuple{Float64,Float64}. Returning false.") [0.2,0.6] ∉ d1
        else
            @test_logs (:warn, "`in`: incompatible combination of point: SVector{2, Float64} and domain eltype: Tuple{Float64, Float64}. Returning false.") SA[0.2,0.6] ∉ d1
            @test_logs (:warn, "`in`: incompatible combination of point: Vector{Float64} and domain eltype: Tuple{Float64, Float64}. Returning false.") [0.2,0.6] ∉ d1
        end
        @test convert(Domain{Tuple{BigFloat,BigFloat}}, d1) == d1
        d1big = convert(Domain{Tuple{BigFloat,BigFloat}}, d1)
        @test eltype(d1big) == Tuple{BigFloat,BigFloat}
        @test eltype(component(d1big,1)) == BigFloat
        @test eltype(component(d1big,2)) == BigFloat

        d2 = TupleProductDomain(['a','b'], 0..1)
        @test d2.domains isa Tuple
        @test dimension(d2) == 2
        @test eltype(d2) == Tuple{Char,Int}
        @test ('a',0.4) ∈ d2
        @test ('b',1.5) ∉ d2
        @test ('c',0.5) ∉ d2

        # other constructor calls
        @test TupleProductDomain(0..1, 1..2.0) isa TupleProductDomain{Tuple{Int,Float64}}
        @test TupleProductDomain([0..1, 1..2.0]) isa TupleProductDomain{Tuple{Float64,Float64}}
        @test TupleProductDomain{Tuple{Int,Float64}}(0..1, 1..2.0) isa TupleProductDomain{Tuple{Int,Float64}}
        @test TupleProductDomain{Tuple{Float64,Float64}}(0..1, 1..2.0) isa TupleProductDomain{Tuple{Float64,Float64}}
        @test TupleProductDomain{Tuple{Float64,Float64}}([0..1, 1..2.0]) isa TupleProductDomain{Tuple{Float64,Float64}}

        bnd = boundary(d1)
        @test eltype(bnd) == Tuple{Float64,Float64}
        @test bnd isa UnionDomain
        @test dimension(bnd) == 2
        @test (0.0, 0.5) ∈ bnd
        @test (0.5, 0.5) ∈ bnd
        @test (0.3, 0.2) ∉ bnd
        @test (0.3, 0.0) ∈ bnd
        @test (0.3, 0.7) ∈ bnd
    end
    @testset "cube" begin
        @test volume(UnitCube()) == 1
        @test EuclideanUnitCube{2}() == EuclideanUnitCube{2,Float64}()
        @test UnitSquare() == UnitSquare{Float64}()
        @test UnitCube(Val(2)) isa EuclideanUnitCube{2,Float64}
        @test UnitCube{SVector{2,BigFloat}}(Val(2)) isa EuclideanUnitCube{2,BigFloat}

        @test Set(corners((0..1)^2)) == Set([ [1,0], [0,1], [1,1], [0,0]])
        @test Set(corners(UnitCube())) == Set([ [0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,1,0], [1,0,1], [0,1,1], [1,1,1]])
        @test boundary(boundary(UnitSquare())) == UnionDomain(map(t->Point(SVector{2}(t)), corners(UnitSquare())))
        @test boundary(boundary(boundary(UnitCube()))) == UnionDomain(map(t->Point(SVector{3}(t)), corners(UnitCube())))
        @test UnitCube() == UnitCube(3)
        @test boundary(UnitCube()) == boundary(UnitCube(3))
        @test boundary(boundary(UnitCube())) == boundary(boundary(UnitCube(3)))
        @test boundary(boundary(boundary(UnitCube()))) == boundary(boundary(boundary(UnitCube(3))))
        @test UnitSquare() == UnitCube(2)
        @test boundary(UnitSquare()) == boundary(UnitCube(2))
        @test boundary(boundary(UnitSquare())) == boundary(boundary(UnitCube(2)))

        @test StaticUnitCube() == UnitCube()
        @test StaticUnitCube(Val(3)) == UnitCube()
        @test convert(Domain{SVector{2,Float64}}, UnitCube()) == StaticUnitCube(Val(2))
        @test convert(Domain{Vector{Float64}}, UnitCube()) === VectorUnitCube(3)

        @test DynamicUnitCube{SVector{3,Float64}}(3) isa DynamicUnitCube
        @test_throws AssertionError DynamicUnitCube{SVector{3,Float64}}(2)

        d1 = VectorUnitCube{Float64}(4)
        @test VectorUnitCube(4) == d1
        @test UnitCube(4) == d1
        @test UnitCube{Vector{Float64}}(4) == d1
        @test dimension(d1) == 4
        @test component(d1, 1) == 0..1
        @test SA[0.9,0.9,0.4,0.2] ∈ d1
        @test [1.2,0.3,0.4,0.6] ∉ d1
        @test convert(Domain{SVector{4,Float64}}, d1) isa StaticUnitCube
        @test convert(Domain{Vector{BigFloat}}, d1) === VectorUnitCube{BigFloat}(4)
        @test DomainSets.VectorUnitSquare() === UnitCube(2)

        @test ProductDomain([UnitInterval(),UnitInterval()]) isa VectorUnitCube{Float64}
        @test ProductDomain{Vector{Float64}}([UnitInterval(),UnitInterval()]) isa VectorUnitCube{Float64}
        @test ProductDomain{Vector{BigFloat}}([UnitInterval(),UnitInterval()]) isa VectorUnitCube{BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(UnitInterval(),UnitInterval()) isa EuclideanUnitCube{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}((UnitInterval(),UnitInterval())) isa EuclideanUnitCube{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(SVector(UnitInterval(),UnitInterval())) isa EuclideanUnitCube{2,BigFloat}

        @test repr(UnitCube()) == "UnitCube()"
        @test repr(UnitCube(Val(4))) == "UnitCube(Val(4))"
        @test repr(UnitSquare()) == "UnitSquare()"
        @test repr(UnitSquare{BigFloat}()) == "UnitSquare{BigFloat}()"
        @test repr(UnitCube(4)) == "UnitCube(4)"

        D = UnitInterval()^2
        @test D isa EuclideanUnitCube{2,Float64}
        @test SA[0.9, 0.9] ∈ D
        @test SA[1.1, 1.1] ∉ D
        @test !isempty(D)
        @test isclosedset(D)
        @test !isopenset(D)
        @test convert(DomainSets.HyperRectangle, VcatDomain(UnitInterval(), UnitInterval())) isa StaticUnitCube

        @test approx_in(SA[-0.1,-0.1], D, 0.1)
        @test !approx_in(SA[-0.1,-0.1], D, 0.09)

        @test intersect(ChebyshevInterval()^3, UnitInterval()^3) isa EuclideanUnitCube{3,Float64}

        #Cube
        D = (-1.5 .. 2.2) × (0.5 .. 0.7) × (-3.0 .. -1.0)
        @test SA[0.9, 0.6, -2.5] ∈ D
        @test SA[0.0, 0.6, 0.0] ∉ D

        @test issubset( (0..1)^3, (-1..2)^3 )
    end
    @testset "Rectangle" begin
        d1 = (-1.0..1.0) × (-1.0..1.0)

        d4 = d1 × (-1.0..1.0)
        @test d4 isa Rectangle
        @test SA[0.5,0.5,0.8] ∈ d4
        @test SA[-1.1,0.3,0.1] ∉ d4
        @test choice(d4) ∈ d4
        @test center(d4) == [0.0,0.0,0.0]

        @test d1[Component(1)] == -1..1
        @test d1[Component(2)] == -1..1
        @test_throws BoundsError d1[Component(3)]

        d5 = (-1.0..1.)×d1
        @test d5 isa Rectangle
        @test SA[0.,0.5,0.5] ∈ d5
        @test SA[0.,-1.1,0.3] ∉ d5
        @test choice(d5) ∈ d5

        d6 = d1 × d1
        @test d6 isa Rectangle
        @test SA[0.,0.,0.5,0.5] ∈ d6
        @test SA[0.,0.,-1.1,0.3] ∉ d6
        @test choice(d6) ∈ d6

        @test Rectangle( SA[1,2], SA[2.0,3.0]) isa Rectangle{SVector{2,Float64}}
        @test Rectangle([0..1, 2..3]) isa Rectangle{Vector{Int}}
        @test Rectangle((0..1, 2..3)) isa Rectangle{SVector{2,Int}}
        @test Rectangle{SVector{2,Float64}}((0..1, 2..3)) isa Rectangle{SVector{2,Float64}}

        @test_throws ArgumentError Rectangle(UnitCircle(), UnitDisk())
        @test_throws ArgumentError Rectangle(OpenInterval(1,2), 3..4)
        @test_throws ArgumentError Rectangle{SVector{2,Float64}}(UnitCircle(), UnitDisk())
        @test_throws ArgumentError Rectangle(0.0, 1.0)

        bnd = boundary(Rectangle([1,2],[3,4]))
        @test [1,3] ∈ bnd
        @test [1,2.5] ∈ bnd
        @test [1.5,4] ∈ bnd
        @test [1.5,3.5] ∉ bnd

        @test 2 .* UnitCube() isa Rectangle{SVector{3,Float64}}
        @test 2 .* UnitCube() .+ SA[1,2,3.0] isa Rectangle{SVector{3,Float64}}
        @test 2 .* UnitCube() .+ SA[1,2,3.0] == ProductDomain(1..3, 2..4, 3..5)
    end
    @testset "fixed product domains" begin
        d1 = ProductDomain(ChebyshevInterval(), ChebyshevInterval())
        @test d1 isa DomainSets.FixedIntervalProduct
        @test d1 isa DomainSets.ChebyshevProductDomain
        @test volume(d1) == 4
        @test component(d1,1) isa ChebyshevInterval{Float64}
        @test components(d1) isa NTuple{2,ChebyshevInterval{Float64}}
        @test convert(Domain{SVector{2,BigFloat}}, d1) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
        @test d1 == DomainSets.ChebyshevProductDomain(Val(2))
        @test d1 == DomainSets.ChebyshevProductDomain{2}()
        @test ProductDomain(components(d1)) == d1
        @test ProductDomain{SVector{2,BigFloat}}(components(d1)) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
        @test ProductDomain{SVector{2,BigFloat}}(components(d1)...) isa DomainSets.ChebyshevProductDomain{2,BigFloat}
    end
    @testset "ProductDomain" begin
        d1 = 0..1.0
        d2 = UnitBall{Float64}()
        d3 = UnitCircle()
        @test ProductDomain{SVector{2,Float64}}(d1, d2) isa VcatDomain
        @test ProductDomain{Tuple{Float64,Float64}}(d1, d2) isa TupleProductDomain
        @test ProductDomain{Vector{Float64}}([d1; d2]) isa VectorProductDomain

        @test ProductDomain((d1,d3)) isa VcatDomain
        @test ProductDomain((d1,d1)) isa Rectangle

        @test volume(ProductDomain(0..1,1..3)) == 2

        @test ProductDomain(SVector(0..1, 0..2)) isa Rectangle{SVector{2,Int}}
        @test ProductDomain(1.05 * UnitDisk(), -1.0 .. 1.0) isa VcatDomain{3,Float64}
        @test ProductDomain(['a','b'], 0..1) isa TupleProductDomain

        @test ProductDomain{Tuple{Float64,Float64}}(0..0.5, 0..0.7) isa Rectangle{Tuple{Float64,Float64}}
        # Some conversions
        @test convert(Domain{Vector{Float64}}, productdomain(d1,d2)) isa VectorProductDomain{Vector{Float64}}
        @test convert(Domain{SVector{2,Float64}}, ProductDomain([d1,d2])) isa VcatDomain{2,Float64}
        @test convert(Domain{Vector{Float64}}, TupleProductDomain(d1,d2)) isa VectorDomain{Float64}

        # intersection of product domains
        @test ProductDomain([0..1.0, 0..2.0]) ∩ ProductDomain([0..1, 0..3]) isa Rectangle{Vector{Float64}}
        @test ProductDomain(0..1.0, 0..2.0) ∩ ProductDomain(0..1, 0..3) isa Rectangle{SVector{2,Float64}}

        # Generic functionality
        long_domain = ProductDomain([0..i for i in 1:20])
        show(io, long_domain)
        @test String(take!(io)) == "($(0..1)) × ($(0..2)) × ($(0..3)) × ... × ($(0..20))"
        @test isopenset(interior(UnitCube()))
        @test isclosedset(closure(interior(UnitCube())))
    end
    @testset "bounding box" begin
        @test boundingbox([0.2, -0.4, 1.0]) == -0.4..1.0
        @test boundingbox(Set([0.2, -0.4, 1.0])) == -0.4..1.0


        @test unionbox(0..1) == 0..1
        @test unionbox(0..1, 2..3) == 0..3
        @test unionbox(0..1, 2..3, 5..6.0) === 0..6.0
        @test unionbox(ChebyshevInterval(), ChebyshevInterval()) === ChebyshevInterval()

        @test unionbox(EmptySpace{Float64}(), 2..3) === 2..3.0
        @test unionbox(2..3, EmptySpace{Float64}()) === 2..3.0
        @test unionbox(FullSpace{Int}(), 2..3.0) === FullSpace{Float64}()
        @test unionbox(2..3.0, FullSpace{Int}()) === FullSpace{Float64}()

        @test DomainSets.unionbox2(0..1, 0..1) == FullSpace()

        d1 = ProductDomain(0..1.0, 2..4.0)
        d2 = ProductDomain(0.5..1.5, 2.5..4.5)
        @test unionbox(d1, d2) == ProductDomain(0..1.5, 2..4.5)

        @test intersectbox(0..1) == 0..1
        @test intersectbox(0..1.5, 1..3) === 1.0..1.5
        @test intersectbox(0..1, 0.5..3, 0.0..4.0) === 0.5..1.0
        @test intersectbox(ChebyshevInterval(), ChebyshevInterval()) === ChebyshevInterval()

        @test intersectbox(EmptySpace{Int}(), 2..3.0) === EmptySpace{Float64}()
        @test intersectbox(2..3.0, EmptySpace{Int}()) === EmptySpace{Float64}()
        @test intersectbox(FullSpace{Float64}(), 2..3) === 2.0..3.0
        @test intersectbox(2..3, FullSpace{Float64}()) === 2.0..3.0

        @test DomainSets.intersectbox2(0..1, 0..1) == FullSpace()

        d1 = ProductDomain(0..1.0, 2..4.0)
        d2 = ProductDomain(0.5..1.5, 2.5..4.5)
        @test intersectbox(d1, d2) == ProductDomain(0.5..1.0, 2.5..4.0)

        @test boundingbox(MappedDomain(LinearMap(1/2), 2..3)) == 4.0..6.0
        @test boundingbox(MappedDomain(LinearMap(1/2), (2..3)^2)) == (4.0..6.0)^2
    end
    @testset "product mapto" begin
        d1 = ProductDomain(1.0..2.0, 1.0..2.0)
        d2 = ProductDomain(2.0..4.0, 2.0..4.0)
        @test mapto(d1, d2) isa DomainSets.VcatMap
        @test jacobian(mapto(d1,d2)) isa ConstantMap
        @test jacdet(mapto(d1,d2)) == DomainSets.ConstantMap{SVector{2,Float64}}(4.0)
        d1v = ProductDomain([1.0..2.0, 1.0..2.0])
        d2v = ProductDomain([2.0..4.0, 2.0..4.0])
        @test mapto(d1v, d2v) isa DomainSets.VectorProductMap
        @test jacobian(mapto(d1v,d2v)) isa DomainSets.ConstantMap
        @test jacdet(mapto(d1v,d2v)) == DomainSets.ConstantMap{Vector{Float64}}(4.0)
    end

    @testset "VcatDomain == bug" begin
        x = VcatDomain(0.0:0.5:2.0, [1,3,4])
        y = VcatDomain(0.0:0.5:2.0, [1,3,4])
        @test x == y
        @test isequaldomain(5, [5.0])
        @test isequaldomain(Set([5.0]), [5.0])
    end
end
