function test_balls()
    # Test UnitBall constructor
    @test UnitBall(2) isa VectorUnitBall{Float64}
    @test UnitBall(Val(2)) isa EuclideanUnitBall{2,Float64}
    @test UnitBall{Float64}() isa StaticUnitBall{Float64}
    @test UnitBall{Float64}(1) isa StaticUnitBall{Float64}
    @test_throws AssertionError UnitBall{Float64}(2)
    @test UnitBall{SVector{2,Float64}}(Val(2)) isa EuclideanUnitBall{2,Float64}
    @test_throws AssertionError UnitBall{SVector{2,Float64}}(Val(3))
    @test UnitBall{SVector{2,Float64}}(2) isa EuclideanUnitBall{2,Float64}
    @test_throws AssertionError UnitBall{SVector{2,Float64}}(3)
    @test UnitBall{SVector{2,Float64}}() isa EuclideanUnitBall{2,Float64}
    @test UnitBall{Vector{Float64}}(2) isa VectorUnitBall{Float64}
    @test_throws MethodError UnitBall{Vector{Float64}}()

    @test UnitBall{Float64,:open}() isa StaticUnitBall{Float64,:open}
    @test UnitBall{Float64,:closed}(1) isa StaticUnitBall{Float64,:closed}
    @test_throws AssertionError UnitBall{Float64,:closed}(2)
    @test UnitBall{SVector{2,Float64},:open}(Val(2)) isa EuclideanUnitBall{2,Float64,:open}
    @test_throws AssertionError UnitBall{SVector{2,Float64},:closed}(Val(3))
    @test UnitBall{SVector{2,Float64},:open}(2) isa EuclideanUnitBall{2,Float64,:open}
    @test_throws AssertionError UnitBall{SVector{2,Float64},:closed}(3)
    @test UnitBall{SVector{2,Float64},:open}() isa EuclideanUnitBall{2,Float64,:open}
    @test UnitBall{Vector{Float64},:closed}(2) isa VectorUnitBall{Float64,:closed}
    @test_throws MethodError UnitBall{Vector{Float64},:open}()

    @test DynamicUnitBall{Float64}(1) isa DynamicUnitBall{Float64}
    @test_throws AssertionError DynamicUnitBall{Float64}(2)

    # and then the Ball constructor
    @test Ball() isa UnitBall
    @test Ball(1.0, 2.0) isa DomainSets.GenericBall{Float64}
    @test Ball(1.0, Point(2.0)) isa DomainSets.GenericBall{Float64}
    @test Ball{Float64}() isa UnitBall{Float64}
    @test Ball{Float64,:open}() isa UnitBall{Float64,:open}
    @test Ball(1.0) == DomainSets.GenericBall{SVector{3,Float64},:closed}(1.0)
    @test hash(Ball(1.0)) == hash(DomainSets.GenericBall{SVector{3,Float64},:closed}(1.0))
    @test Ball{BigFloat}(1) == DomainSets.GenericBall{BigFloat,:closed}(big(1))
    @test hash(Ball{BigFloat}(1)) == hash(DomainSets.GenericBall{BigFloat,:closed}(big(1)))
    @test Ball{BigFloat,:open}(1) == DomainSets.GenericBall{BigFloat,:open}(big(1))
    @test hash(Ball{BigFloat,:open}(1)) == hash(DomainSets.GenericBall{BigFloat,:open}(big(1)))
    @test_throws MethodError Ball{Vector{Float64}}(1.0)
    @test Ball(1.0, [1,2,3]) isa DomainSets.GenericBall{Vector{Float64},:closed,Float64}
    @test Ball{Vector{Float64}}(1.0, [1,2,3]) isa DomainSets.GenericBall{Vector{Float64},:closed,Float64}
    @test Ball{Vector{Float64}}(1.0, Point([1,2,3])) isa DomainSets.GenericBall{Vector{Float64},:closed,Float64}

    # the Disk constructor
    @test typeof(Disk()) == EuclideanUnitBall{2,Float64,:closed}
    @test Disk(2) isa Ball{SVector{2,Float64}}
    @test Disk(2.0) isa Ball{SVector{2,Float64}}
    @test Disk(2, SA[0.0,1.0]) isa Ball{SVector{2,Float64}}

    @test map_domain(AffineMap(2, [1,2,3]), UnitBall()) isa DomainSets.GenericBall
    @test map_domain(AffineMap(2, 3), UnitBall{Float64}()) isa DomainSets.GenericBall
    @test map_domain(LinearMap{SVector{3,Float64}}(2), UnitBall()) isa DomainSets.GenericBall
    @test map_domain(LinearMap(2), UnitBall{Float64}()) isa DomainSets.GenericBall
    @test map_domain(LinearMap(2), UnitBall()) isa DomainSets.GenericBall
    @test map_domain(Translation([1,2,3]), UnitBall()) isa DomainSets.GenericBall
    @test map_domain(Translation(1), UnitBall{Float64}()) isa DomainSets.GenericBall

    @test StaticUnitBall() isa StaticUnitBall{Float64}
    @test StaticUnitBall(Val(2)) isa StaticUnitBall{SVector{2,Float64}}

    @test GenericBall() == GenericBall(1.0)
    @test hash(GenericBall()) == hash(GenericBall(1.0))
    @test GenericBall(2.0, 1:5) isa GenericBall{Vector{Float64},:closed,Float64}
    @test GenericBall(2, 1.0:5.0) isa GenericBall{Vector{Float64},:closed,Float64}

    # set difference with mixed vector types
    @test eltype(Disk(0.8)\Rectangle([-0.15,-1.0],[0.15,0.0])) == AbstractVector{Float64}

    @test issubset(Ball(1.0, 2.0), Ball(1.0, 2.0))
    @test issubset(Ball(1.0, 2.0), Ball(2.0, 2.0))
    @test !issubset(Ball(2.0, 2.0), Ball(1.0, 2.0))
    @test issubset(Ball{Float64,:open}(2.0, 2.0), Ball(2.0, 2.0))
    @test !issubset(Ball(2.0, 2.0), Ball{Float64,:open}(2.0, 2.0))
    @test issubset(Ball(1.0, 0.0), Ball(2.5, 1.0))
    @test !issubset(Ball(1.0, 0.0), Ball(0.5, 1.0))
    @test !issubset(UnitBall(Val(2)), UnitBall(Val(3)))

    D = UnitDisk()
    @test SA[1.,0.] ∈ D
    @test SA[1.,1.] ∉ D
    @test approx_in(SA[1.0,0.0+1e-5], D, 1e-4)
    @test !isempty(D)
    @test isclosedset(D)
    @test !isopenset(D)
    D2 = convert(Domain{SVector{2,BigFloat}}, D)
    @test eltype(D2) == SVector{2,BigFloat}
    @test boundary(D) == UnitCircle()
    @test dimension(D) == 2
    @test boundingbox(D) == ProductDomain(ChebyshevInterval(), ChebyshevInterval())
    @test normal(UnitDisk(), [sqrt(2)/2, sqrt(2)/2]) ≈ [sqrt(2)/2, sqrt(2)/2]

    @test boundingbox(UnitBall{Float64}()) == ChebyshevInterval()
    @test boundingbox(Ball(1.0, 1.0)) == 0.0..2.0 # issue #113
    @test boundingbox(Ball(1.0, SA[BigFloat(1), BigFloat(2)])) == Rectangle([0.0, 1.0], [2.0, 3.0]) # issue #114

    @test convert(SublevelSet, UnitDisk()) isa SublevelSet{SVector{2,Float64},:closed}
    @test convert(SublevelSet{SVector{2,Float64}}, UnitDisk()) isa SublevelSet{SVector{2,Float64},:closed}
    @test convert(SublevelSet, EuclideanUnitBall{2,Float64,:open}()) isa SublevelSet{SVector{2,Float64},:open}

    @test convert(Interval, UnitBall{Float64}()) === ChebyshevInterval()
    @test convert(Interval, UnitBall{Float64,:open}()) === OpenInterval(-1.0, 1.0)
    @test UnitBall{Float64}() == ChebyshevInterval()

    @test convert(Domain{SVector{2,Float64}}, UnitBall(2)) isa StaticUnitBall
    @test convert(Domain{Vector{Float64}}, UnitBall(Val(2))) isa DynamicUnitBall

    @test repr(UnitBall()) == "UnitBall()"
    @test repr(UnitBall(Val(4))) == "UnitBall(Val(4))"
    @test repr(EuclideanUnitBall{3,Float64,:open}()) == "UnitBall()  (open)"
    @test repr(EuclideanUnitBall{4,Float64,:open}()) == "UnitBall(Val(4))  (open)"
    @test repr(UnitDisk()) == "UnitDisk()"
    @test repr(UnitDisk{BigFloat}()) == "UnitDisk{BigFloat}()"
    @test repr(UnitBall{Float64}()) == "UnitBall{Float64}()"
    @test repr(UnitBall{Float64,:open}()) == "UnitBall{Float64}()  (open)"
    @test repr(VectorUnitBall{Float64}(4)) == "UnitBall(4)"
    @test repr(VectorUnitBall{Float64,:open}(4)) == "UnitBall(4)  (open)"
    @test repr(Ball(1.0,2.0)) == "Ball(1.0, 2.0)"

    D = EuclideanUnitBall{2,Float64,:open}()
    @test !in(SA[1.0,0.0], D)
    @test in(SA[1.0-eps(Float64),0.0], D)
    @test approx_in(SA[1.1,0.0], D, 0.2)
    @test !approx_in(SA[1.1,0.0], D, 0.01)
    @test SA[0.2,0.2] ∈ D
    @test !isclosedset(D)
    @test isopenset(D)
    @test boundary(D) == UnitCircle()

    D = 2UnitDisk()
    @test D isa DomainSets.GenericBall
    @test SA[1.4, 1.4] ∈ D
    @test SA[1.5, 1.5] ∉ D
    @test typeof(1.2 * D)==typeof(D * 1.2)
    @test SA[1.5,1.5] ∈ 1.2 * D
    @test SA[1.5,1.5] ∈ D * 1.2
    @test !isempty(D)

    @test canonicaldomain(D) == UnitDisk()
    @test affinematrix(mapfrom_canonical(D)) == [2 0; 0 2]
    @test affinevector(mapfrom_canonical(D)) == [0; 0]
    @test parameterdomain(D) == UnitSquare()
    @test parameterdomain(UnitBall(2)) == UnitSquare()
    @test parameterdomain(UnitBall(3)) == UnitBall(3)
    @test mapfrom_parameterdomain(D)([0.2,0.4]) ∈ D
    @test mapfrom_parameterdomain(UnitBall(2))([0.2,0.4]) ∈ UnitBall(2)
    @test boundingbox(D) == (-2.0..2.0)^2

    D = 2UnitDisk() .+ SA[1.0,1.0]
    @test D isa DomainSets.GenericBall
    @test boundary(D) isa DomainSets.GenericSphere
    @test SA[2.4, 2.4] ∈ D
    @test SA[3.5, 2.5] ∉ D
    @test !isempty(D)
    @test isopenset(interior(D))
    @test isclosedset(closure(D))
    @test canonicaldomain(D) isa UnitDisk

    B = UnitBall()
    @test SA[1.,0.0,0.] ∈ B
    @test SA[1.,0.1,0.] ∉ B
    @test !isempty(B)
    @test isclosedset(B)
    @test !isopenset(B)
    @test boundary(B) == UnitSphere()
    @test isopenset(interior(B))
    @test isclosedset(closure(B))

    B = 2UnitBall()
    @test D isa DomainSets.GenericBall
    @test SA[1.9,0.0,0.0] ∈ B
    @test SA[0,-1.9,0.0] ∈ B
    @test SA[0.0,0.0,-1.9] ∈ B
    @test SA[1.9,1.9,0.0] ∉ B
    @test !isempty(B)

    B = 2.0UnitBall() .+ SA[1.0,1.0,1.0]
    @test D isa DomainSets.GenericBall
    @test SA[2.9,1.0,1.0] ∈ B
    @test SA[1.0,-0.9,1.0] ∈ B
    @test SA[1.0,1.0,-0.9] ∈ B
    @test SA[2.9,2.9,1.0] ∉ B
    @test !isempty(B)

    C = VectorUnitBall(4)
    @test [1, 0, 0, 0] ∈ C
    @test [0.0,0.1,0.2,0.1] ∈ C
    @test SA[0.0,0.1,0.2,0.1] ∈ C
    @test_logs (:warn, "`in`: incompatible combination of vector with length 2 and domain 'UnitBall(4)' with dimension 4. Returning false.") [0.0,0.1] ∉ C
    @test [0.0,1.1,0.2,0.1] ∉ C
    @test !isempty(C)
    @test isclosedset(C)
    @test !isopenset(C)
    @test boundary(C) == VectorUnitSphere(4)
    cheb = ChebyshevInterval()
    @test boundingbox(C) == ProductDomain([cheb, cheb, cheb, cheb])
    @test dimension(C) == 4
    @test isopenset(interior(C))

    D = VectorUnitBall{Float64,:open}(4)
    @test !in([1, 0, 0, 0], D)
    @test in([1-eps(Float64), 0, 0, 0], D)
    @test approx_in([1.1, 0, 0, 0], D, 0.2)
    @test !approx_in([1.1, 0, 0, 0], D, 0.01)
    @test !isempty(D)
    @test approx_in(SA[1.01,0.0,0.0,0.0], D, 0.05)

    E = Ball{Float64,:open}(2.0)
    @test 0.5 ∈ E
    @test approx_in(0.5, E)
    @test approx_in(0.5, closure(E))
    @test isclosedset(StaticUnitBall{SVector{2,Float64}}())
    @test EuclideanUnitBall{2}() isa EuclideanUnitBall{2,Float64}

    show(io, EuclideanUnitBall{2,Float64,:open}())
    @test String(take!(io)) == "UnitBall(Val(2))  (open)"
    show(io, UnitCircle())
    @test String(take!(io)) == "UnitCircle()"
end

function test_spheres()
    # test UnitSphere constructor
    @test UnitSphere(2) isa VectorUnitSphere{Float64}
    @test UnitSphere(Val(2)) isa EuclideanUnitSphere{2,Float64}
    @test UnitSphere{Float64}() isa StaticUnitSphere{Float64}
    @test UnitSphere{Float64}(1) isa StaticUnitSphere{Float64}
    @test_throws AssertionError UnitSphere{Float64}(2)
    @test UnitSphere{SVector{2,Float64}}(Val(2)) isa EuclideanUnitSphere{2,Float64}
    @test_throws AssertionError UnitSphere{SVector{2,Float64}}(Val(3))
    @test UnitSphere{SVector{2,Float64}}(2) isa EuclideanUnitSphere{2,Float64}
    @test_throws AssertionError UnitSphere{SVector{2,Float64}}(3)
    @test UnitSphere{SVector{2,Float64}}() isa EuclideanUnitSphere{2,Float64}
    @test UnitSphere{Vector{Float64}}(2) isa VectorUnitSphere{Float64}
    @test_throws MethodError UnitSphere{Vector{Float64}}()

    @test DynamicUnitSphere{Float64}(1) isa DynamicUnitSphere{Float64}
    @test_throws AssertionError DynamicUnitSphere{Float64}(2)

    # and the generic sphere constructor
    @test Sphere() isa UnitSphere
    @test Sphere(1.0, 2.0) isa DomainSets.GenericSphere{Float64}
    @test Sphere(1.0, Point(2.0)) isa DomainSets.GenericSphere{Float64}
    @test Sphere{Float64}() isa UnitSphere{Float64}
    @test Sphere(1.0) == DomainSets.GenericSphere{SVector{3,Float64}}(1.0)
    @test hash(Sphere(1.0)) == hash(DomainSets.GenericSphere{SVector{3,Float64}}(1.0))
    @test Sphere{BigFloat}(1) == DomainSets.GenericSphere{BigFloat}(big(1))
    @test hash(Sphere{BigFloat}(1)) == hash(DomainSets.GenericSphere{BigFloat}(big(1)))
    @test_throws MethodError Sphere{Vector{Float64}}(1.0)
    @test Sphere(1.0, [1,2,3]) isa DomainSets.GenericSphere{Vector{Float64},Float64}
    @test Sphere{Vector{Float64}}(1.0, [1,2,3]) isa DomainSets.GenericSphere{Vector{Float64},Float64}
    @test Sphere{Vector{Float64}}(1.0, Point([1,2,3])) isa DomainSets.GenericSphere{Vector{Float64},Float64}

    @test map_domain(AffineMap(2, [1,2,3]), UnitSphere()) isa DomainSets.GenericSphere
    @test map_domain(AffineMap(2, 3), UnitSphere{Float64}()) isa DomainSets.GenericSphere
    @test map_domain(LinearMap{SVector{3,Float64}}(2), UnitSphere()) isa DomainSets.GenericSphere
    @test map_domain(LinearMap(2), UnitSphere{Float64}()) isa DomainSets.GenericSphere
    @test map_domain(LinearMap(2), UnitSphere()) isa DomainSets.GenericSphere
    @test map_domain(Translation([1,2,3]), UnitSphere()) isa DomainSets.GenericSphere
    @test map_domain(Translation(1), UnitSphere{Float64}()) isa DomainSets.GenericSphere

    @test StaticUnitSphere() isa StaticUnitSphere{Float64}
    @test StaticUnitSphere(Val(2)) isa StaticUnitSphere{SVector{2,Float64}}

    @test GenericSphere() == GenericSphere(1.0)
    @test hash(GenericSphere()) == hash(GenericSphere(1.0))
    @test GenericSphere(2.0, 1:5) isa GenericSphere{Vector{Float64},Float64}
    @test GenericSphere(2, 1.0:5.0) isa GenericSphere{Vector{Float64},Float64}

    @test repr(UnitSphere()) == "UnitSphere()"
    @test repr(UnitSphere(Val(4))) == "UnitSphere(Val(4))"
    @test repr(UnitCircle()) == "UnitCircle()"
    @test repr(UnitCircle{BigFloat}()) == "UnitCircle{BigFloat}()"
    @test repr(VectorUnitSphere{Float64}(4)) == "UnitSphere(4)"
    @test repr(UnitSphere{Float64}()) == "UnitSphere{Float64}()"
    @test repr(Sphere(1.0,2.0)) == "Sphere(1.0, 2.0)"

    @test DomainSets.UnitCircleMap() == DomainSets.UnitCircleMap{Float64}()
    @test DomainSets.AngleMap() == DomainSets.AngleMap{Float64}()

    C = UnitCircle()
    @test SA[1.,0.] ∈ C
    @test SA[1.,1.] ∉ C
    @test approx_in(SA[1.,0.], C)
    @test !approx_in(SA[1.,1.], C)
    @test !isempty(C)
    @test isclosedset(C)
    @test !isopenset(C)
    @test parameterdomain(C) == UnitInterval()
    @test parameterdomain(UnitSphere(2)) == UnitInterval()
    @test parameterdomain(UnitSphere(3)) == UnitSphere(3)
    @test hasparameterization(C)
    @test hasparameterization(UnitSphere(2))

    p = parameterization(C)
    @test mapsize(p) == (2,)
    x = applymap(p, 1/2)
    @test jacobian(p, 0.4) ≈ SA[-2pi*sin(2pi*0.4), 2pi*cos(2pi*0.4)]
    @test diffvolume(p, 0.4) ≈ 2*pi
    @test diffvolume(p)(0.4) ≈ 2*pi
    @test approx_in(x, C)
    q = leftinverse(p)
    @test applymap(q, x) ≈ 1/2
    @test applymap(q, x) ≈ leftinverse(p, x)
    @test applymap(q, -x) ≈ 1
    @test rightinverse(q) == p
    @test rightinverse(q)(0) ≈ rightinverse(q, 0)
    @test jacobian(q) isa DomainSets.LazyJacobian
    @test jacobian(q, x) isa LinearAlgebra.Transpose{Float64,SVector{2,Float64}}

    @test boundingbox(C) == ProductDomain(ChebyshevInterval(), ChebyshevInterval())
    @test boundingbox(UnitSphere{Float64}()) == ChebyshevInterval()

    @test convert(LevelSet, UnitCircle()) isa LevelSet{SVector{2,Float64}}
    @test convert(LevelSet{SVector{2,BigFloat}}, UnitCircle()) isa LevelSet{SVector{2,BigFloat}}
    @test pseudolevel(UnitCircle(), 0.1) isa SublevelSet
    @test SA[1.05,0] ∈ pseudolevel(UnitCircle(), 0.1)
    @test SA[1.15,0] ∉ pseudolevel(UnitCircle(), 0.1)

    C2 = convert(Domain{SVector{2,BigFloat}}, C)
    @test eltype(C2) == SVector{2,BigFloat}

    @test convert(Domain{SVector{2,Float64}}, UnitSphere(2)) isa StaticUnitSphere
    @test convert(Domain{Vector{Float64}}, UnitSphere(Val(2))) isa DynamicUnitSphere

    C = 2UnitCircle() .+ SA[1.,1.]
    @test C isa DomainSets.GenericSphere
    @test approx_in(SA[3.,1.], C)

    @test canonicaldomain(C) == UnitCircle()
    @test affinematrix(mapfrom_canonical(C)) == [2 0; 0 2]
    @test affinevector(mapfrom_canonical(C)) == [1; 1]
    @test parameterdomain(C) == UnitInterval()
    @test mapfrom_parameterdomain(C) isa ComposedMap
    @test mapfrom_parameterdomain(C)(0.5) ≈ [-1; 1]
    @test mapfrom_parameterdomain(C, 0.5) ≈ [-1; 1]
    @test mapto_parameterdomain(C)([-1; 1]) ≈ 0.5
    @test mapto_parameterdomain(C, [-1; 1]) ≈ 0.5
    @test boundingbox(C) == (-1.0..3.0)^2

    C = UnitCircle() .+ SA[1,1]
    @test C isa DomainSets.GenericSphere
    @test approx_in(SA[2,1], C)

    S = UnitSphere()
    @test SA[1.,0.,0.] ∈ S
    @test SA[1.,0.,1.] ∉ S
    @test approx_in(SA[cos(1.),sin(1.),0.], S)
    @test !isempty(S)
    S2 = convert(Domain{SVector{3,BigFloat}}, S)
    @test eltype(S2) == SVector{3,BigFloat}

    @test Basis3Vector() in S

    @test issubset(UnitSphere(), UnitBall())

    S = 2 * UnitSphere() .+ SA[1.,1.,1.]
    @test S isa DomainSets.GenericSphere
    @test approx_in(SA[1. + 2*cos(1.),1. + 2*sin(1.),1.], S)
    @test !approx_in(SA[4.,1.,5.], S)

    D = UnitCircle()
    @test convert(Domain{SVector{2,BigFloat}}, D) ≡ UnitCircle{BigFloat}()
    @test SVector(1,0) in D
    @test SVector(nextfloat(1.0),0) ∉ D

    D = UnitSphere()
    @test convert(Domain{SVector{3,BigFloat}}, D) ≡ UnitSphere{SVector{3,BigFloat}}()
    @test SVector(1,0,0) in D
    @test SVector(nextfloat(1.0),0,0) ∉ D

    D4 = UnitSphere(4)
    @test D4 isa DynamicUnitSphere
    @test D4 isa VectorUnitSphere
    @test dimension(D4) == 4
    cheb = ChebyshevInterval()
    @test boundingbox(D4) == ProductDomain([cheb, cheb, cheb, cheb])

    ## sphere points
    x_sphere = [0.1,0.2,sqrt(1-(0.1)^2-(0.2)^2)]
    p1 = DomainSets.EuclideanSpherePoint(x_sphere)
    @test domain(p1) == UnitSphere(3)
    @test DomainSets.point(p1) == x_sphere
    @test p1 ∈ UnitSphere(3)

    p2 = DomainSets.EuclideanSpherePoint(SVector{3}(x_sphere))
    @test domain(p2) === UnitSphere()
    @test DomainSets.point(p2) == x_sphere
    @test p2 ∈ UnitSphere(Val(3))

    p3 = DomainSets.SphericalCoordinate(0.4, 0.5)
    @test domain(p3) === UnitSphere()
    @test p3 ∈ UnitSphere()
    @test approx_in(DomainSets.point(p3), UnitSphere())
end
