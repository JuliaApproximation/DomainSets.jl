# test_specific_domains.jl

const v = TypeFactory{SVector}()

const io = IOBuffer()
function test_specific_domains()
    @testset "$(rpad("Specific domains",80))" begin
        test_emptyspace()
        test_fullspace()
        test_point()
        test_interval()
        test_unitball()
        test_derived_unitball()
        test_cube()
        test_circle()
        test_mapped_domain()
        test_simplex()
        test_sphere()
        test_arithmetics()
        test_cartesianproduct_domain()
    end

    @testset "$(rpad("Set operations",80))" begin
      test_set_operations()
    end

end

function test_emptyspace()
    println("- an empty space")
    d1 = EmptySpace()
    show(io,d1)
    @test isempty(d1)
    @test String(take!(io)) == "the empty space with eltype Float64"
    @test eltype(d1) == Float64
    @test 0.5 ∉ d1
    @test d1 ∩ d1 == d1
    @test d1 ∪ d1 == d1
    d2 = interval()
    @test d1 ∩ d2 == d1
    @test d1 ∪ d2 == d2
    @test d2 \ d1 == d2
    @test d2 \ d2 == d1

    d2 = EmptySpace(SVector{2,Float64})
    @test isempty(d2)
    @test v[0.1,0.2] ∉ d2
end

function test_fullspace()
    println("- a full Euclidean space")
    d1 = FullSpace()
    show(io,d1)
    @test String(take!(io)) == "the full space with eltype Float64"
    @test 0.5 ∈ d1
    @test d1 ∪ d1 == d1
    @test d1 ∩ d1 == d1
    d2 = interval()
    @test d1 ∪ d2 == d1
    @test d1 ∩ d2 == d2
    @test typeof(fullspace(interval())+1) <: FullSpace
    @test typeof(fullspace(interval())*3) <: FullSpace

    d2 = FullSpace(SVector{2,Float64})
    @test v[0.1,0.2] ∈ d2

    @test d2 == Domain(SVector{2,Float64})
    @test d2 == convert(Domain,SVector{2,Float64})
    @test d2 == convert(Domain{SVector{2,Float64}}, SVector{2,Float64})
end

function test_point()
    println("- points")
    d = Domain(1.0)
    @test d isa Point
    @test 1 ∈ d
    @test 1.1 ∉ d

    @test d+1 == Domain(2.0)
    @test 1+d == Domain(2.0)
    @test 1-d == Domain(0.0)
    @test d-1 == Domain(0.0)
    @test 2d  == Domain(2.0)
    @test d*2 == Domain(2.0)
    @test d/2 == Domain(0.5)
    @test 2\d == Domain(0.5)

    d1 = Domain(Set([1,2,3]))
    d2 = Point(1) ∪ Point(2) ∪ Point(3)

    @test d1 == d2
end


function test_interval(T = Float64)
    println("- intervals")

    d = interval(zero(T), one(T))
    @test T(0.5) ∈ d
    @test T(1.1) ∉ d
    @test 0.5f0 ∈ d
    @test 1.1f0 ∉ d
    @test BigFloat(0.5) ∈ d
    @test BigFloat(1.1) ∉ d

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
    @test 1. ∈ d
    @test -1. ∉ d
    @test similar_interval(d, T(0), T(Inf)) == d



    d = negative_halfline(T)
    @test leftendpoint(d) == -T(Inf)
    @test rightendpoint(d) == zero(T)
    @test !isclosed(d)
    @test isopen(d)
    @test !iscompact(d)
    @test -1. ∈ d
    @test 1. ∉ d
    @test similar_interval(d, T(-Inf), T(0)) == d

    d = Domains.open_interval()
    @test isopen(d)
    @test !isclosed(d)
    @test leftendpoint(d)∉d
    @test rightendpoint(d)∉d
    d = Domains.closed_interval()
    @test !isopen(d)
    @test isclosed(d)
    @test leftendpoint(d) ∈ d
    @test rightendpoint(d) ∈ d
    d = HalfOpenLeftInterval()
    @test !isopen(d)
    @test !isclosed(d)
    @test leftendpoint(d) ∉ d
    @test rightendpoint(d) ∈ d
    d = HalfOpenRightInterval()
    @test !isopen(d)
    @test !isclosed(d)
    @test leftendpoint(d) ∈ d
    @test rightendpoint(d) ∉ d

    @test typeof(UnitInterval{Float64}(interval(0.,1.))) <: UnitInterval
    @test typeof(ChebyshevInterval{Float64}(interval(-1,1.))) <: ChebyshevInterval

    ## Some mappings preserve the interval structure
    # Translation
    d = interval(zero(T), one(T))
    @test d == +d

    d2 = d + one(T)
    @test typeof(d2) == typeof(d)
    @test leftendpoint(d2) == one(T)
    @test rightendpoint(d2) == 2*one(T)

    d2 = one(T) + d
    @test typeof(d2) == typeof(d)
    @test leftendpoint(d2) == one(T)
    @test rightendpoint(d2) == 2*one(T)

    d2 = d - one(T)
    @test typeof(d2) == typeof(d)
    @test leftendpoint(d2) == -one(T)
    @test rightendpoint(d2) == zero(T)

    d2 = -d
    @test typeof(d2) == typeof(d)
    @test leftendpoint(d2) == -one(T)
    @test rightendpoint(d2) == zero(T)

    d2 = one(T) - d
    @test d2 == d

    # translation for UnitInterval
    # Does a shifted unit interval return an interval?
    d = UnitInterval{T}()
    d2 = d + one(T)
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == one(T)
    @test rightendpoint(d2) == 2*one(T)

    d2 = one(T) + d
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == one(T)
    @test rightendpoint(d2) == 2*one(T)

    d2 = d - one(T)
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == -one(T)
    @test rightendpoint(d2) == zero(T)

    d2 = -d
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == -one(T)
    @test rightendpoint(d2) == zero(T)

    d2 = one(T) - d
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == zero(T)
    @test rightendpoint(d2) == one(T)


    # translation for ChebyshevInterval
    d = ChebyshevInterval{T}()
    d2 = d + one(T)
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == zero(T)
    @test rightendpoint(d2) == 2*one(T)

    d2 = one(T) + d
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == zero(T)
    @test rightendpoint(d2) == 2*one(T)

    d2 = d - one(T)
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == -2one(T)
    @test rightendpoint(d2) == zero(T)

    @test -d == d

    d2 = one(T) - d
    @test typeof(d2) <: AbstractInterval
    @test leftendpoint(d2) == zero(T)
    @test rightendpoint(d2) == 2one(T)


    # Scaling
    d = interval(zero(T), one(T))
    d3 = T(2) * d
    @test typeof(d3) == typeof(d)
    @test leftendpoint(d3) == zero(T)
    @test rightendpoint(d3) == T(2)

    d3 = d * T(2)
    @test typeof(d3) == typeof(d)
    @test leftendpoint(d3) == zero(T)
    @test rightendpoint(d3) == T(2)

    d = interval(zero(T), one(T))
    d4 = d / T(2)
    @test typeof(d4) == typeof(d)
    @test leftendpoint(d4) == zero(T)
    @test rightendpoint(d4) == T(1)/T(2)

    d4 = T(2) \ d
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

    # - setdiff of intervals
    d1 = interval(-2one(T), 2one(T))
    @test d1 \ interval(3one(T), 4one(T)) == d1
    @test d1 \ interval(zero(T), one(T)) == interval(-2one(T),zero(T)) ∪ interval(one(T), 2one(T))
    @test d1 \ interval(zero(T), 3one(T)) == interval(-2one(T),zero(T))
    @test d1 \ interval(-3one(T),zero(T)) == interval(zero(T),2one(T))
    @test d1 \ interval(-4one(T),-3one(T)) == d1
    @test d1 \ interval(-4one(T),4one(T)) == EmptySpace{T}()

    d1 \ (-3one(T)) == d1
    d1 \ (-2one(T)) == Interval{:open,:closed}(-2one(T),2one(T))
    d1 \ (2one(T)) == Interval{:closed,:open}(-2one(T),2one(T))
    d1 \ zero(T) == Interval{:closed,:open}(-2one(T),zero(T)) ∪ Interval{:open,:closed}(zero(T),2one(T))

    # - empty interval
    @test isempty(interval(one(T),zero(T)))
    @test zero(T) ∉ interval(one(T),zero(T))
    @test isempty(Interval{:open,:open}(zero(T),zero(T)))
    @test zero(T) ∉ Interval{:open,:open}(zero(T),zero(T))
    @test isempty(Interval{:open,:closed}(zero(T),zero(T)))
    @test zero(T) ∉ Interval{:open,:closed}(zero(T),zero(T))
    @test isempty(Interval{:closed,:open}(zero(T),zero(T)))
    @test zero(T) ∉ Interval{:closed,:open}(zero(T),zero(T))
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

function test_disk_ball()
    C = disk(2.0, v[1.0,1.0])
    @test in(v[2.4, 2.4], C)
    @test !in(v[3.5, 2.5], C)

    S = ball(2.0, v[1.0,1.0,1.0])
    @test v[2.9,0.0,0.0] ∈ S
    @test in(v[1,-0.9,1.0],S)
    @test in(v[1.0,1.0,-0.9],S)
    @test !in(v[2.9,2.9,1.0],S)
end
struct DerivedUnitBall<: DerivedDomain{SVector{2,Float64}}
    superdomain :: Domain

    DerivedUnitBall() = new(disk(2.0))
end
function test_derived_unitball()
    C = DerivedUnitBall()
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

    @test Domains.supereltype(C) == eltype(disk(2.0))
end

function test_circle()
    C = circle()
    v[1.,0.] ∈ C
    v[1.,1.] ∉ C
    C = circle(2., v[1.,1.])
    v[2.,1.] ∈ C
    v[2.,1.] ∉ C
    S = sphere()
    v[1.,0.,0.] ∈ S
    v[1.,0.,1.] ∉ S
    S = sphere(2., v[1.,1.,1.])
    v[1.,2.,1.] ∈ S
    v[2.,2.,1.] ∉ S

    C = disk()
    v[1.,0.] ∈ C
    v[1.,1.] ∉ C
    C = disk(2., v[1.,1.])
    v[2.,1.] ∈ C
    v[2.,1.] ∉ C
    S = ball()
    v[1.,0.,0.] ∈ S
    v[1.,0.,1.] ∉ S
    S = ball(2., v[1.,1.,1.])
    v[1.,2.,1.] ∈ S
    v[2.,2.,1.] ∉ S
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

function test_mapped_domain()
  D = cube(Val{2})
  show(io,rotate(D,1.))
  @test String(take!(io)) == "A mapped domain based on the interval [0.0, 1.0] x the interval [0.0, 1.0]"

  D = rotate(cube(Val{2}),pi)
  @test v[-0.9, -0.9] ∈ D
  @test v[-1.1, -1.1] ∉ D

  D = rotate(cube(Val{2}),pi,v[-.5,-.5])
  @test v[0.9, 0.9] ∈ D
  @test v[1.1, 1.1] ∉ D

  D = rotate(cube(Val{3})+v[-.5,-.5,-.5], pi, pi, pi)
  @test v[0.4, 0.4, 0.4] ∈ D
  @test v[0.6, 0.6, 0.6] ∉ D

  D = rotate(cube(-1.5, 2.2, 0.5, 0.7, -3.0, -1.0),pi,pi,pi,v[.35, .65, -2.])
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
    DS = D\S
    @test v[0.1, 0.1, 0.1] ∉ DS

    D1 = 2*D
    D2 = D*2
    D3 = D/2

    @test v[2., 2., 2.] ∈ D1
    @test v[0.9, 0.6,-2.5] ∉ D1
    @test v[2., 2., 2.] ∈ D2
    @test v[0.9, 0.6,-2.5] ∉ D2
    @test v[.5, .4, .45] ∈ D3
    @test v[.3, 0.6,-.2] ∉ D3

end

function test_cartesianproduct_domain()
    # ProductDomain 1
    T1 = interval(-1.0, 1.0)^2
    @test v[0.5,0.5] ∈ T1
    @test v[-1.1,0.3] ∉ T1

    T1 = cartesianproduct(interval(-1.0, 1.0), 2)
    @test v[0.5,0.5] ∈ T1
    @test v[-1.1,0.3] ∉ T1

    T2 = interval(-1.0, 1.0) × interval(-1.5, 2.5)
    @test v[0.5,0.5] ∈ T2
    @test v[-1.1,0.3] ∉ T2

    # ProductDomain 2
    T3 = ProductDomain(disk(1.05), interval(-1.0, 1.0))
    @test v[0.5,0.5,0.8] ∈ T3
    @test v[-1.1,0.3,0.1] ∉ T3

    T4 = T1×interval(-1.,1.)
    @test v[0.5,0.5,0.8] ∈ T4
    @test v[-1.1,0.3,0.1] ∉ T4

    T5 = interval(-1.,1.)×T1
    @test v[0.,0.5,0.5] ∈ T5
    @test v[0.,-1.1,0.3] ∉ T5

    T6 = T1×T1
    @test v[0.,0.,0.5,0.5] ∈ T6
    @test v[0.,0.,-1.1,0.3] ∉ T6

    io = IOBuffer()
    show(io,T1)
    @test String(take!(io)) == "the interval [-1.0, 1.0] x the interval [-1.0, 1.0]"

end


function test_set_operations()
  d1 = disk()
  d2 = interval(-.9,.9)^2
  d3 = rectangle(-.5,-.1,.5,.1)

  println("- union")
  u1 = d1 ∪ d2
  u2 = u1 ∪ d3

  u3 = d3 ∪ u1
  u4 = u1 ∪ u2
  x = SVector(0.,.15)
  y = SVector(1.1,.75)
  @test x∈u3
  @test x∈u4

  @test y∉u3
  @test y∉u4

  ũ1 = UnionDomain(d1,d2)
  @test u1 == ũ1
  ũ1 = UnionDomain((d1,d2))
  @test u1 == ũ1
  ũ2 = UnionDomain([d1,d2])
  @test ũ2 == ũ2
  @test u1 == ũ2

  # ordering doesn't matter
  @test UnionDomain(d1,d2) == UnionDomain(d2,d1)

  show(io,u1)
  @test String(take!(io)) == "a union of 2 domains:\n\t1.\t: the 2-dimensional unit ball\n\t2.\t: the interval [-0.9, 0.9] x the interval [-0.9, 0.9]\n"

  println("- intersection")
  d1 = disk()
  d2 = interval(-.4,.4)^2
  d3 = rectangle(-.5,.5,-.1,.1)
  # intersection of productdomains
  i1 = d2 & d3
  show(io,i1)
  @test String(take!(io)) == "the interval [-0.4, 0.4] x the interval [-0.1, 0.1]"
  i2 = d1 & d2
  show(io,i2)
  @test String(take!(io)) == "the intersection of 2 domains:\n\t1.\t: the 2-dimensional unit ball\n\t2.\t: the interval [-0.4, 0.4] x the interval [-0.4, 0.4]\n"

  i3 = d3 & i2
  i4 = i2 & d3
  i5 = i3 & i2

  x = SVector(0.,.05)
  y = SVector(0.,.75)
  @test x∈i3
  @test x∈i4

  @test y∉i3
  @test y∉i4

  println("- difference")
  d1 = disk()
  d2 = rectangle(-.5,.5,-.1,.1)
  # intersection of productdomains
  d = d1\d2
  show(io,d)
  @test String(take!(io)) == "the difference of 2 domains:\n\t1.\t: the 2-dimensional unit ball\n\t2.\t: the interval [-0.5, 0.5] x the interval [-0.1, 0.1]\n"

  x = SVector(0.,.74)
  y = SVector(0.,.25)
  @test x∈d
  @test x∈d

  println("- arithmetic")
  d1 = interval(0,1)
  d2 = interval(2,3)
  d = d1 ∪ d2

  @test d+1 == (d1+1) ∪ (d2+1)
  @test d-1 == (d1-1) ∪ (d2-1)
  @test 2d  == (2d1)  ∪ (2d2)
  @test d*2 == (d1*2) ∪ (d2*2)
  @test d/2 == (d1/2) ∪ (d2/2)
  @test 2\d == (2\d1) ∪ (2\d2)


  println("- different types")
  d̃1 = interval(0,1)
  d1 = interval(0f0, 1f0)
  d2 = interval(2,3)

  @test d1 ∪ d2 == d̃1 ∪ d2
end
