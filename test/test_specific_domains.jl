using StaticArrays, DomainSets, Test

using DomainSets:
    MappedDomain,
    similar_interval,
    GenericBall, GenericSphere

struct Basis3Vector <: StaticVector{3,Float64} end

Base.getindex(::Basis3Vector, k::Int) = k == 1 ? 1.0 : 0.0

struct NamedBall <: DomainSets.DerivedDomain{SVector{2,Float64}}
    domain  ::  Domain{SVector{2,Float64}}

    NamedBall() = new(2UnitDisk())
end

include("test_domain_interval.jl")
include("test_domain_product.jl")
include("test_domain_ball.jl")
include("test_domain_simplex.jl")

@testset "specific domains" begin
    @testset "empty space" begin
        d1 = EmptySpace()
        show(io,d1)
        @test isempty(d1)
        @test String(take!(io)) == "{} (empty domain)"
        @test eltype(d1) == Float64
        @test convert(Domain{BigFloat}, d1) === EmptySpace{BigFloat}()
        @test 0.5 ∉ d1
        @test !approx_in(0.5, d1)
        @test d1 ∩ d1 == d1
        @test d1 ∪ d1 == d1
        @test d1 \ d1 == d1
        @test boundary(d1) == d1
        @test dimension(d1) == 1
        @test isclosedset(d1)
        @test isopenset(d1)
        @test interior(d1) == d1
        @test closure(d1) == d1
        @test boundingbox(d1) == d1
        @test distance_to(d1, 0) == Inf
        @test d1 == 2..1
        @test 2..1 == d1
        d2 = 0..1
        @test d1 ∩ d2 == d1
        @test d2 ∩ d1 == d1
        @test d1 ∪ d2 == d2
        @test d2 ∪ d1 == d2
        @test d1 \ d2 == d1
        @test d2 \ d1 == d2
        @test d2 \ d2 == d1
        # Test some promotions
        @test EmptySpace{Float64}() ∪ (0..1) isa AbstractInterval{Float64}
        @test EmptySpace{Int}() ∩ (0..1.0) isa EmptySpace{Float64}
        @test EmptySpace{Int}() \ (0..1.0) isa EmptySpace{Float64}
        @test (0..1) \ EmptySpace{Float64}() isa AbstractInterval{Float64}

        d2 = EmptySpace(SVector{2,Float64})
        @test isempty(d2)
        @test SA[0.1,0.2] ∉ d2
        @test [0.1,0.2] ∉ d2
        @test 1:2 ∉ d2
        @test !approx_in(SA[0.1,0.2], d2)
        @test boundary(d2) == d2
        @test dimension(d2) == 2

        @test EmptySpace{Int}() == EmptySpace{Float64}()
        @test hash(EmptySpace{Int}()) == hash(EmptySpace{Float64}())
        @test emptyspace(0..1) == EmptySpace{Int}()
        @test emptyspace([1,2]) == EmptySpace{Int}()

        m = LinearMap(2)
        @test map_domain(m, emptyspace(Int)) == EmptySpace{Int}()
        @test mapped_domain(m, emptyspace(Int)) == EmptySpace{Int}()

        @test mapto(EmptySpace(), OpenInterval(0,0)) isa IdentityMap
        @test mapto(OpenInterval(0,0), EmptySpace()) isa IdentityMap
        @test_throws ArgumentError mapto(EmptySpace(), 0..1)
        @test_throws ArgumentError mapto(0..1, EmptySpace())
    end

    @testset "full space" begin
        d1 = FullSpace()
        @test d1 == FullSpace{Float64}()
        show(io,d1)
        @test String(take!(io)) == "{x} (full space)"
        @test convert(Domain{BigFloat}, d1) === FullSpace{BigFloat}()
        @test DomainSets.euclideanspace(Val{2}()) == FullSpace{SVector{2,Float64}}()
        @test 0.5 ∈ d1
        @test choice(d1) == 0
        @test d1 ∪ d1 == d1
        @test d1 ∩ d1 == d1
        @test isempty(d1) == false
        @test boundary(d1) == EmptySpace{Float64}()
        @test isclosedset(d1)
        @test isopenset(d1)
        @test interior(d1) == d1
        @test closure(d1) == d1
        @test dimension(d1) == 1
        @test boundingbox(d1) == d1
        @test distance_to(d1, 0) == 0
        @test DomainSets.isequaldomain1(d1, d1)
        @test DomainSets.isequaldomain2(d1, d1)
        d2 = 0..1
        @test d1 ∪ d2 == d1
        @test d1 ∩ d2 == d2
        @test d2 ∩ d1 == d2
        @test (0..1.0) \ FullSpace{Int}() isa EmptySpace{Float64}
        @test typeof(FullSpace(0..1) .+ 1) <: FullSpace
        @test typeof(FullSpace(0..1) * 3) <: FullSpace
        @test infimum(d1) == typemin(Float64)
        @test supremum(d1) == typemax(Float64)
        @test FullSpace{Int}() == FullSpace{Float64}()

        d2 = FullSpace{SVector{2,Float64}}()
        @test SA[0.1,0.2] ∈ d2
        @test approx_in(SA[0.1,0.2], d2)
        @test !isempty(d2)
        @test boundary(d2) == EmptySpace{SVector{2,Float64}}()

        @test d2 == Domain(SVector{2,Float64})
        @test d2 == convert(Domain,SVector{2,Float64})
        @test d2 == convert(Domain{SVector{2,Float64}}, SVector{2,Float64})

        @test FullSpace{Int}() == FullSpace{Float64}()
        @test hash(FullSpace{Int}()) == hash(FullSpace{Float64}())
        @test fullspace(0..1) == FullSpace{Int}()
        @test fullspace([1,2]) == FullSpace{Int}()

        @test uniondomain(UnitDisk(), FullSpace{SVector{2,Float64}}()) == FullSpace{SVector{2,Float64}}()

        # type domains
        @test typedomain(0..1) == TypeDomain{Int}()
        td1 = TypeDomain{Float64}()
        @test 2 ∈ td1
        @test 2.5 ∈ td1
        @test π ∈ td1
        @test approx_in(π, td1)
        @test π+im ∉ td1
        @test [1,2,3] ∉ td1
        @test 5 ∈ TypeDomain{Any}()
        @test TypeDomain{Int}() == TypeDomain{Int}()
    end

    @testset "number sets" begin
        # natural numbers
        @test 2 ∈ ℕ
        @test 0 ∈ ℕ
        @test -1 ∉ ℕ
        @test 1.5 ∉ ℕ
        @test π ∉ ℕ
        @test [1] ∉ ℕ
        @test 2+im ∉ ℕ
        @test 4//2 ∈ ℕ
        @test (-4)//2 ∉ ℕ
        @test 5//2 ∉ ℕ
        @test approx_in(1.01, ℕ, 0.05)
        @test !approx_in(1.01, ℕ, 0.001)
        @test !approx_in(-1.01, ℕ, 0.05)
        @test ℕ == ℕ
        # integers
        @test 2 ∈ ℤ
        @test 0 ∈ ℤ
        @test -1 ∈ ℤ
        @test 1.5 ∉ ℤ
        @test π ∉ ℤ
        @test [1] ∉ ℤ
        @test 2+im ∉ ℤ
        @test 4//2 ∈ ℤ
        @test (-4)//2 ∈ ℤ
        @test 5//2 ∉ ℤ
        @test approx_in(1.01, ℤ, 0.05)
        @test !approx_in(1.01, ℤ, 0.001)
        @test approx_in(-1.01, ℤ, 0.05)
        @test ℤ == ℤ
        # rationals
        @test 2 ∈ ℚ
        @test 0 ∈ ℚ
        @test -1 ∈ ℚ
        @test 1.5 ∉ ℚ
        @test π ∉ ℚ
        @test [1] ∉ ℚ
        @test 2+im ∉ ℚ
        @test 4//2 ∈ ℚ
        @test (-4)//2 ∈ ℚ
        @test 5//2 ∈ ℚ
        @test ℚ == ℚ
        # real numbers
        @test 2 ∈ ℝ
        @test 0 ∈ ℝ
        @test -1 ∈ ℝ
        @test 1.5 ∈ ℝ
        @test π ∈ ℝ
        @test [1] ∉ ℝ
        @test 2+im ∉ ℝ
        @test 4//2 ∈ ℝ
        @test (-4)//2 ∈ ℝ
        @test 5//2 ∈ ℝ
        @test [π] ∈ ℝ1
        @test [π,exp(1)] ∈ ℝ2
        @test [1,2,3] ∈ ℝ3
        @test [1,2,3,4] ∈ ℝ4
        @test approx_in(1.0+0.01im, ℝ, 0.05)
        @test !approx_in(1.0+0.01im, ℝ, 0.005)
        @test ℝ == ℝ
        # complex numbers
        @test 2 ∈ ℂ
        @test 0 ∈ ℂ
        @test -1 ∈ ℂ
        @test 1.5 ∈ ℂ
        @test π ∈ ℂ
        @test [1] ∉ ℂ
        @test 2+im ∈ ℂ
        @test 4//2+im ∈ ℂ
        @test 4//2 ∈ ℂ
        @test (-4)//2 ∈ ℂ
        @test 5//2 ∈ ℂ
        @test ℂ == ℂ
    end

    @testset "point" begin
        d = Domain(1.0)
        @test d isa Point
        @test 1 ∈ d
        @test 1.1 ∉ d
        @test approx_in(1.1, d, 0.2)
        @test !approx_in(1.2, d, 0.1)
        @test !isempty(d)
        @test boundary(d) == d
        @test boundingbox(d) == 1.0..1.0
        @test infimum(d) == d.x
        @test supremum(d) == d.x
        @test isclosedset(d)
        @test !isopenset(d)
        @test dimension(d) == 1
        @test isempty(interior(d))
        @test closure(d) == d
        @test canonicaldomain(d) == Point(0.0)
        @test mapfrom_canonical(d) == Translation(1.0)

        @test distance_to(d, 0.5) == abs(0.5-d.x)

        @test d .+ 1 == Domain(2.0)
        @test 1 .+ d == Domain(2.0)
        @test 1 .- d == Domain(0.0)
        @test d .- 1 == Domain(0.0)
        @test 2d  == Domain(2.0)
        @test d * 2 == Domain(2.0)
        @test d / 2 == Domain(0.5)
        @test 2 \ d == Domain(0.5)

        d1 = Domain(Set([1,2,3]))
        d2 = Point(1) ∪ Point(2) ∪ Point(3)

        @test Point(0.5) ∩ (0..1) == Point(0.5)
        @test (0..1) ∩ Point(0.5) == Point(0.5)
        @test isempty(Point(0.5) ∩ (1..2))
        @test isempty((1..2) ∩ Point(0.5))
        @test Point(0.5) ∩ DomainRef(0.5) == Point(0.5)
        @test isempty(Point(0.5) ∩ DomainRef(0))

        @test d1 == d2

        @test convert(Domain{Float64}, Point(1)) ≡ Point(1.0)
        @test Number(Point(1)) ≡ convert(Number, Point(1)) ≡ convert(Int, Point(1)) ≡ 1
        @test convert(Domain{Float64}, 1) isa Point{Float64}
        @test Int(Point(2)) == 2
        @test Float64(Point(2.0)) == 2.0
        @test Float64(Point(2)) isa Float64

        @test choice(Point(1)) == 1

        @test Point(1) + Point(2) == Point(3)
        @test Point(1) - Point(2) == Point(-1)

        @test 0.5 ∉ (0..1)\Point(0.5)
        @test (0..1) \ Point(0.5) isa  UnionDomain{Float64}
        @test (0..1) \ Point(0.0) == Interval{:open,:closed,Float64}(0,1)
        @test (0..1) \ Point(1.0) == Interval{:closed,:open,Float64}(0,1)
        @test (0..1) \ Point(2.0) == Interval{:closed,:closed,Float64}(0,1)
        @test (0..1) \ DomainRef(2.0) == (0..1) \ Point(2.0)
        @test issubset(Point(1), (0..2))
        @test Point(0.5) \ (0..1) == EmptySpace{Float64}()
        @test Point(0.5) \ (1..2) == Point(0.5)
        @test isequaldomain(Point(1), Point(1.0))
        @test isequaldomain(Point(1), 1)

        pv = Point([1,2,3])
        @test dimension(pv)==3
        @test canonicaldomain(pv) == Point([0,0,0])
        @test mapfrom_canonical(pv) == Translation(pv.x)

        @test intersectdomain(Point(1), Point(1.0)) == Point(1)
        @test isempty(intersectdomain(Point(1), Point(2.0)))

        # show
        @test repr(Point(3)) == "Point(3)"
        @test repr(Point(3.0)) == "Point(3.0)"
    end

    # @testset "intervals" begin
    #     test_intervals()
    # end

    @testset "balls" begin
        test_balls()
    end

    @testset "spheres" begin
        test_spheres()
    end

    @testset "wrapped domains" begin
        B = NamedBall()
        @test SA[1.4, 1.4] ∈ B
        @test SA[1.5, 1.5] ∉ B
        @test typeof(1.2 * B)==typeof(B * 1.2)
        @test SA[1.5,1.5] ∈ 1.2 * B
        @test SA[1.5,1.5] ∈ B * 1.2
        @test eltype(B) == eltype(2UnitDisk())

        @test hascanonicaldomain(B)
        @test canonicaldomain(B) == UnitDisk()
        @test mapfrom_canonical(B) == LinearMap{eltype(B)}(2)
        @test DomainSets.simplifies(B)
        @test canonicaldomain(DomainSets.Equal(), B) === superdomain(B)
        @test B == superdomain(B)

        @test Domain([1,2,3]) isa DomainSets.WrappedDomain{Int,Vector{Int}}
        @test Domain([1,2,3]) == Domain([1.0,2.0,3.0])
        @test canonicaldomain(Domain([1,2,3])) == [1,2,3]

        d = DomainSets.ExampleNamedDomain(UnitBall())
        @test !DomainSets.isinterval(d)
        @test superdomain(d) == UnitBall()
        @test hascanonicaldomain(d)
        @test DomainSets.simplifies(d)
        @test canonicaldomain(DomainSets.Equal(), d) === superdomain(d)
        @test d == superdomain(d)
        @test canonicaldomain(DomainSets.Isomorphic(), d) === superdomain(d)
    end

    @testset "complex unit circle/disk" begin
        C = ComplexUnitCircle()
        @test eltype(C) == Complex{Float64}
        @test isclosedset(C)
        @test !isopenset(C)
        @test 1 ∈ C
        @test 1im ∈ C
        @test 1.1im ∉ C
        @test 0.2+0.5im ∉ C
        @test 1.2+0.5im ∉ C
        @test parameterdomain(C) == UnitInterval()

        D = ComplexUnitDisk()
        @test eltype(D) == Complex{Float64}
        @test isclosedset(D)
        @test !isopenset(D)
        @test 1 ∈ D
        @test 1im ∈ D
        @test 1.1im ∉ D
        @test 0.2+0.5im ∈ D
        @test 1.2+0.5im ∉ D

        D2 = ComplexUnitDisk{BigFloat,:open}()
        @test eltype(D2) == Complex{BigFloat}
        @test isopenset(D2)
        @test 1im ∉ D2
        @test 0.999 ∈ D2

        @test repr(ComplexUnitCircle()) == "ComplexUnitCircle()"
        @test repr(ComplexUnitDisk()) == "ComplexUnitDisk()"
        @test repr(ComplexUnitDisk{Float64,:open}()) == "ComplexUnitDisk()  (open)"
        @test repr(ComplexUnitDisk{BigFloat}()) == "ComplexUnitDisk{BigFloat}()"
        @test repr(ComplexUnitDisk{BigFloat,:open}()) == "ComplexUnitDisk{BigFloat}()  (open)"

        @test pseudolevel(ComplexUnitCircle(), 0.1) isa SublevelSet{Complex{Float64},:open}
        p = pseudolevel(ComplexUnitCircle(), 0.1)
        @test 0.8 ∉ p
        @test 0.95 ∈ p
        @test 1+0.1im ∈ p
        @test 1.1+0.2im ∉ p
    end

    @testset "derived types" begin
        # Create an ellipse, the curve
        E = ellipse(2.0, 4.0)
        @test SA[2.0,0.0] ∈ E
        @test SA[0.0,4.0] ∈ E
        @test SA[2.0+1e-10,0.0] ∉ E
        @test SA[0.0,0.0] ∉ E
        E = ellipse(1, 2.0)
        @test eltype(E) == SVector{2,Float64}

        # Create an ellipse, the domain with volume
        E2 = ellipse_shape(2.0, 4.0)
        @test SA[2.0,0.0] ∈ E2
        @test SA[0.0,4.0] ∈ E2
        @test SA[2.0+1e-10,0.0] ∉ E2
        @test SA[0.0,0.0] ∈ E2
        @test SA[1.0,1.0] ∈ E2

        E2 = ellipse_shape(1, 2.0)
        @test eltype(E) == SVector{2,Float64}

        C = cylinder()
        @test eltype(C) == SVector{3,Float64}
        C2 = cylinder(1.0, 2)
        @test SA[0.5,0.2,1.5] ∈ C2
    end

    @testset "mapped_domain" begin
        @test MappedDomain(cos, 0..1.0) isa MappedDomain{Float64}
        @test MappedDomain{Float64}(cos, 0..1.0) isa MappedDomain{Float64}
        @test cos.(0..1.0) isa MappedDomain
        @test isempty(MappedDomain(LinearMap(2.0), EmptySpace()))
        @test map_domain(cos, 0..1.0) isa MappedDomain
        @test mapped_domain(acos, 0..1.0) isa MappedDomain

        # Test chaining of maps
        D = UnitCircle()
        D1 = MappedDomain(inverse(LinearMap(2)), D)
        @test typeof(D1) <: MappedDomain
        @test typeof(superdomain(D1)) <: UnitSphere
        @test isclosedset(D1)
        @test !isopenset(D1)
        @test convert(Domain{SVector{2,BigFloat}}, D1) isa MappedDomain{SVector{2,BigFloat}}
        D2 = 2 * D1
        @test typeof(superdomain(D2)) <: UnitSphere

        D = UnitInterval()^2
        show(io, textmime, rotate(D,1.))
        @test String(take!(io))[1:17] == "A .* UnitSquare()"

        D = rotate(UnitInterval()^2, π/2)
        @test SA[-0.9, 0.9] ∈ D
        @test SA[-1.1, -1.1] ∉ D
        x = choice(D)
        @test forward_map(D)(x) ≈ forward_map(D, x)
        @test DomainSets.toexternalpoint(D, x) ≈ forward_map(D, x)

        D = rotate(UnitInterval()^2, π, SA[-.5,-.5])
        @test SA[-1.5, -1.5] ∈ D
        @test SA[-0.5, -0.5] ∉ D

        D = rotate(UnitInterval()^3 .+ SA[-.5,-.5,-.5], pi, pi, pi)
        @test SA[0.4, 0.4, 0.4] ∈ D
        @test SA[0.6, 0.6, 0.6] ∉ D

        D = rotate((-1.5.. 2.2) × (0.5 .. 0.7) × (-3.0 .. -1.0), π, π, π, SA[.35, .65, -2.])
        @test SA[0.9, 0.6, -2.5] ∈ D
        @test SA[0.0, 0.6, 0.0] ∉ D

        B = mapped_domain(inverse(LinearMap(2.0)), VectorUnitBall(10))
        @test dimension(B) == 10
        @test mapped_domain(inverse_map(B), superdomain(B)) == B
        @test isopenset(interior(B))
        @test B == closure(interior(B))
        @test DomainSets.superdomain(boundary(B)) isa UnitSphere
        @test canonicaldomain(B) == VectorUnitBall(10)
        @test mapfrom_canonical(B) == forward_map(B)
        @test mapto_canonical(B) == inverse_map(B)
        @test parameterdomain(B) == canonicaldomain(B)
        @test mapfrom_parameterdomain(B) == mapfrom_canonical(B)
        @test mapto_parameterdomain(B) == mapto_canonical(B)

        # Test parametric domain
        using DomainSets: ParametricDomain
        m = AffineMap(ones(2), [4; 5])
        pd = ParametricDomain(m, UnitInterval())
        @test pd isa Domain{Vector{Float64}}
        @test forward_map(pd) == m
        @test forward_map(pd, 0.4) ≈ m(0.4)
        @test mapfrom_canonical(pd) == m
        @test canonicaldomain(pd) == 0..1
        @test boundary(pd) isa UnionDomain
        @test interior(pd) isa ParametricDomain
        @test closure(pd) isa ParametricDomain
        @test corners(pd) == [ [4.0; 5.0], [5.0; 6.0]]
        @test DomainSets.similardomain(pd, eltype(pd)) isa ParametricDomain

        @test convert(MappedDomain, 3..4.0) isa MappedDomain
        @test convert(MappedDomain{Float64}, 3..4.0) isa MappedDomain
    end

    @testset "simplex" begin
        test_simplex()
    end

    @testset "level sets" begin
        d1 = LevelSet(cos, 1.0)
        @test d1 isa LevelSet{Float64}
        @test convert(Domain{ComplexF64}, d1) isa LevelSet{ComplexF64}
        show(io,d1)
        @test String(take!(io)) == "level set f(x) = 1.0 with f = cos"
        @test 0.0 ∈ d1
        @test 0im ∈ d1
        @test 0.1 ∉ d1
        @test 0.1+1im ∉ d1

        # prod yields the function (x,y) -> x*y
        d2 = ZeroSet{SVector{2,Float64}}(prod)
        @test d2 isa ZeroSet{SVector{2,Float64}}
        @test SA[0.1,0.3] ∉ d2
        @test SA[0.0,0.3] ∈ d2
        @test SA[0.1,0.0] ∈ d2
        @test ZeroSet(cos) isa ZeroSet{Float64}
        @test convert(Domain{BigFloat}, ZeroSet(cos)) isa ZeroSet{BigFloat}
        @test convert(LevelSet, ZeroSet{BigFloat}(cos)) isa LevelSet{BigFloat}
        @test convert(LevelSet{BigFloat}, ZeroSet{Float64}(cos)) isa LevelSet{BigFloat}

        d3 = SublevelSet(cos, 0.5)
        d3_open = SublevelSet{Float64,:open}(cos,0.5)
        @test d3 isa SublevelSet{Float64,:closed}
        @test interior(d3) == d3_open
        @test closure(d3_open) == d3
        @test closure(d3) == d3
        @test interior(d3_open) == d3_open
        @test boundary(d3) == LevelSet(cos, 0.5)
        @test 3.0 ∈ d3
        @test 0.0 ∉ d3
        @test 0.0 ∉ d3_open
        show(io, d3)
        @test String(take!(io)) == "sublevel set f(x) <= 0.5 with f = cos"
        show(io, d3_open)
        @test String(take!(io)) == "sublevel set f(x) < 0.5 with f = cos"
        @test convert(Domain{BigFloat}, d3) isa SublevelSet{BigFloat,:closed}
        @test convert(Domain{BigFloat}, d3_open) isa SublevelSet{BigFloat,:open}

        d4 = SubzeroSet{SVector{2,Float64}}(prod)
        d4_open = SubzeroSet{SVector{2,Float64},:open}(prod)
        @test d4 isa SubzeroSet{SVector{2,Float64},:closed}
        @test interior(d4) == d4_open
        @test closure(d4_open) == d4
        @test closure(d4) == d4
        @test interior(d4_open) == d4_open
        @test boundary(d4) == ZeroSet{SVector{2,Float64}}(prod)
        @test SA[0.1,0.3] ∉ d4
        @test SA[-0.1,0.3] ∈ d4
        @test SA[-0.1,-0.3] ∉ d4
        @test SA[-0.1,0.3] ∈ d4_open
        convert(Domain{SVector{2,BigFloat}}, d4) isa SubzeroSet{SVector{2,BigFloat},:closed}
        convert(Domain{SVector{2,BigFloat}}, d4_open) isa SubzeroSet{SVector{2,BigFloat},:open}
        @test SubzeroSet(cos) == SubzeroSet{Float64}(cos)

        d5 = SuperlevelSet(cos, 0.5)
        d5_open = SuperlevelSet{Float64,:open}(cos, 0.5)
        @test d5 isa SuperlevelSet{Float64,:closed}
        @test interior(d5) == d5_open
        @test closure(d5_open) == d5
        @test closure(d5) == d5
        @test interior(d5_open) == d5_open
        @test boundary(d5) == LevelSet(cos, 0.5)
        @test 3.0 ∉ d5
        @test 0.0 ∈ d5
        @test 0.0 ∈ d5
        @test 0.0 ∈ d5_open
        show(io, d5)
        @test String(take!(io)) == "superlevel set f(x) >= 0.5 with f = cos"
        show(io, d5_open)
        @test String(take!(io)) == "superlevel set f(x) > 0.5 with f = cos"
        @test convert(Domain{BigFloat}, d5) isa SuperlevelSet{BigFloat}
        @test convert(Domain{BigFloat}, d5_open) isa SuperlevelSet{BigFloat,:open}

        d6 = SuperzeroSet{SVector{2,Float64}}(prod)
        d6_open = SuperzeroSet{SVector{2,Float64},:open}(prod)
        @test d6 isa SuperzeroSet{SVector{2,Float64},:closed}
        @test interior(d6) == d6_open
        @test closure(d6_open) == d6
        @test closure(d6) == d6
        @test interior(d6_open) == d6_open
        @test boundary(d6) == ZeroSet{SVector{2,Float64}}(prod)
        @test SA[0.1,0.3] ∈ d6
        @test SA[-0.1,0.3] ∉ d6
        @test SA[-0.1,-0.3] ∈ d6
        @test SuperzeroSet(cos) isa SuperzeroSet{Float64}
        @test convert(Domain{SVector{2,BigFloat}}, d6) isa SuperzeroSet{SVector{2,BigFloat},:closed}
        @test convert(Domain{SVector{2,BigFloat}}, d6_open) isa SuperzeroSet{SVector{2,BigFloat},:open}
    end

    @testset "indicator functions" begin
        using DomainSets: BoundedIndicatorFunction
        ispositive(x) = x >= 0
        d1 = IndicatorFunction(ispositive)
        @test d1 isa IndicatorFunction{Float64}
        @test DomainSets.indicatorfunction(d1) == ispositive
        show(io,d1)
        @test String(take!(io)) == "indicator domain defined by function f = ispositive"
        @test 0 ∈ d1
        @test big(0) ∈ d1
        @test -1 ∉ d1
        @test d1 ∩ ChebyshevInterval() isa BoundedIndicatorFunction
        @test ChebyshevInterval() ∩ d1 isa BoundedIndicatorFunction

        @test convert(IndicatorFunction, 0..1) isa IndicatorFunction
        @test convert(IndicatorFunction, d1) == d1
        @test convert(Domain{BigFloat}, d1) isa IndicatorFunction{BigFloat}
        @test 0.5 ∈ convert(IndicatorFunction, 0..1)

        d2 = BoundedIndicatorFunction(x -> prod(x)>0, UnitDisk())
        @test d2 isa BoundedIndicatorFunction
        @test eltype(d2) == SVector{2,Float64}
        @test SA[0.1,0.2] ∈ d2
        @test SA[-0.1,0.2] ∉ d2
        @test SA[0.1,-0.2] ∉ d2
        @test SA[-0.1,-0.2] ∈ d2
    end

    @testset "generator domains" begin
        d1 = Domain(x for x in 0..1)
        @test d1 == 0..1

        d2 = Domain(x>0 for x in -1..1)
        @test -0.5 ∉ d2
        @test 0.5 ∈ d2
        show(io, d2)
        @test String(take!(io)) == "indicator function bounded by: $(-1..1)"

        d3 = Domain(x*y>0 for (x,y) in UnitDisk())
        @test eltype(d3) == SVector{2,Float64}
        @test d3 isa DomainSets.BoundedIndicatorFunction
        @test d3.domain == UnitDisk()
        @test [0.4,0.2] ∈ d3
        @test [0.4,-0.2] ∉ d3

        d4 = Domain( x+y+z > 0 for (x,y) in UnitDisk(), z in 0..1.0)
        @test d4 isa DomainSets.BoundedIndicatorFunction{T,F,<:TupleProductDomain} where {T,F}
        @test eltype(d4) == Tuple{SVector{2, Float64}, Float64}
        @test DomainSets.indicatorfunction(d4) isa Function
        @test DomainSets.boundingdomain(d4) == TupleProductDomain(UnitDisk(), 0..1.0)
        @test boundingbox(d4) == boundingbox(DomainSets.boundingdomain(d4))
        @test ( [0.5,0.2], 0.5) ∈ d4
        @test ( [0.5,0.2], 1.5) ∉ d4
        @test ( [-0.5,-0.2], 0.1) ∉ d4

        # conditional comprehension syntax in 1d
        g5 = (x > 2 for x in 0..5.0 if x < 3)
        d5 = Domain(g5)
        @test d5 isa DomainSets.BoundedIndicatorFunction
        @test eltype(d5) == Float64
        @test -1 ∉ d5
        @test 1 ∉ d5
        @test 2 ∉ d5
        @test 2.5 ∈ d5
        @test 3 ∉ d5
        @test 3.1 ∉ d5
        @test 5.5 ∉ d5

        # conditional comprehension syntax in 2d
        g6 = (x<4 for x in 0..5.0, y in 1..2.0 if x > 3)
        d6 = Domain(g6)
        @test d6 isa DomainSets.BoundedIndicatorFunction
        @test eltype(d6) == Tuple{Float64,Float64}
        # test values for y smaller than 1, in 1..2, and larger than 2
        # test values for x smaller than 0, 3, 4, 5 and larger
        @test (-1, 1.5) ∉ d6
        @test (1, 1.5) ∉ d6
        @test (3.5, 1.5) ∈ d6
        @test (4.5, 1.5) ∉ d6
        @test (5.5, 1.5) ∉ d6
        @test (-1, 0.5) ∉ d6
        @test (1, 0.5) ∉ d6
        @test (3.5, 0.5) ∉ d6
        @test (4.5, 0.5) ∉ d6
        @test (5.5, 0.5) ∉ d6
        @test (-1, 2.5) ∉ d6
        @test (1, 2.5) ∉ d6
        @test (3.5, 2.5) ∉ d6
        @test (4.5, 2.5) ∉ d6
        @test (5.5, 2.5) ∉ d6

        g7 = (x+y+z<4 for (x,y) in UnitDisk(), z in ChebyshevInterval() if z < 0)
        d7 = Domain(g7)
        @test eltype(d7) == Tuple{SVector{2, Float64}, Float64}
        @test (SVector(0.1,0.2),-0.3) ∈ d7
        @test (SVector(0.1,0.2),0.3) ∉ d7
        @test (SVector(1,2),3) ∉ d7
    end
end

@testset "cartesian product" begin
    test_product_domains()
end

@testset "intervals" begin
    test_intervals()
end

@testset "balls" begin
    test_balls()
end

@testset "simplex" begin
    test_simplex()
end
