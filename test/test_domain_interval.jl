using DomainSets:
    PositiveRealLine,
    NonnegativeRealLine,
    NegativeRealLine,
    NonpositiveRealLine,
    isinterval

using IntervalSets

function test_intervals()
    T = Float64
    @testset "ClosedInterval" begin
        d = zero(T)..one(T)
        @test approx_in(-0.1, d, 0.2)
        @test approx_in(1.1, d, 0.2)
        @test !approx_in(-0.2, d, 0.1)
        @test !approx_in(1.2, d, 0.1)
        @test isclosedset(d)
        @test !isopenset(d)
        @test similar_interval(d, big(0), big(1)) == Interval{:closed,:closed,BigInt}(0,1)
        @test closure(d) == d
        @test isopenset(interior(d))
        @test center(d) == one(T)/2

        @test iscompact(d)
        @test isinterval(d)
        @test typeof(similar_interval(d, one(T), 2*one(T))) == typeof(d)

        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)

        @test boundary(d) == uniondomain(Point(zero(T)), Point(one(T)))
        @test corners(d) == [0,1]
        @test boundingbox(d) == d
        @test closure(d) == d

        @test similar_interval(0..1, 0, big(1.0)) isa ClosedInterval{BigFloat}
        @test convert(AbstractInterval{Float64}, 0..1) isa ClosedInterval{Float64}

        @test canonicaldomain(d) isa ChebyshevInterval{Float64}
        @test canonicaldomain(0..1) isa ChebyshevInterval{Int}

        @test intersectdomain(0..1, 1..3) isa Point{Int}
        @test 1..1 == Point(1)
        @test 1..1 != Point(2)

        @test normal(0..1, -0.3) == -1
        @test normal(0..1, 0.7) == 1
        @test distance_to(0..1, 1.5) == 0.5
        @test distance_to(0..1, 0.5) == 0
        @test distance_to(OpenInterval(0,0.0), 0.5) == Inf
    end
    @testset "UnitInterval" begin
        d = UnitInterval{T}()
        @test leftendpoint(d) == zero(T)
        @test rightendpoint(d) == one(T)
        @test minimum(d) == infimum(d) == leftendpoint(d)
        @test maximum(d) == supremum(d) == rightendpoint(d)

        @test d ∩ d === d
        @test d ∪ d === d
        @test d \ d === EmptySpace{T}()

        @test isclosedset(d)
        @test !isopenset(d)
        @test iscompact(d)

        @test convert(Domain, d) ≡ d
        @test convert(Domain{T}, d) ≡ d
        @test convert(AbstractInterval, d) ≡ d
        @test convert(AbstractInterval{T}, d) ≡ d
        @test convert(UnitInterval, d) ≡ d
        @test convert(UnitInterval{T}, d) ≡ d
        @test convert(Domain{Float64}, d) ≡ UnitInterval()
        @test convert(AbstractInterval{Float64}, d) ≡ UnitInterval()
        @test convert(UnitInterval{Float64}, d) ≡ UnitInterval()

        @test canonicaldomain(d) === d
    end
    @testset "ChebyshevInterval" begin
        d = ChebyshevInterval{T}()
        @test leftendpoint(d) == -one(T)
        @test rightendpoint(d) == one(T)
        @test minimum(d) == infimum(d) == leftendpoint(d)
        @test maximum(d) == supremum(d) == rightendpoint(d)

        # test for promotion
        @test intersect(0..1, ChebyshevInterval()) == 0..1

        @test d ∩ d === d
        @test d ∪ d === d
        @test d \ d === EmptySpace{T}()
        unit = UnitInterval{T}()
        @test d ∩ unit === unit
        @test unit ∩ d === unit
        @test d ∪ unit === d
        @test unit ∪ d === d
        @test unit \ d === EmptySpace{T}()
        @test -d == d

        @test isclosedset(d)
        @test !isopenset(d)
        @test iscompact(d)

        @test convert(Domain, d) ≡ d
        @test convert(Domain{T}, d) ≡ d
        @test convert(AbstractInterval, d) ≡ d
        @test convert(AbstractInterval{T}, d) ≡ d
        @test convert(ChebyshevInterval, d) ≡ d
        @test convert(ChebyshevInterval{T}, d) ≡ d
        @test convert(Domain{Float64}, d) ≡ ChebyshevInterval()
        @test convert(AbstractInterval{Float64}, d) ≡ ChebyshevInterval()
        @test convert(ChebyshevInterval{Float64}, d) ≡ ChebyshevInterval()

        @test canonicaldomain(d) === d

        show(io, ChebyshevInterval())
        @test String(take!(io)) == "$(-1.0..1.0) (Chebyshev)"
    end
    @testset "HalfLine" begin
        d = HalfLine{T}()
        d_open = PositiveRealLine{T}()
        @test d_open === HalfLine{T,:open}()
        @test leftendpoint(d) == zero(T)
        @test rightendpoint(d) == T(Inf)
        @test DomainSets.endpoints(d) == DomainSets.endpoints(d_open)
        @test minimum(d) == infimum(d) == leftendpoint(d)
        @test infimum(d_open) == leftendpoint(d_open)
        @test_throws ArgumentError minimum(d_open)
        @test supremum(d) == rightendpoint(d)
        @test supremum(d_open) == rightendpoint(d_open)
        @test_throws ArgumentError maximum(d)
        @test_throws ArgumentError maximum(d_open)

        @test d ∩ d === d
        @test d ∪ d === d
        @test d \ d == EmptySpace{T}()
        @test d_open ∩ d_open === d_open
        @test d_open ∪ d_open === d_open
        @test d_open \ d_open == EmptySpace{T}()
        @test d ∩ d_open === d_open
        @test d_open ∩ d === d_open
        @test d ∪ d_open === d
        @test d_open ∪ d === d
        unit = UnitInterval{T}()
        cheb = ChebyshevInterval{T}()
        @test d ∩ unit === unit
        @test unit ∩ d === unit
        @test d ∩ cheb === unit
        @test cheb ∩ d === unit
        @test d ∪ unit === d
        @test unit ∪ d === d
        @test unit \ d === EmptySpace{T}()

        @test !isclosedset(d)
        @test !isopenset(d)
        @test !iscompact(d)
        @test 1. ∈ d
        @test -1. ∉ d
        @test approx_in(-0.1, d, 0.5)
        @test !approx_in(-0.5, d, 0.1)
        @test similar_interval(d, T(0), T(Inf)) == d

        @test 2d isa MappedDomain
        @test -2d isa MappedDomain

        @test boundary(d) == Point(0)
        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∉ ∂(d)
        @test choice(d) == 0
        @test choice(d_open) == 1

        @test isopenset(d_open)
        @test !isclosedset(d_open)
        @test boundary(d_open) == boundary(d)
        @test 1 ∈ d_open
        @test 0 ∉ d_open
        @test interior(d) == d_open
        @test closure(d_open) == d

        show(io, HalfLine())
        @test String(take!(io)) == "0.0 .. Inf (closed-open) (HalfLine)"
    end
    @testset "NegativeHalfLine" begin
        d = NegativeHalfLine{T}()
        d_closed = NonpositiveRealLine{T}()
        @test d_closed == NegativeHalfLine{T,:closed}()
        @test leftendpoint(d) == -T(Inf)
        @test rightendpoint(d) == zero(T)
        @test infimum(d) == leftendpoint(d)
        @test supremum(d) == rightendpoint(d)
        @test_throws ArgumentError minimum(d)
        @test_throws ArgumentError maximum(d)
        @test_throws ArgumentError minimum(d_closed)
        @test maximum(d_closed) == supremum(d_closed) == rightendpoint(d_closed)

        @test d ∩ d === d
        @test d ∪ d === d
        @test d \ d == EmptySpace{T}()
        halfline = HalfLine{T}()
        halfline_open = HalfLine{T,:open}()
        @test isempty(d ∩ halfline)
        @test d ∩ halfline === EmptySpace{T}()
        @test d ∩ halfline_open === EmptySpace{T}()
        @test d_closed ∩ halfline == Point(zero(T))
        @test halfline ∩ d === EmptySpace{T}()
        @test d ∪ halfline === RealLine{T}()
        @test halfline ∪ d === RealLine{T}()
        @test d_closed ∪ halfline === RealLine{T}()
        @test halfline ∪ d_closed === RealLine{T}()
        @test halfline_open ∪ d_closed === RealLine{T}()
        @test d_closed ∪ halfline_open === RealLine{T}()
        @test setdiffdomain(halfline, d) === halfline
        @test setdiffdomain(halfline, d_closed) === halfline_open
        unit = UnitInterval{T}()
        cheb = ChebyshevInterval{T}()
        @test unit ∩ d === EmptySpace{T}()
        @test d ∩ unit === EmptySpace{T}()
        @test unit \ d === unit
        @test cheb \ d === unit
        @test halfline \ d === halfline
        @test d \ unit === d
        @test d \ halfline === d

        @test !isclosedset(d)
        @test isopenset(d)
        @test !iscompact(d)
        @test -1. ∈ d
        @test 0 ∉ d
        @test 1. ∉ d
        @test approx_in(0.5, d, 1.)
        @test !approx_in(0.5, d, 0.4)
        @test similar_interval(d, T(-Inf), T(0)) == d

        @test 2d isa MappedDomain
        @test -2d isa MappedDomain

        @test boundary(d) == Point(0)
        @test leftendpoint(d) ∉ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)
        @test choice(d) == -1
        @test choice(d_closed) == 0

        # additional tests for the right-closed negative half line
        @test !isrightopen(d_closed)
        @test !isopenset(d_closed)
        @test boundary(d_closed) == boundary(d)
        @test -1 ∈ d_closed
        @test 0 ∈ d_closed
        @test interior(d_closed) == d
        @test closure(d) == d_closed

        show(io, NegativeHalfLine())
        @test String(take!(io)) == "-Inf .. 0.0 (open) (NegativeHalfLine)"
    end

    @testset "RealLine" begin
        @test RealLine() isa RealLine{Float64}
        d = RealLine{T}()
        @test isopenset(d)
        @test DomainSets.isfullspace(d)
        @test 1 ∈ d
        @test choice(d) ∈ d
        @test choice(d) isa T
        @test isempty(boundary(d))
        @test interior(d) == d
        @test_throws ArgumentError minimum(d)
        @test_throws ArgumentError maximum(d)
        @test interior(d) == d
        @test DomainSets.similardomain(d, widen(T)) isa RealLine{widen(T)}
        @test setdiffdomain(d, NonnegativeRealLine{T}()) === NegativeRealLine{T}()
        @test setdiffdomain(d, PositiveRealLine{T}()) === NonpositiveRealLine{T}()
        @test setdiffdomain(d, NonpositiveRealLine{T}()) === PositiveRealLine{T}()
        @test setdiffdomain(d, NegativeRealLine{T}()) === NonnegativeRealLine{T}()

        @test similar_interval(d, T(-Inf), T(Inf)) == d

        @test (0..1) ∩ d == 0..1
        @test d ∩ (0..1) == 0..1
        @test (0..1) ∪ d == d
        @test d ∪ (0..1) == d
    end

    @testset "OpenInterval" begin
        d = OpenInterval(0.0,1.0)
        @test isopenset(d)
        @test closure(d) == UnitInterval()

        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)

        @test canonicaldomain(d) == OpenInterval(-1,1)
        @test canonicaldomain(OpenInterval(0,1)) isa OpenInterval{Int}
        @test canonicaldomain(OpenInterval(0,0)) == EmptySpace()
        @test canonicaldomain(OpenInterval(2,Inf)) == DomainSets.PositiveRealLine()
        @test mapfrom_canonical(OpenInterval(2,Inf)) isa AffineMap
        @test canonicaldomain(OpenInterval(-Inf,Inf)) == DomainSets.RealLine()
        @test canonicaldomain(OpenInterval(Inf,-Inf)) == DomainSets.RealLine()
        @test canonicaldomain(OpenInterval(Inf,Inf)) === EmptySpace{Float64}()

        @test isempty(boundary(OpenInterval(0,0)))
        @test length(corners(OpenInterval(0,0))) == 0
    end
    @testset "Halfopen intervals" begin
        d1 = Interval{:closed,:open}(0.0,1.0)
        @test closure(d1) == 0..1
        @test leftendpoint(d1) ∈ ∂(d1)
        @test rightendpoint(d1) ∈ ∂(d1)
        d2 = Interval{:open,:closed}(0,1)
        @test closure(d2) == 0..1
        @test leftendpoint(d2) ∈ ∂(d2)
        @test rightendpoint(d2) ∈ ∂(d2)

        @test canonicaldomain(Interval{:closed,:open}(0,1)) isa Interval{:closed,:open,Int}
        @test canonicaldomain(Interval{:open,:closed}(0,1)) isa Interval{:open,:closed,Int}

        @test canonicaldomain(Interval{:closed,:open}(0,1)) == Interval{:closed,:open}(-1,1)
        @test canonicaldomain(Interval{:open,:closed}(0,1)) == Interval{:open,:closed}(-1,1)
        @test canonicaldomain(Interval{:closed,:open}(0,0)) == EmptySpace()
        @test canonicaldomain(Interval{:open,:closed}(0,0)) == EmptySpace()
        @test canonicaldomain(Interval{:closed,:open}(0,Inf)) == HalfLine()
        @test canonicaldomain(Interval{:open,:closed}(-Inf,0)) == HalfLine()
        @test_throws ArgumentError canonicaldomain(Interval{:closed,:open}(Inf,0))
        @test_throws ArgumentError canonicaldomain(Interval{:open,:closed}(0,Inf))
    end

    @testset "Integer intervals" begin
        d1 = 0..1
        @test leftendpoint(d1) ∈ ∂(d1)
        @test rightendpoint(d1) ∈ ∂(d1)

        d2 = Interval{:open,:closed}(0,1)
        @test leftendpoint(d2) ∈ ∂(d2)
        @test rightendpoint(d2) ∈ ∂(d2)
        @test closure(d2) == 0..1
        @test choice(d2) == 1

        d3 = Interval{:closed,:open}(0,1)
        @test leftendpoint(d3) ∈ ∂(d3)
        @test rightendpoint(d3) ∈ ∂(d3)
        @test closure(d3) == 0..1
        @test choice(d3) == 0

        d4 = Interval{:open,:open}(0,1)
        @test_throws BoundsError choice(d4)
    end

    @testset "approximate in for open and closed intervals" begin
        @test !approx_in(1.0, Interval{:open,:open}(0,1), 0)
        @test !approx_in(1.0, Interval{:closed,:open}(0,1), 0)
        @test !approx_in(0.0, Interval{:open,:open}(0,1), 0)
        @test !approx_in(0.0, Interval{:open,:closed}(0,1), 0)
        @test approx_in(0.0, Interval{:closed,:closed}(0,1), 0)
        @test approx_in(1.0, Interval{:closed,:closed}(0,1), 0)

        @test (0..1) ≈ (1e-16..1)
        @test (-1..0) ≈ (-1..1e-16)
        @test_skip isapprox(HalfLine(), 0..Inf)
    end

    @testset "mapping between intervals" begin
        @test canonicaldomain(UnitInterval()) == UnitInterval()
        m = mapto(2..3, ChebyshevInterval())
        @test isaffinemap(m)
        @test m(2) ≈ -1
        @test m(3) ≈ 1
        m2 = mapto(4.0..6, 2..3)
        @test isaffinemap(m2)
        @test m2(4) ≈ 2
        @test m2(6) ≈ 3
    end

    @test DomainSets.isinterval(0..1)
    @test !DomainSets.isinterval(UnitBall())

    @test typeof(UnitInterval{Float64}(0.0..1.0)) <: UnitInterval
    @test typeof(ChebyshevInterval{Float64}(-1.0..1.0)) <: ChebyshevInterval

    @testset "some mappings that preserve interval structure" begin
        # Translation
        d = zero(T)..one(T)

        @test -Interval{:closed,:open}(2,3) isa Interval{:open,:closed}

        d2 = d .+ one(T)
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = one(T) .+ d
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = d .- one(T)
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = -d
        @test typeof(d2) == typeof(d)
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = one(T) .- d
        @test d2 == d

        # translation for UnitInterval
        # Does a shifted unit interval return an interval?
        d = UnitInterval{T}()
        d2 = d .+ one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = one(T) .+ d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == one(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = d .- one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = -d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == -one(T)
        @test rightendpoint(d2) == zero(T)

        d2 = one(T) .- d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == one(T)


        # translation for ChebyshevInterval
        d = ChebyshevInterval{T}()
        d2 = d .+ one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = one(T) .+ d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == 2*one(T)

        d2 = d .- one(T)
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == -2one(T)
        @test rightendpoint(d2) == zero(T)

        @test -d == d

        d2 = one(T) .- d
        @test typeof(d2) <: AbstractInterval
        @test leftendpoint(d2) == zero(T)
        @test rightendpoint(d2) == 2one(T)
    end

    # Scaling
    d = zero(T)..one(T)
    d3 = T(2) * d
    @test typeof(d3) == typeof(d)
    @test leftendpoint(d3) == zero(T)
    @test rightendpoint(d3) == T(2)

    d3 = d * T(2)
    @test typeof(d3) == typeof(d)
    @test leftendpoint(d3) == zero(T)
    @test rightendpoint(d3) == T(2)

    d = zero(T)..one(T)
    d4 = d / T(2)
    @test typeof(d4) == typeof(d)
    @test leftendpoint(d4) == zero(T)
    @test rightendpoint(d4) == T(1)/T(2)

    d4 = T(2) \ d
    @test typeof(d4) == typeof(d)
    @test leftendpoint(d4) == zero(T)
    @test rightendpoint(d4) == T(1)/T(2)

    # Equality
    @test ChebyshevInterval() == ClosedInterval(-1.0,1.0) == ClosedInterval(-1,1)
    @test ChebyshevInterval() ≠ Interval{:closed,:open}(-1.0,1.0)
    @test ChebyshevInterval() ≈ ClosedInterval(-1.0,1.0)

    @testset "union and intersection of intervals" begin
        # Union and intersection of intervals
        i1 = zero(T)..one(T)
        i2 = one(T)/3 .. one(T)/2
        i3 = one(T)/2 .. 2*one(T)
        i4 = T(2) .. T(3)
        # - union of completely overlapping intervals
        du1 = uniondomain(i1, i2)
        @test typeof(du1) <: AbstractInterval
        @test leftendpoint(du1) == leftendpoint(i1)
        @test rightendpoint(du1) == rightendpoint(i1)
        @test uniondomain(0..1, 0.1..1.5) isa AbstractInterval{Float64}

        # - intersection of completely overlapping intervals
        du2 = intersectdomain(i1, i2)
        @test typeof(du2) <: AbstractInterval
        @test leftendpoint(du2) == leftendpoint(i2)
        @test rightendpoint(du2) == rightendpoint(i2)

        # - union of partially overlapping intervals
        du3 = uniondomain(i1, i3)
        @test typeof(du3) <: AbstractInterval
        @test leftendpoint(du3) == leftendpoint(i1)
        @test rightendpoint(du3) == rightendpoint(i3)

        @test uniondomain(OpenInterval(0,1), 0..2) == 0..2
        @test uniondomain(OpenInterval(0,1), OpenInterval(0,2)) == OpenInterval(0,2)
        @test uniondomain(1..2, 0..1.5) == 0..2.0
        @test uniondomain(1..2.5, 0.8..1.5) == 0.8..2.5
        @test uniondomain(1..2.5, 0.8..2.5) == 0.8..2.5
        @test uniondomain(OpenInterval(1,2.5), OpenInterval(0.8,2.5)) == OpenInterval(0.8,2.5)

        # - intersection of partially overlapping intervals
        du4 = intersectdomain(i1, i3)
        @test typeof(du4) <: AbstractInterval
        @test leftendpoint(du4) == leftendpoint(i3)
        @test rightendpoint(du4) == rightendpoint(i1)

        # - union of non-overlapping intervals
        du5 = UnionDomain((i1,)) ∪ UnionDomain((i4,))
        @test typeof(du5) <: UnionDomain

        # - intersection of non-overlapping intervals
        du6 = intersectdomain(i1, i4)
        @test isempty(du6)

        # - setdiff of intervals
        @test setdiffdomain(1.0..0.0, 2.0..3.0) == 1..0.0
        @test setdiffdomain(1.0..0.0, -3.0 .. -2.0) == 1..0.0
        @test setdiffdomain(1.0..1.0, 2.0..3) == 1.0..1.0
        @test setdiffdomain(1.0..1.0, 0.0..3) isa EmptySpace{Float64}
        @test setdiffdomain(1.0..2.0, 0.0..0.0) == 1.0..2.0
        @test setdiffdomain(1.0..2.0, 1.5..1.5) ==
            uniondomain(Interval{:closed,:open}(1,1.5), Interval{:open,:closed}(1.5,2))
        @test setdiffdomain(1.0..2.0, 0.0..0.5) == 1.0..2.0
        @test setdiffdomain(1.0..2.0, 2.5..3.0) == 1.0..2.0
        @test setdiffdomain(1.0..2.0, 0.0..3.0) == EmptySpace{Float64}()
        @test setdiffdomain(1.0..2.0, 2.0..3.0) == Interval{:closed,:open}(1.0..2.0)
        @test setdiffdomain(1.0..2.0, OpenInterval(2.0,3.0)) == Interval{:closed,:closed}(1.0..2.0)
        @test setdiffdomain(1.0..2.0, 1.5..3.0) == Interval{:closed,:open}(1,1.5)
        @test setdiffdomain(1.0..2.0, OpenInterval(1.5,3.0)) == 1.0..1.5
        @test setdiffdomain(1.0..2.0, 1.5..2.0) == Interval{:closed,:open}(1.0,1.5)
        @test setdiffdomain(1.0..2.0, OpenInterval(1.5,2.0)) == (1.0..1.5) ∪ Point(2.0)
        @test setdiffdomain(1.0..2.0, Interval{:closed,:open}(1.5,2.0)) == Interval{:closed,:open}(1.0..1.5) ∪ Point(2.0)
        @test setdiffdomain(0.0..3.0, 1.0..2.0) ==
            uniondomain(Interval{:closed,:open}(0.0,1.0), Interval{:open,:closed}(2.0,3.0))
        @test setdiffdomain(0.0..1.0, 0.0..2.0) isa EmptySpace{Float64}
        @test setdiffdomain(0.0..1.0, OpenInterval(0.0,2.0)) == Point(0.0)
        @test setdiffdomain(0.0..2.0, 0.0..1.0) == Interval{:open,:closed}(1.0, 2.0)
        @test setdiffdomain(0.0..2.0, OpenInterval(0.0,1.0)) ==
            uniondomain(Point(0.0), Interval{:closed,:closed}(1.0,2.0))
        @test setdiffdomain(0.0..1.0, 0.0..1.0) isa EmptySpace{Float64}
        @test setdiffdomain(0.0..1.0, OpenInterval(0.0,1.0)) == Point(0.0) ∪ Point(1.0)
        @test setdiffdomain(0.0..1.0, Interval{:closed,:open}(0.0,1.0)) == Point(1.0)
        @test setdiffdomain(0.0..1.0, Interval{:open,:closed}(0.0,1.0)) == Point(0.0)
        @test setdiffdomain(1.0..2.0, 0.0..1.0) == Interval{:open,:closed}(1.0,2.0)
        @test setdiffdomain(1.0..2.0, OpenInterval(0.0,1.0)) == 1.0..2.0
        @test setdiffdomain(1.0..2.0, 0.0..1.5) == Interval{:open,:closed}(1.5,2.0)
        @test setdiffdomain(1.0..2.0, OpenInterval(0.0,1.5)) == 1.5..2.0
        @test setdiffdomain(1.0..2.0, 0.0..2.0) == EmptySpace{Float64}()
        @test setdiffdomain(1.0..2.0, OpenInterval(0.0,2.0)) == Point(2.0)
    end

    @testset "conversion to domain" begin
        d = zero(T).. one(T)
        @test convert(Domain, d) ≡ d
        @test convert(Domain{Float32}, d) ≡ (0f0 .. 1f0)
        @test convert(Domain{Float64}, d) ≡ (0.0 .. 1.0)
        @test convert(Domain, zero(T)..one(T)) ≡ d
        @test convert(Domain{T}, zero(T)..one(T)) ≡ d
    end

    @testset "operations on special intervals" begin
        @test isempty(UnitInterval() \ UnitInterval())
        @test ChebyshevInterval() ∩ UnitInterval() === UnitInterval()
        @test UnitInterval() ∩ ChebyshevInterval() === UnitInterval()
        @test UnitInterval() ∩ NonnegativeRealLine() === UnitInterval()
        @test NonnegativeRealLine() ∩ UnitInterval() === UnitInterval()
        @test isempty(UnitInterval() ∩ NegativeRealLine())
        @test isempty(NegativeRealLine() ∩ UnitInterval())
        @test ChebyshevInterval() ∩ NonnegativeRealLine() === UnitInterval()
        @test NonnegativeRealLine() ∩ ChebyshevInterval() === UnitInterval()
        @test isempty(HalfLine() ∩ NegativeRealLine())
        @test isempty(NegativeRealLine() ∩ HalfLine())
        @test NonpositiveRealLine() ∩ NonnegativeRealLine() == Point(0)
        @test NonnegativeRealLine() ∩ NonpositiveRealLine() == Point(0)
        @test NonnegativeRealLine() ∩ PositiveRealLine() === PositiveRealLine()
        @test PositiveRealLine() ∩ NonnegativeRealLine() === PositiveRealLine()
        @test (0..1) ∩ RealLine() == (0..1)
        @test OpenInterval(0,1) ∩ RealLine() == OpenInterval(0,1)
        @test RealLine() ∩ (0..1) == (0..1)
        @test RealLine() ∩ OpenInterval(0,1) == OpenInterval(0,1)
        @test RealLine() ∩ RealLine() == RealLine()

        @test UnitInterval() ∪ ChebyshevInterval() === ChebyshevInterval()
        @test ChebyshevInterval() ∪ UnitInterval() === ChebyshevInterval()
        @test UnitInterval() ∪ NonnegativeRealLine() === NonnegativeRealLine()
        @test NonnegativeRealLine() ∪ UnitInterval() === NonnegativeRealLine()
        @test NonnegativeRealLine() ∪ NonpositiveRealLine() === RealLine()
        @test NonpositiveRealLine() ∪ NonnegativeRealLine() === RealLine()
        @test PositiveRealLine() ∪ NonpositiveRealLine() === RealLine()
        @test NonpositiveRealLine() ∪ PositiveRealLine() === RealLine()
        @test NegativeRealLine() ∪ NonnegativeRealLine() === RealLine()
        @test NonnegativeRealLine() ∪ NegativeRealLine() === RealLine()
        @test NonnegativeRealLine() ∪ PositiveRealLine() === NonnegativeRealLine()
        @test PositiveRealLine() ∪ NonnegativeRealLine() === NonnegativeRealLine()
        @test NonpositiveRealLine() ∪ NegativeHalfLine() === NonpositiveRealLine()
        @test NegativeHalfLine() ∪ NonpositiveRealLine() === NonpositiveRealLine()
        @test NonpositiveRealLine() ∪ NonpositiveRealLine() === NonpositiveRealLine()
        @test RealLine() ∪ RealLine() === RealLine()
        @test uniondomain(0..1, RealLine()) === RealLine()
        @test uniondomain(OpenInterval(0,1), RealLine()) === RealLine()
        @test uniondomain(RealLine(), 0..1) === RealLine()
        @test uniondomain(RealLine(), OpenInterval(0,1)) === RealLine()

        @test isempty(UnitInterval() \ ChebyshevInterval())
        @test isempty(UnitInterval() \ NonnegativeRealLine())
        @test UnitInterval() \ NegativeRealLine() === UnitInterval()
        @test ChebyshevInterval() \ NegativeRealLine() === UnitInterval()
        @test HalfLine() \ NegativeRealLine() === HalfLine()
        @test HalfLine() \ NonpositiveRealLine() === HalfLine{T,:open}()
        @test NegativeRealLine() \ UnitInterval() === NegativeRealLine()
        @test NegativeHalfLine() \ PositiveRealLine() === NegativeHalfLine()
        @test NegativeHalfLine() \ NonnegativeRealLine() === NegativeRealLine()
        @test RealLine() \ NonnegativeRealLine() === NegativeRealLine()
        @test RealLine() \ PositiveRealLine() === NonpositiveRealLine()
        @test RealLine() \ NonpositiveRealLine() === PositiveRealLine()
        @test RealLine() \ NegativeRealLine() === NonnegativeRealLine()

        @test issubset(UnitInterval{T}(), ChebyshevInterval{T}())
    end
end
