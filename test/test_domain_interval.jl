function test_intervals()
    T = Float64
    @testset "ClosedInterval{$T}" begin
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

        @test iscompact(d)
        @test typeof(similar_interval(d, one(T), 2*one(T))) == typeof(d)

        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)

        @test boundary(d) == uniondomain(Point(zero(T)), Point(one(T)))
        @test corners(d) == [0,1]
        @test boundingbox(d) == d

        @test similar_interval(0..1, 0, big(1.0)) isa ClosedInterval{BigFloat}

        @test canonicaldomain(d) isa ChebyshevInterval{Float64}
        @test canonicaldomain(0..1) isa ChebyshevInterval{Float64}

        @test intersectdomain(0..1, 1..3) isa Point{Int}
        @test 1..1 == Point(1)
        @test 1..1 != Point(2)
    end
    @testset "UnitInterval{$T}" begin
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
    @testset "ChebyshevInterval{$T}" begin
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
        @test String(take!(io)) == "-1.0..1.0 (Chebyshev)"
    end
    @testset "HalfLine{$T}" begin
        d = HalfLine{T}()
        @test leftendpoint(d) == zero(T)
        @test rightendpoint(d) == T(Inf)
        @test minimum(d) == infimum(d) == leftendpoint(d)
        @test supremum(d) == rightendpoint(d)
        @test_throws ArgumentError maximum(d)

        @test d ∩ d === d
        @test d ∪ d === d
        @test d \ d == EmptySpace{T}()
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
    end
    @testset "NegativeHalfLine{$T}" begin
        d = NegativeHalfLine{T}()
        @test leftendpoint(d) == -T(Inf)
        @test rightendpoint(d) == zero(T)
        @test infimum(d) == leftendpoint(d)
        @test supremum(d) == rightendpoint(d)
        @test_throws ArgumentError minimum(d)
        @test_throws ArgumentError maximum(d)

        @test d ∩ d === d
        @test d ∪ d === d
        @test d \ d == EmptySpace{T}()
        unit = UnitInterval{T}()
        cheb = ChebyshevInterval{T}()
        halfline = HalfLine{T}()
        @test unit ∩ d === EmptySpace{T}()
        @test d ∩ unit === EmptySpace{T}()
        @test d ∩ halfline === EmptySpace{T}()
        @test halfline ∩ d === EmptySpace{T}()
        @test d ∪ halfline === FullSpace{T}()
        @test halfline ∪ d === FullSpace{T}()
        @test unit \ d === unit
        @test cheb \ d === unit
        @test halfline \ d === halfline
        @test d \ unit === d
        @test d \ halfline === d

        @test !isclosedset(d)
        @test isopenset(d)
        @test !iscompact(d)
        @test -1. ∈ d
        @test 1. ∉ d
        @test approx_in(0.5, d, 1.)
        @test !approx_in(0.5, d, 0.4)
        @test similar_interval(d, T(-Inf), T(0)) == d

        @test 2d isa MappedDomain
        @test -2d isa MappedDomain

        @test boundary(d) == Point(0)
        @test leftendpoint(d) ∉ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)
    end

    @testset "OpenInterval{$T}" begin
        d = OpenInterval(0,1)
        @test isopenset(d)
        @test closure(d) == UnitInterval()

        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)
    end

    @testset "Integer intervals" begin
        d = 0..1
        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)

        d = Interval{:open,:closed}(0,1)
        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)
        @test closure(d) == 0..1

        d = Interval{:closed,:open}(0,1)
        @test leftendpoint(d) ∈ ∂(d)
        @test rightendpoint(d) ∈ ∂(d)
        @test closure(d) == 0..1
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
        @test isapprox(HalfLine(), 0..Inf, atol=100eps())
    end

    @testset "mapping between intervals" begin
        @test canonicaldomain(UnitInterval()) == UnitInterval()
        m = mapto(2..3, ChebyshevInterval())
        @test isaffine(m)
        @test m(2) ≈ -1
        @test m(3) ≈ 1
        m2 = mapto(4.0..6, 2..3)
        @test isaffine(m2)
        @test m2(4) ≈ 2
        @test m2(6) ≈ 3
    end

    @test DomainSets.isinterval(0..1)
    @test !DomainSets.isinterval(UnitBall())

    @test typeof(UnitInterval{Float64}(0.0..1.0)) <: UnitInterval
    @test typeof(ChebyshevInterval{Float64}(-1.0..1.0)) <: ChebyshevInterval

    ## Some mappings preserve the interval structure
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
    @test ChebyshevInterval() ≈ ClosedInterval(-1.0,1.0) ≈ Interval{:closed,:open}(-1.0,1.0)

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
    du5 = UnionDomain(i1) ∪ UnionDomain(i4)
    @test typeof(du5) <: UnionDomain

    # - intersection of non-overlapping intervals
    du6 = intersectdomain(i1, i4)
    @test isempty(du6)

    # - setdiff of intervals
    d1 = -2one(T).. 2one(T)
    @test d1 \ (3one(T) .. 4one(T)) == d1
    @test d1 \ (zero(T) .. one(T)) == UnionDomain((-2one(T)..zero(T))) ∪ UnionDomain((one(T).. 2one(T)))
    @test d1 \ (zero(T) .. 3one(T)) == (-2one(T) .. zero(T))
    @test d1 \ (-3one(T) .. zero(T)) == (zero(T) .. 2one(T))
    @test d1 \ (-4one(T) .. -3one(T)) == d1
    @test d1 \ (-4one(T) .. 4one(T)) == EmptySpace{T}()
    @test setdiffdomain(d1, (3one(T) .. 4one(T))) == d1
    @test setdiffdomain(d1, (zero(T) .. one(T))) == UnionDomain((-2one(T)..zero(T))) ∪ UnionDomain((one(T).. 2one(T)))
    @test setdiffdomain(d1, (zero(T) .. 3one(T))) == (-2one(T) .. zero(T))
    @test setdiffdomain(d1, (-3one(T) .. zero(T))) == (zero(T) .. 2one(T))
    @test setdiffdomain(d1, (-4one(T) .. -3one(T))) == d1
    @test setdiffdomain(d1, (-4one(T) .. 4one(T))) == EmptySpace{T}()

    # mixed types
    @test setdiffdomain(0..1, 0.0..0.5) == 0.5..1

    @test setdiffdomain(d1, -3) == d1
    @test setdiffdomain(d1, -2) == Interval{:open,:closed}(-2one(T),2one(T))
    @test setdiffdomain(d1, 2one(T)) == Interval{:closed,:open}(-2one(T),2one(T))
    @test setdiffdomain(d1, zero(T)) == UnionDomain(Interval{:closed,:open}(-2one(T),zero(T))) ∪ UnionDomain(Interval{:open,:closed}(zero(T),2one(T)))

    # - empty interval
    @test isempty(one(T)..zero(T))
    @test zero(T) ∉ (one(T)..zero(T))
    @test isempty(Interval{:open,:open}(zero(T),zero(T)))
    @test zero(T) ∉ Interval{:open,:open}(zero(T),zero(T))
    @test isempty(Interval{:open,:closed}(zero(T),zero(T)))
    @test zero(T) ∉ Interval{:open,:closed}(zero(T),zero(T))
    @test isempty(Interval{:closed,:open}(zero(T),zero(T)))
    @test zero(T) ∉ Interval{:closed,:open}(zero(T),zero(T))

    d = one(T) .. zero(T)
    @test_throws ArgumentError minimum(d)
    @test_throws ArgumentError maximum(d)
    @test_throws ArgumentError infimum(d)
    @test_throws ArgumentError supremum(d)

    # Subset relations of intervals
    @test issubset((zero(T)..one(T)), (zero(T).. 2*one(T)))
    @test issubset((zero(T)..one(T)), (zero(T).. one(T)))
    @test issubset(OpenInterval(zero(T),one(T)), zero(T) .. one(T))
    @test !issubset(zero(T) .. one(T), OpenInterval(zero(T), one(T)))
    @test issubset(UnitInterval{T}(), ChebyshevInterval{T}())

    # - convert
    d = zero(T).. one(T)
    @test d ≡ Interval(zero(T), one(T))
    @test d ≡ ClosedInterval(zero(T), one(T))

    @test convert(Domain, d) ≡ d
    @test convert(Domain{Float32}, d) ≡ (0f0 .. 1f0)
    @test convert(Domain{Float64}, d) ≡ (0.0 .. 1.0)
    @test convert(Domain, zero(T)..one(T)) ≡ d
    @test convert(Domain{T}, zero(T)..one(T)) ≡ d
    @test convert(AbstractInterval, zero(T)..one(T)) ≡ d
    @test convert(AbstractInterval{T}, zero(T)..one(T)) ≡ d
    @test convert(Interval, zero(T)..one(T)) ≡ d
    @test Interval(zero(T)..one(T)) ≡ d
    @test convert(ClosedInterval, zero(T)..one(T)) ≡ d
    @test ClosedInterval(zero(T)..one(T)) ≡ d
    @test convert(ClosedInterval{T}, zero(T)..one(T)) ≡ d
    @test ClosedInterval{T}(zero(T)..one(T)) ≡ d


    @testset "conversion from other types" begin
        @test convert(Domain{T}, 0..1) ≡ d
        @test convert(AbstractInterval{T}, 0..1) ≡ d
        @test convert(ClosedInterval{T}, 0..1) ≡ d
        @test ClosedInterval{T}(0..1) ≡ d
    end
end
