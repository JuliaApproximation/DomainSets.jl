@testset "set operations" begin
    @testset "union" begin
        d1 = UnitDisk()
        d2 = (-.9..0.9)^2
        d3 = ProductDomain(-.5 .. -.1, ChebyshevInterval())
        d4 = (0.0..1.5)
        d5 = [1.0,3.0]

        @test convert(Domain{SVector{2,Float64}}, d3) isa Domain{SVector{2,Float64}}

        u1 = UnitDisk() ∪ (-.9..0.9)^2
        u2 = u1 ∪ d3
        @test dimension(u1) == 2
        @test dimension(u2) == 2
        @test boundingbox(u1) == ChebyshevInterval()^2

        u3 = d3 ∪ u1
        u4 = u1 ∪ u2
        @test SVector(0.,.15) ∈ u3
        @test SVector(0.,.15) ∈ u4
        @test SVector(1.1,.75) ∉ u3
        @test SVector(1.1,.75) ∉ u4

        ũ1 = UnionDomain(d1,d2)
        @test u1 == ũ1
        @test UnionDomain{SVector{2,Float64}}(d1,d2) == ũ1
        ũ1 = UnionDomain((d1,d2))
        @test u1 == ũ1
        ũ2 = UnionDomain([d1,d2])
        @test ũ2 == ũ2
        @test u1 == ũ2
        @test UnionDomain{SVector{2,Float64}}(d1) isa UnionDomain

        # Don't create a union with two identical elements
        @test UnitDisk() ∪ UnitDisk() isa UnitDisk

        # union with non-Domain type that implements domain interface
        u45 = (0.0..1.5) ∪ [1.0,3.0]
        @test u45 isa Domain{Float64}
        @test u45 isa UnionDomain
        @test eltype(component(u45,1)) == Float64
        @test eltype(component(u45,2)) == Float64
        @test 0.2 ∈ u45
        @test 1.2 ∈ u45
        @test 3 ∈ u45
        @test -1.2 ∉ u45
        @test convert(Domain{BigFloat}, u45) isa Domain{BigFloat}

        u45b = (0.0..1.5) ∪ [1,3]
        @test u45b isa Domain{Float64}
        @test component(u45b,2) isa AbstractArray{Float64}
        @test [1,3] ∪ (0.0..1.5) isa Domain{Float64}

        @test issubset([0,1], 0..1)
        @test !issubset([0,1,2], 0..1)
        @test issubset(Set([0,1]), 0..1)
        @test !issubset(Set([0,2]), 0..1)

        @test uniondomain() == EmptySpace{Any}()
        @test uniondomain(0..1) == 0..1
        @test uniondomain(0..1, [0,1]) == 0..1
        @test uniondomain([0,1], 0..1) == 0..1

        @test uniondomain([0,1], [0.0,1.0]) == [0,1]
        @test uniondomain([0,2], [0.0,1.0]) isa UnionDomain

        # larger union expressions
        @test uniondomain(0..1, 1..3, Point(0.4), 2..5, FullSpace(), Point(-0.2)) isa FullSpace
        @test uniondomain(0..1, 1..3, Point(0.4)) == 0..3
        @test uniondomain(0..1, Point(0.4), 1..3) == 0..3
        @test uniondomain(Point(0.4), 0..1, 1..3) == 0..3

        # ordering doesn't matter
        @test UnionDomain(d1,d2) == UnionDomain(d2,d1)

        @test UnionDomain((d1,d2)) == UnionDomain(d1,d2)
        @test UnionDomain(d1) isa UnionDomain
        @test UnionDomain(UnionDomain(d1,d2),d3) == UnionDomain(d3,UnionDomain(d1,d2))

        @test convert(Domain, [0..1,2..3,3..4]) isa UnionDomain{Int}
        @test convert(Domain{Float64}, [0..1,2..3,3..4]) isa UnionDomain{Float64}
        @test convert(Domain, Set([0..1,2..3,3..4])) isa UnionDomain{Int}
        @test convert(Domain{Float64}, Set([0..1,2..3,3..4])) isa UnionDomain{Float64}
        @test convert(Domain, Set([1,2,3])) isa UnionDomain{Int}
        @test 2 ∈ convert(Domain, Set([1,2,3]))
        @test 2 ∈ convert(Domain{Float64}, Set([1,2,3]))
        @test 4 ∉ convert(Domain, Set([1,2,3]))

        @test interior(uniondomain(0..1, 2..3)) == uniondomain(OpenInterval(0,1),OpenInterval(2,3))
        @test closure(uniondomain(OpenInterval(0,1),OpenInterval(2,3))) == uniondomain(0..1, 2..3)

        @test uniondomain(0..1, 2..3) \ (0.5..2.5) == uniondomain(Interval{:closed,:open}(0..0.5), Interval{:open,:closed}(2.5..3.0))
        @test uniondomain(0..1, 2..3) \ uniondomain(0.5..2.5, -2..(-1)) == uniondomain(Interval{:closed,:open}(0..0.5), Interval{:open,:closed}(2.5,3.0))

        @test !isempty(u1)
        show(io, textmime, u1)
        @test String(take!(io)) == "UnitDisk() ∪ ((-0.9..0.9) × (-0.9..0.9))"

        # repeated union
        @test ncomponents(uniondomain(UnitBall{Float64}(), UnitInterval(), UnitInterval())) == 2
        @test ncomponents(uniondomain(UnitInterval(), UnitBall{Float64}(), UnitInterval())) == 2
        @test ncomponents(uniondomain(UnitInterval(), UnitInterval(), UnitBall{Float64}())) == 2
    end

    @testset "intersect" begin
        @test IntersectDomain(0..1) == IntersectDomain((0..1))
        @test IntersectDomain{Float64}(0..1) == IntersectDomain(0.0..1.0)
        @test IntersectDomain{Float64}(0..1, 1..2) == IntersectDomain((0..1, 1..2))
        @test intersectdomain(0..1, 1..2) == Point(1)

        @test intersectdomain(0..1, 0.5..1.5) == (0..1) & (0.5..1.5)

        # intersection of productdomains
        i1 = intersectdomain((-.4..0.4)^2, (-.5 .. 0.5) × (-.1.. 0.1))
        @test i1 == productdomain(-0.4..0.4, -0.1..0.1)
        show(io,i1)
        @test String(take!(io)) == "(-0.4..0.4) × (-0.1..0.1)"
        @test intersectdomain(productdomain(UnitDisk(),-1..1), productdomain(-1..1, UnitDisk())) isa IntersectDomain
        i2 = UnitDisk() & (-.4..0.4)^2
        show(io, textmime, i2)
        @test String(take!(io)) == "UnitDisk() ∩ ((-0.4..0.4) × (-0.4..0.4))"
        @test dimension(i1) == 2
        @test dimension(i2) == 2

        i3 = ((-.5 .. 0.5) × (-.1.. 0.1)) & i2
        i4 = i2 & ((-.5 .. 0.5) × (-.1.. 0.1))
        i5 = i3 & i2
        @test SVector(0.,.05) ∈ i3
        @test SVector(0.,.05) ∈ i4
        @test SVector(0.,.75) ∉ i3
        @test SVector(0.,.75) ∉ i4

        d45 = IntersectDomain(0.0..1.5, [1.0,3.0])
        @test d45 isa Domain{Float64}
        @test 1.0 ∈ d45
        @test 1.1 ∉ d45
        @test convert(Domain{BigFloat}, d45) isa Domain{BigFloat}

        @test intersectdomain() == EmptySpace{Any}()
        @test intersectdomain(UnitDisk()) == UnitDisk()

        @test (0..1) ∩ [1.5] isa IntersectDomain{Float64}
        @test [0.5] ∩ (1..2) isa IntersectDomain{Float64}

        @test intersectdomain(0..1, [0,1]) == [0,1]
        @test intersectdomain([0,1], 0..1) == [0,1]
        @test intersectdomain([0,1], [0.0,1.0]) == [0,1]
        @test intersectdomain([0,2], [0.0,1.0]) isa IntersectDomain

        @test IntersectDomain(UnitDisk(), UnitSquare()) == IntersectDomain(UnitSquare(), UnitDisk())

        @test intersectdomain(uniondomain(0..1, 2..3), uniondomain(0..0.5, 3.5..4.5)) == 0..0.5
        @test intersectdomain(uniondomain(0..1, 2..3), 2..4) == 2..3
        @test intersectdomain(2..4, uniondomain(0..1, 2..3)) == 2..3

        # repeated intersection
        @test ncomponents(intersectdomain(UnitBall{Float64}(), UnitInterval(), UnitInterval())) == 2
        @test ncomponents(intersectdomain(UnitInterval(), UnitBall{Float64}(), UnitInterval())) == 2
        @test ncomponents(intersectdomain(UnitInterval(), UnitInterval(), UnitBall{Float64}())) == 2
        # larger intersection expressions
        @test intersectdomain(0..1, 1..3, Point(0.4), 2..5, FullSpace(), Point(-0.2)) isa EmptySpace
        @test intersectdomain(0..1, 1..3, Point(1.0)) == Point(1.0)
        @test intersectdomain(0..1, Point(1.0), 1..3) == Point(1.0)
        @test intersectdomain(Point(1.0), 0..1, 1..3) == Point(1.0)

        @test boundingbox(UnitSphere() ∩ 2UnitBall()) == (-1..1)^3
    end

    @testset "setdiff" begin
        d1 = UnitDisk() \ ProductDomain(-.5..0.5, -.1..0.1)
        @test dimension(d1) == 2
        @test SVector(0.,.74) ∈ d1
        @test SVector(0.,.25) ∈ d1
        @test SVector(0.,-.05) ∉ d1
        @test SVector(0.5, 0.1) ∉ d1
        @test SVector(1.01, 0.1) ∉ d1
        @test approx_in(SVector(1.01, 0.1), d1, 0.1)
        show(io, textmime, d1)
        @test String(take!(io)) == "UnitDisk() \\ ((-0.5..0.5) × (-0.1..0.1))"
        @test setdiff(d1, ProductDomain(-0.5..0.5, -.1..0.1)) == d1

        d2 = SetdiffDomain(0.0..3.0, [1.0, 2.5])
        @test d2 isa Domain{Float64}
        @test d2 isa SetdiffDomain
        @test 0.99 ∈ d2
        @test_throws MethodError approx_in(3.01, d2, 0.1)
        @test 1.0 ∉ d2
        @test convert(Domain{BigFloat}, d2) isa Domain{BigFloat}

        @test (0..1) \ [0.5] isa SetdiffDomain{Float64}
        d3 = [0,5] \ (0..3)
        @test d3 isa SetdiffDomain{Int}
        @test 0 ∉ d3
        @test 5 ∈ d3

        @test setdiff(0..1, 2..3) == setdiffdomain(0..1, 2..3)
        @test setdiff(0..1, 0.5) == setdiffdomain(0..1, 0.5)
        @test setdiff(0.5, 0..1) == setdiffdomain(0.5, 0..1)

        @test setdiff(0..1, EmptySpace()) == 0..1
        @test setdiff(0..1, 0.0..1.0) == EmptySpace()

        @test (0..1)^2 \ UnitCircle() == UnitInterval()^2 \ UnitCircle()
    end

    @testset "arithmetic" begin
        d1 = (0..1)
        d2 = (2..3)
        d = UnionDomain(d1) ∪ UnionDomain(d2)

        @test d .+ 1 == UnionDomain(d1 .+ 1) ∪ (d2 .+ 1)
        @test d .- 1 == UnionDomain(d1 .- 1) ∪ (d2 .- 1)
        @test 2 * d  == UnionDomain(2 * d1)  ∪ (2 * d2)
        @test d * 2 == UnionDomain(d1 * 2) ∪ (d2 * 2)
        @test d / 2 == UnionDomain(d1 / 2) ∪ (d2 / 2)
        @test 2 \ d == UnionDomain(2 \ d1) ∪ (2 \ d2)

        @test infimum(d) == minimum(d) == 0
        @test supremum(d) == maximum(d) == 3
    end

    @testset "different types" begin
        d̃1 = (0..1)
        d1 = (0f0.. 1f0)
        d2 = (2..3)

        @test UnionDomain(d1) ∪ d2 == UnionDomain(d̃1) ∪ d2
    end

    @testset "disk × interval" begin
        d = (0.0..1) × UnitDisk()
        @test SVector(0.1, 0.2, 0.3) ∈ d

        d = UnitDisk() × (0.0..1)
        @test SVector(0.1, 0.2, 0.3) ∈ d
    end

    @testset "Wrapped domain" begin
        d = 0..1.0
        w = DomainSets.WrappedDomain(d)
        @test w isa Domain
        @test 0.5 ∈ w

        d1 = DomainSets.WrappedDomain{Any}([0,5])
        @test 0 ∈ d1
        @test 1 ∉ d1

        @test convert(Domain{Float64}, [0.4,0.5]) isa Domain{Float64}
    end

    @testset "more set operations" begin
        D = UnitInterval()^3
        S = 2 * UnitBall()
        @testset "joint domain" begin
            DS = D ∪ S
            @test SA[0.0, 0.6, 0.0] ∈ DS
            @test SA[0.9, 0.6,-2.5] ∉ DS
        end

        @testset "domain intersection" begin
            DS = D ∩ S
            @test SA[0.1, 0.1, 0.1] ∈ DS
            @test SA[0.1, -0.1, 0.1] ∉ DS
        end
        @testset "domain difference" begin
            DS = D\S
            @test SA[0.1, 0.1, 0.1] ∉ DS

            D1 = 2 * D
            D2 = D * 2
            D3 = D / 2

            @test SA[2., 2., 2.] ∈ D1
            @test SA[0.9, 0.6,-2.5] ∉ D1
            @test SA[2., 2., 2.] ∈ D2
            @test SA[0.9, 0.6,-2.5] ∉ D2
            @test SA[.5, .4, .45] ∈ D3
            @test SA[.3, 0.6,-.2] ∉ D3
        end
    end
end
