
struct MyDomain{T} <: Domain{T}
    somefield::T
end

@testset "canonical domains" begin
    d1 = 2..3.0

    equal = DomainSets.Equal()
    @test !hascanonicaldomain(equal, d1)
    @test canonicaldomain(equal, d1) === d1
    @test mapfrom_canonical(equal, d1) isa IdentityMap
    @test mapto_canonical(equal, d1) isa IdentityMap
    @test mapfrom_canonical(equal, d1, 2.5) == 2.5
    @test mapto_canonical(equal, d1, 2.5) == 2.5

    @test hascanonicaldomain(d1)
    @test canonicaldomain(d1) isa ChebyshevInterval
    @test canonicaldomain(2..3) == ChebyshevInterval{Float64}()
    @test mapfrom_canonical(d1) isa AffineMap
    @test mapto_canonical(d1) isa AffineMap
    @test mapfrom_canonical(d1, 0.4) ≈ 2.7
    @test mapto_canonical(d1, 2.7) ≈ 0.4

    iso = DomainSets.Isomorphic()
    @test !hascanonicaldomain(iso, d1)
    @test canonicaldomain(iso, d1) === d1
    @test mapfrom_canonical(iso, d1) isa IdentityMap
    @test mapto_canonical(iso, d1) isa IdentityMap
    @test mapfrom_canonical(iso, d1, 2.5) == 2.5
    @test mapto_canonical(iso, d1, 2.5) == 2.5

    d2 = UnitBall(Val(1))
    @test canonicaldomain(iso, d2) == UnitBall{Float64}()
    @test mapfrom_canonical(iso, d2) isa DomainSets.NumberToVector
    @test mapto_canonical(iso, d2) isa DomainSets.VectorToNumber

    @test_throws ArgumentError mapto(UnitCircle(), UnitDisk())
    @test mapto(UnitInterval(), UnitInterval()) isa IdentityMap
    @test mapto(UnitInterval(), ChebyshevInterval()) isa AffineMap

    @test canonicaldomain(2.0..Inf) isa ChebyshevInterval  # for now

    @test canonicaldomain(OpenInterval(2.0,Inf)) === HalfLine{Float64,:open}()
    @test canonicaldomain(OpenInterval(-Inf,2.0)) === HalfLine{Float64,:open}()
    @test canonicaldomain(OpenInterval(-Inf,Inf)) === RealLine{Float64}()
    @test canonicaldomain(OpenInterval(Inf,-Inf)) === RealLine{Float64}()
    @test canonicaldomain(Interval{:closed,:open}(2.0,Inf)) === HalfLine{Float64,:closed}()
    @test canonicaldomain(Interval{:closed,:open}(2.0,4.0)) == Interval{:closed,:open}(-1.0,1.0)
    @test canonicaldomain(Interval{:open,:closed}(-Inf,2.0)) === HalfLine{Float64,:closed}()
    @test canonicaldomain(Interval{:open,:closed}(1.0,2.0)) ==  Interval{:open,:closed}(-1.0,1.0)
    @test_throws ArgumentError canonicaldomain(Interval{:closed,:open}(-Inf,Inf))
    @test_throws ArgumentError canonicaldomain(Interval{:closed,:open}(-Inf,2.0))
    @test_throws ArgumentError canonicaldomain(Interval{:open,:closed}(2.0,Inf))

    @test mapfrom_canonical(OpenInterval(2.0,Inf))(0.0) == 2.0
    @test mapfrom_canonical(OpenInterval(-Inf,2.0))(0.0) == 2.0
    @test DomainSets.isidentity(mapfrom_canonical(OpenInterval(-Inf,Inf)))
    @test mapfrom_canonical(OpenInterval(Inf,-Inf))(2.0) == -2.0
    @test mapfrom_canonical(Interval{:closed,:open}(2.0,Inf))(0.0) == 2.0
    @test mapfrom_canonical(Interval{:closed,:open}(2.0,4.0))(-1.0) == 2.0
    @test mapfrom_canonical(Interval{:open,:closed}(-Inf,2.0))(0.0) == 2.0
    @test mapfrom_canonical(Interval{:open,:closed}(1.0,2.0))(1.0) == 2.0
    @test_throws ArgumentError mapfrom_canonical(Interval{:closed,:open}(-Inf,Inf))
    @test_throws ArgumentError mapfrom_canonical(Interval{:closed,:open}(-Inf,2.0))
    @test_throws ArgumentError mapfrom_canonical(Interval{:open,:closed}(2.0,Inf))

    @test DomainSets.interval_map(-Inf,Inf,-Inf,Inf) isa StaticIdentityMap
    @test DomainSets.interval_map(Inf,-Inf,Inf,-Inf) isa StaticIdentityMap
    @test DomainSets.interval_map(-Inf,Inf,Inf,-Inf) == LinearMap(-1)
    @test DomainSets.interval_map(Inf,-Inf,-Inf,Inf) == LinearMap(-1)
    @test DomainSets.interval_map(Inf,Inf,Inf,Inf) isa StaticIdentityMap
    @test DomainSets.interval_map(-Inf,-Inf,-Inf,-Inf) isa StaticIdentityMap

    @testset "canonical types" begin
        struct MyCanonicalType <: DomainSets.CanonicalType
        end
        ctype = MyCanonicalType()
        d = 2..3
        @test canonicaldomain(ctype, d) == d
        @test mapfrom_canonical(ctype, d) isa IdentityMap
        @test mapfrom_canonical(ctype, d, 2.5) == 2.5
        @test mapto_canonical(ctype, d, 2.5) == 2.5
    end

    @testset "domainhash" begin
        mydomain = MyDomain(2.0)
        @test hash(mydomain) == DomainSets.domainhash(mydomain)
    end
end
