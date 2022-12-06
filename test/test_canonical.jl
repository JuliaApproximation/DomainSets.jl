@testset "canonical domains" begin
    @test canonicaldomain(2..3) isa ChebyshevInterval
    @test canonicaldomain(2..3) == ChebyshevInterval{Float64}()
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
end
