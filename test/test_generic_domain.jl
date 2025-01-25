using DomainSets: factors, nfactors, factor


widen_eltype(::Type{T}) where {T<:Number} = widen(T)
widen_eltype(::Type{SVector{N,T}}) where {N,T<:Number} = SVector{N,widen(T)}
widen_eltype(::Type{Vector{T}}) where {T<:Number} = Vector{widen(T)}


# We test the generic functionality of a domain.
# These tests check whether the given domain correctly implements the
# interface of a domain.
function test_generic_domain(d)
    @test isrealdomain(d) == isrealtype(domaineltype(d))
    @test isrealdomain(d) == isrealtype(domain_numtype(d))

    if d isa Domain
        @test convert(Domain{eltype(d)}, d) == d
        @test convert(Domain{widen_eltype(eltype(d))}, d) == d
    elseif DomainStyle(d) isa IsDomain
        @test convert_eltype(eltype(d), d) == d
        @test convert_eltype(widen_eltype(eltype(d)), d) == d
    end
    @test prectype(convert_prectype(BigFloat, 2)) == BigFloat

    if !isempty(d)
        x = choice(d)
        @test x ∈ d
        @test approx_in(x, d, 0.01)
        @test_throws ErrorException approx_in(x, d, -1)
        @test in.([x,x], d) == in.([x,x], Ref(d))
        @test approx_in.([x,x], d, 0.01) == approx_in.([x,x], Ref(d), 0.01)
    else
        try
            x = choice(d)
            @test false
        catch
        end
    end
    @test isequaldomain(canonicaldomain(DomainSets.Equal(), d), d)
    @test hash(d) == hash(DomainSets.simplify(d))
    if hascanonicaldomain(d)
        cd = canonicaldomain(d)
        @test mapfrom_canonical(d) == mapto(cd, d)
        @test mapto_canonical(d) == mapto(d, cd)
        x1 = choice(cd)
        @test mapfrom_canonical(d, x1) ∈ d
        @test mapto_canonical(d, x) ∈ cd
    else
        @test mapto_canonical(d) == IdentityMap{eltype(d)}(dimension(d))
        @test mapfrom_canonical(d) == IdentityMap{eltype(d)}(dimension(d))
    end
    if hasparameterization(d)
        par = parameterdomain(d)
        @test mapfrom_parameterdomain(d) == mapto(par, d)
        xp = choice(par)
        @test approx_in(mapfrom_parameterdomain(d, xp), d)
    end
    if iscomposite(d)
        @test ncomponents(d) == length(components(d))
        els = components(d)
        @test all([component(d,i) == els[i] for i in 1:ncomponents(d)])
        if d isa ProductDomain
            @test factors(d) == components(d)
            @test nfactors(d) == ncomponents(d)
            if nfactors(d) > 0
                @test factor(d, 1) == component(d, 1)
            end
        end
    end
end

@testset "generic domain interface" begin
    domains = [
        0..1,
        UnitInterval(),
        ChebyshevInterval(),
        HalfLine(),
        NegativeHalfLine(),
        UnitInterval()^3,
        (0..1) × (2.0..3) × (3..4),
        UnitBall(),
        VectorUnitBall(),
        VectorUnitBall(8),
        UnitDisk(),
        VectorUnitDisk(),
        UnitCircle(),
        VectorUnitCircle(),
        UnitSphere(),
        VectorUnitSphere(),
        UnitSimplex(Val(2)),
        UnitSimplex(2),
        DomainSets.WrappedDomain(0..2.0)
    ]

    @testset "generic domains" begin
        for domain in domains
            test_generic_domain(domain)
        end
    end

    @testset "generic functionality" begin
        struct SomeDomain <: Domain{Float64}
        end
        @test_throws MethodError 0.5 ∈ SomeDomain()
        @test_throws MethodError approx_in(0.5, SomeDomain())

        if VERSION < v"1.6-"
            @test_logs (:warn, "`in`: incompatible combination of point: Tuple{Float64,Float64} and domain eltype: SArray{Tuple{2},Float64,1,2}. Returning false.") (0.5,0.2) ∈ UnitCircle()
            @test_logs (:warn, "`in`: incompatible combination of point: Tuple{Float64,Float64} and domain eltype: SArray{Tuple{2},Float64,1,2}. Returning false.") approx_in((0.5,0.2), UnitCircle())
        else
            @test_logs (:warn, "`in`: incompatible combination of point: Tuple{Float64, Float64} and domain eltype: SVector{2, Float64}. Returning false.") (0.5,0.2) ∈ UnitCircle()
            @test_logs (:warn, "`in`: incompatible combination of point: Tuple{Float64, Float64} and domain eltype: SVector{2, Float64}. Returning false.") approx_in((0.5,0.2), UnitCircle())
        end

        # functionality using broadcast
        @test 2 * (1..2) == 2 .* (1..2)
        @test (1..2) * 2 == (1..2) .* 2
        @test (1..2) / 2 ≈ 0.5..1
        @test 2 \ (1..2) ≈ 0.5..1
        @test all(rand(4) .∈ (-1..1))
        @test all(approx_in.(rand(4), (-1..1)))
        @test all(approx_in.([1.005,1.0005], (-1..1), [1e-2,1e-3]))
        @test all(rand(4) .∉ (2..3))
        @test_throws MethodError (0..1) + 0.4

        # promotion
        @test DomainSets.promote_domains() == ()
        @test DomainSets.promote_domains(0..1, 2..4.0) isa Tuple{ClosedInterval{Float64},ClosedInterval{Float64}}
        @test DomainSets.promote_domains(0..1.0, [1,2,3]) isa Tuple{Interval,Vector{Float64}}
        @test DomainSets.promote_domains([1,2,3], 0..1.0) isa Tuple{Vector{Float64},Interval}

        # compatible point-domain pairs
        @test DomainSets.iscompatiblepair(0.5, 0..1)
        @test DomainSets.iscompatiblepair(0.5, 0..1.0)
        @test !DomainSets.iscompatiblepair(0.5, UnitBall())
        @test DomainSets.iscompatiblepair(0.5, [0.3])
        @test DomainSets.iscompatiblepair(0, [0.3])
        @test DomainSets.iscompatiblepair(0.5, [1, 2, 3])
        @test DomainSets.iscompatiblepair(0, Set([1, 2, 3]))
        @test DomainSets.iscompatiblepair(0.5, Set([1, 2, 3]))
        @test DomainSets.iscompatiblepair(0, Set([1.0, 2, 3]))

        @test DomainSets.convert_eltype(Float64, Set([1,2])) isa Set{Float64}
        @test_throws ArgumentError DomainSets.convert_eltype(Float64, (1,2)) isa NTuple{2,Float64}
        @test DomainSets.convert_eltype(Int, (1,2)) == (1,2)
    end

    @test choice(Set([1,2,3])) ∈ Set([1,2,3])
    @test choice([1,2,3]) == 1
end
