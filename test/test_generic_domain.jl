
widen_eltype(::Type{T}) where {T<:Number} = widen(T)
widen_eltype(::Type{SVector{N,T}}) where {N,T<:Number} = SVector{N,widen(T)}
widen_eltype(::Type{Vector{T}}) where {T<:Number} = Vector{widen(T)}


# We test the generic functionality of a domain.
# These tests check whether the given domain correctly implements the
# interface of a domain.
function test_generic_domain(d::Domain)
    @test isreal(d) == isreal(eltype(d))
    @test isreal(d) == isreal(numtype(d))

    if !isempty(d)
        x = point_in_domain(d)
        @test x ∈ d
        @test approx_in(x, d, 0.01)
        @test_throws ErrorException approx_in(x, d, -1)
        @test in.([x,x], d) == in.([x,x], Ref(d))
        @test approx_in.([x,x], d, 0.01) == approx_in.([x,x], Ref(d), 0.01)
    else
        try
            x = point_in_domain(d)
            @test false
        catch
        end
    end
    if canonicaldomain(d) == d
        @test tocanonical(d) == IdentityMap{eltype(d)}()
        @test fromcanonical(d) == IdentityMap{eltype(d)}()
    else
        cd = canonicaldomain(d)
        @test fromcanonical(d) == bijection(cd, d)
        @test tocanonical(d) == bijection(d, cd)
    end
    @test convert(Domain{eltype(d)}, d) == d
    @test convert(Domain{widen_eltype(eltype(d))}, d) == d
    @test prectype(convert_prectype(d, BigFloat)) == BigFloat
end

@testset "generic domains" begin
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
        UnitSimplex{2}(),
        VectorUnitSimplex(2),
        WrappedDomain(0..2.0)
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

        # some functionality in broadcast
        @test 2 * (1..2) == 2 .* (1..2)
        @test (1..2) * 2 == (1..2) .* 2
        @test (1..2) / 2 ≈ (0.5..1)
        @test 2 \ (1..2) ≈ (0.5..1)
    end
end
