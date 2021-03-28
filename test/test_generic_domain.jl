
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
        @test_throws ErrorException approx_in(x, d, -1)
    else
        try
            x = point_in_domain(d)
            @test false
        catch
        end
    end
    @test convert(Domain{eltype(d)}, d) == d
    @test convert(Domain{widen_eltype(eltype(d))}, d) == d
end

@testset "generic domains" begin
    domains = [
        0..1,
        UnitInterval(),
        ChebyshevInterval(),
        HalfLine(),
        NegativeHalfLine(),
        UnitInterval()^3,
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
        WrappedDomain(0..1.0)
    ]

    @testset "$(rpad("generic domains",80))" begin
        for domain in domains
            test_generic_domain(domain)
        end
    end
end
