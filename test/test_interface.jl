struct InheritedDomain <: Domain{Int} end

struct InterfaceDomain end
DomainSets.domaineltype(::InterfaceDomain) = Int
DomainSets.DomainStyle(::Type{InterfaceDomain}) = IsDomain()

struct NonDomain end

@testset "Domain interface" begin
    d1 = InheritedDomain()
    @test d1 isa Domain
    # @test eltype(d1) == Any   # to be changed in breaking release
    @test eltype(d1) == Int
    @test domaineltype(d1) == Int
    @test domaineltype(DomainRef(d1)) == Int
    @test DomainStyle(d1) == IsDomain()
    @test checkdomain(d1) === d1
    @test domain(d1) === d1

    d2 = InterfaceDomain()
    @test !(d2 isa Domain)
    @test domaineltype(d2) == Int
    @test domaineltype(DomainRef(d2)) == Int
    @test domaineltype(DomainRef(d2)) == domaineltype(d2)
    @test DomainStyle(d2) == IsDomain()
    @test DomainRef(d2) isa DomainRef
    @test domain(DomainRef(d2)) == d2
    @test checkdomain(d2) === d2
    @test checkdomain(DomainRef(d2)) === d2

    d3 = NonDomain()
    @test !(d3 isa Domain)
    @test DomainStyle(d3) == NotDomain()
    @test_throws ErrorException checkdomain(d3)

    @test DomainStyle(2.0) == IsDomain()
    @test DomainStyle([1,2]) == IsDomain()
    @test DomainStyle(Set([1,2])) == IsDomain()

    s = "string"
    @test domaineltype(s) == eltype(s)
    p = Point(0.5)
    @test !(p == 0.5)
    @test p == DomainRef(0.5)

    @test union(0.5,0.7) isa Vector{Float64}
    @test uniondomain(0.5, 0.7) == [0.5, 0.7]
    @test uniondomain(0.5, 1) == [0.5, 1]
    # @test_throws MethodError union(Point(0.5), 0.5)
    @test union(Point(0.5), DomainRef(0.5)) == Point(0.5)
    @test union(Point(0.5), DomainRef(0.7)) isa UnionDomain{Float64}
end
