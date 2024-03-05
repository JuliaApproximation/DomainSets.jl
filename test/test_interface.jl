
@testset "Domain interface" begin
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
