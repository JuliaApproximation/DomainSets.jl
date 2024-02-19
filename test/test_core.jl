
@testset "DomainSetsCore" begin
    p = Point(0.5)
    @test !(p == 0.5)
    @test p == DomainRef(0.5)

    @test union(0.5,0.7) isa Vector{Float64}
    @test uniondomain(0.5, 0.7) isa UnionDomain{Float64,Tuple{Float64,Float64}}
    @test uniondomain(0.5, 1) isa UnionDomain{Float64,Tuple{Float64,Float64}}
    @test_throws MethodError union(Point(0.5), 0.5)
    @test union(Point(0.5), DomainRef(0.5)) == Point(0.5)
    @test union(Point(0.5), DomainRef(0.7)) isa UnionDomain{Float64}
end
