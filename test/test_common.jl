using DomainSets:
    convert_eltype,
    domaineltype,
    domain_numtype,
    domain_prectype

function test_components()
    @test components([1,2,3]) == ()
    d = UnionDomain(Point(1),Point(2))
    @test iscomposite(d)
    @test ncomponents(d) == length(components(d))
end

function test_eltype()
    @test convert_eltype(Float64, Point(0)) isa Point{Float64}
    @test convert_eltype(Float64, Point(0)) == Point(0)
    @test convert_eltype(Float64, [1,2]) isa Vector{Float64}
    @test convert_eltype(Float64, [1,2]) == [1,2]
    @test convert_eltype(Float64, Set([1,2])) isa Set{Float64}
    @test convert_eltype(Float64, Set([1,2])) == Set([1,2])
    @test convert_eltype(Float64, 1:5) isa AbstractVector{Float64}
    @test convert_eltype(Float64, 1:5) == 1:5
    @test convert_eltype(Float64, 1) isa Float64
    @test convert_eltype(Float64, 1) == 1
end

@testset "common functionality" begin
    @testset "elements" begin
        test_components()
    end
    @testset "eltype" begin
        test_eltype()
    end
end
