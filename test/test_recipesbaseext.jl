module DomainSetsRecipesBaseTests
using DomainSets, StaticArrays, RecipesBase, Test

@testset "RecipesBase" begin
    @testset "unions" begin
        # each component is a series, with only the first having a legend entry
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(), UnionDomain(0..1, 2..3, 4..5))
        @test [only(r.args) for r in rd] == [0..1, 2..3, 4..5]
        @test [r.plotattributes[:primary] for r in rd] == [true, false, false]

        # attributes are passed on to every component
        u = UnionDomain(UnitDisk(), Ball(0.5, SVector(2.5, 0.0)))
        rd = RecipesBase.apply_recipe(Dict{Symbol,Any}(:linecolor => :red), u)
        @test [only(r.args) for r in rd] == collect(components(u))
        @test all(r -> r.plotattributes[:linecolor] == :red, rd)
    end
end
end # module
