
struct MySimpleMap end
FunctionMaps.MapStyle(::Type{MySimpleMap}) = IsMap()

@testset "map interface" begin
    @test MapStyle(2.0) isa NotMap
    @test MapStyle(2.0) == MapStyle(typeof(2.0))
    @test_throws ArgumentError checkmap(2.0)
    m = LinearMap(2)
    @test MapStyle(m) isa IsMap
    @test functionmap(m) === m
    @test checkmap(m) === m
    m2 = MySimpleMap()
    @test MapStyle(m2) isa IsMap
    @test checkmap(m2) == m2
end
