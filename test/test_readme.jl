# Make sure the examples in the README continue to function

@testset "examples" begin
    using DomainSets, StaticArrays
    @test repr(UnitInterval()) == "$(0.0..1.0) (Unit)"
    @test repr(ChebyshevInterval()) == "$(-1.0..1.0) (Chebyshev)"
    @test repr(HalfLine()) == "$(0.0..Inf) (closed–open) (HalfLine)"

    using DomainSets: ×
    @test repr((-1..1) × (0..3) × (4.0..5.0)) == "($(-1.0..1.0)) × ($(0.0..3.0)) × ($(4.0..5.0))"
    @test SVector(1,2) in (-1..1) × (0..3)

    @test SVector(0,0,1.0) in UnitSphere(Val(3))
    @test [0.0,1.0,0.0,0.0] in UnitSphere(4)
    @test SVector(1,0) in UnitCircle()

    @test SVector(0.1,0.2,0.3) in UnitBall(Val(3))
    @test [0.1,0.2,0.3,-0.1] in UnitBall(4)
    @test SVector(0.1,0.2) in UnitDisk()

    @test 1:5 in ProductDomain([0..i for i in 1:5])
    @test ("a", 0.4) ∈ ProductDomain(["a","b"], 0..1)

    d = UnitCircle() ∪ 2UnitCircle()
    @test in.([SVector(1,0),SVector(0,2), SVector(1.5,1.5)], Ref(d)) == [1,1,0]
    d = UnitCircle() ∩ (2UnitCircle() .+ SVector(1.0,0.0))
    @test !(SVector(1,0) in d)
    @test SVector(-1,0) in d

    d = LevelSet{SVector{2,Float64}}(prod, 1.0)
    @test [0.5,2] ∈ d

    d = IndicatorFunction{Float64}( t ->  cos(t) > 0)
    @test (0.5 ∈ d, 3.1 ∈ d) == (true, false)

    d = Domain(x>0 for x in -1..1)
    @test (0.5 ∈ d, -0.5 ∈ d) == (true, false)

    d = Domain( x*y > 0 for (x,y) in UnitDisk())
    @test ([0.2, 0.3] ∈ d, [0.2, -0.3] ∈ d) == (true, false)

    d = Domain( x+y+z > 0 for (x,y,z) in ProductDomain(UnitDisk(), 0..1))
    @test [0.3,0.2,0.5] ∈ d
end
