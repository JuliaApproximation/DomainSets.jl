# Make sure the examples in the README continue to function

@testset "$(rpad("examples",80))" begin
    using DomainSets, StaticArrays; import DomainSets: ×
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
end
