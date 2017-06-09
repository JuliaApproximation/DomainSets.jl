# test_specific_domains.jl

const v = TypeFactory{SVector}()

function test_specific_domains()
    @testset "$(rpad("Specific domains",80))" begin
        test_emptyspace()
        test_fullspace()
        test_interval()
        test_unitball()
        test_cube()
        test_arithmetics()
        test_tensorproduct_domain()
    end
end

function test_emptyspace()
    println("- an empty space")
    d1 = EmptySpace()
    @test 0.5 ∉ d1

    d2 = EmptySpace(SVector{2,Float64})
    @test v[0.1,0.2] ∉ d2
end

function test_fullspace()
    println("- a full Euclidean space")
    d1 = FullSpace()
    @test 0.5 ∈ d1

    d2 = FullSpace(SVector{2,Float64})
    @test v[0.1,0.2] ∈ d2
end

function test_interval()
    println("- intervals")

    d = Interval(0.0, 1.0)
    @test 0.5 ∈ d
    @test 1.1 ∉ d

    # Interval
    Intervala = Interval(-2.0,1.0)
    Intervala = Intervala+1
    @test leftendpoint(Intervala) == -1.0
    @test leftendpoint(2*Intervala) == -2.0
    @test rightendpoint(Intervala/4) == 0.5
end

function test_unitball()
    C = Disk(2.0)
    @test in(v[1.4, 1.4], C)
    @test !in(v[1.5, 1.5], C)
    @test typeof(1.2*C)==typeof(C*1.2)
    @test in(v[1.5,1.5],1.2*C)
    @test in(v[1.5,1.5],C*1.2)

    S = Ball(2.0)
    @test v[1.9,0.0,0.0] ∈ S
    @test in(v[0,-1.9,0.0],S)
    @test in(v[0.0,0.0,-1.9],S)
    @test !in(v[1.9,1.9,0.0],S)
end

function test_cube()
    #Square
    D = Cube(2)
    @test v[0.9, 0.9] ∈ D
    @test v[1.1, 1.1] ∉ D

    #Cube
    D = Cube((-1.5,0.5,-3.0),(2.2,0.7,-1.0))
    @test v[0.9, 0.6, -2.5] ∈ D
    @test v[0.0, 0.6, 0.0] ∉ D
end

function test_arithmetics()
    # joint domain
    D = Cube(3)
    S = Ball(2.0)
    DS = D ∪ S
    @test v[0.0, 0.6, 0.0] ∈ DS
    @test v[0.9, 0.6,-2.5] ∉ DS

    # domain intersection
    DS = D ∩ S
    @test v[0.1, 0.1, 0.1] ∈ DS
    @test v[0.1, -0.1, 0.1] ∉ DS

    # domain difference
    DS = D-S
    @test v[0.1, 0.1, 0.1] ∉ DS
end

function test_tensorproduct_domain()
    # ProductDomain 1
    T = tensorproduct(Interval(-1.0, 1.0) ,2)
    @test v[0.5,0.5] ∈ T
    @test v[-1.1,0.3] ∉ T

    # ProductDomain 2
    T = ProductDomain(Disk(1.05), Interval(-1.0,1.0))
    @test (v[0.5,0.5],0.8) ∈ T
    @test (v[-1.1,0.3],0.1) ∉ T
end
