# test_specific_domains.jl

function test_specific_domains()
    @testset "$(rpad("Specific domains",80))" begin
        test_emptyspace()
        test_fullspace()
        test_interval()
        test_unitball()
        test_cube()
        test_fractals()
        test_arithmetics()
        test_tensorproduct_domain()
    end
end

function test_emptyspace()
    println("- an empty space")
    d1 = EmptySpace()
    @test 0.5 ∉ d1

    d2 = EmptySpace(Val{2}())
    @test [0.1,0.2] ∉ d2
end

function test_fullspace()
    println("- a full Euclidean space")
    d1 = FullSpace()
    @test 0.5 ∈ d1

    d2 = FullSpace(Val{2}())
    @test [0.1,0.2] ∈ d2
end

function test_interval()
    println("- intervals")

    d = Interval(0, 1)
    @test 0.5 ∈ d
    @test 1.1 ∉ d

    # Interval
    Intervala = Interval(-1.0,1.0)
    Intervala = Intervala+1
    Intervala = 1+Intervala
    @test left(Intervala) == 1
    @test left(2*Intervala) == 2
    @test right(Intervala/4) == 0.75
end

function test_unitball()
    C = Disk(2.0)
    @test in([1.4, 1.4], C)
    @test !in([1.5, 1.5], C)
    @test boundingbox(C) == BBox((-2.0,-2.0),(2.0,2.0))
    @test typeof(1.2*C)==typeof(C*1.2)
    @test in([1.5,1.5],1.2*C)
    @test in([1.5,1.5],C*1.2)

    S = Ball(2.0)
    @test [1.9,0.0,0.0] ∈ S
    @test in([0,-1.9,0.0],S)
    @test in([0.0,0.0,-1.9],S)
    @test !in([1.9,1.9,0.0],S)
    @test boundingbox(S) == BBox((-2.0,-2.0,-2.0),(2.0,2.0,2.0))
end

function test_cube()
    #Square
    D = Cube(2)
    @test [0.9, 0.9] ∈ D
    @test [1.1, 1.1] ∉ D
    @test boundingbox(D) == BBox((0,0),(1,1))

    #Cube
    D = Cube((-1.5,0.5,-3.0),(2.2,0.7,-1.0))
    @test [0.9, 0.6, -2.5] ∈ D
    @test [0.0, 0.6, 0.0] ∉ D
    @test boundingbox(D) == BBox((-1.5,0.5,-3.0),(2.2,0.7,-1.0))
end

function test_fractals()

end

function test_arithmetics()
    # joint domain
    D = Cube(3)
    S = Ball(2.0)
    DS = D ∪ S
    @test [0.0, 0.6, 0.0] ∈ DS
    @test [0.9, 0.6,-2.5] ∉ DS
    @test boundingbox(DS) == BBox((-2.0, -2.0, -2.0),(2.0, 2.0, 2.0))

    # domain intersection
    DS = D ∩ S
    @test [0.1, 0.1, 0.1] ∈ DS
    @test [0.1, -0.1, 0.1] ∉ DS
    @test boundingbox(DS) == BBox((0.0, 0.0, 0.0),(1.0, 1.0, 1.0))

    # domain difference
    DS = D-S
    @test [0.1, 0.1, 0.1] ∉ DS
    @test boundingbox(DS) == boundingbox(D)
end

function test_tensorproduct_domain()
    # ProductDomain 1
    T = tensorproduct(Interval(-1.0, 1.0) ,2)
    @test [0.5,0.5] ∈ T
    @test [-1.1,0.3] ∉ T

    # ProductDomain 2
    T = ProductDomain(Disk(1.05), Interval(-1.0,1.0))
    @test [0.5,0.5,0.8] ∈ T
    @test [-1.1,0.3,0.1] ∉ T
end
