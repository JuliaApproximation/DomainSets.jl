# test_specific_domains.jl

function test_specific_domains()
    @testset "$(rpad("Specific domains",80))" begin
        test_emptydomain()
        test_euclideanspace()
        test_interval()
        test_fractals()
    end
end

function test_emptydomain()
    println("- an empty domain")
    d1 = EmptyDomain()
    @test 0.5 ∉ d1

    d2 = EmptyDomain(Val{2}())
    @test [0.1,0.2] ∉ d2
end

function test_euclideanspace()
    println("- Euclidean space")
    d1 = EuclideanSpace()
    @test 0.5 ∈ d1

    d2 = EuclideanSpace(Val{2}())
    @test [0.1,0.2] ∈ d2
end

function test_interval()
    println("- intervals")

    d = Interval(0, 1)
    @test 0.5 ∈ d
    @test 1.1 ∉ d
end

function test_fractals()

end
