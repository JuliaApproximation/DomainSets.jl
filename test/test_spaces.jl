# test_spaces.jl

function test_spaces()
    @testset "$(rpad("Spaces",80))" begin
        test_basic_spaces()
    end
end

function test_basic_spaces()
    R = VectorSpace{1,Float64}
    R2 = VectorSpace{2,Float64}
    @test embedded(R, R2)
    @test typeof(convert_space(R2, zero(R))) == eltype(R2)

    A = ComplexSpace{BigFloat}
    B = VectorSpace{2,Float64}
    @test embedded(B, A)
    @test !embedded(A, B)
end
