# test_spaces.jl

function test_spaces()
    @testset "$(rpad("Spaces",80))" begin
        test_basic_spaces()
        test_product_spaces()
    end
end

######################################
# Helper functions and generic tests
######################################

nonzero_element(::Type{T}) where {T <: Number} = one(T)
nonzero_element(::Type{SVector{N,T}}) where {N,T<:Number} = ones(SVector{N,T})
nonzero_element(::Type{Tuple{T,S}}) where {T,S} = (nonzero_element(T),nonzero_element(S))
nonzero_element(::Type{Tuple{T,S,U}}) where {T,S,U} = (nonzero_element(T),nonzero_element(S),nonzero_element(U))
nonzero_element(::Type{Tuple{T,S,U,V}}) where {T,S,U,V} = (nonzero_element(T),nonzero_element(S),nonzero_element(U),nonzero_element(V))

# Generic tests for isomorphisms
function test_isomorphism(A, B)
    @test isomorphic(A, B)
    @test isomorphic(B, A)
    test_embedding(A, B)
    test_embedding(B, A)
    T = eltype(A)
    S = eltype(B)
    x = nonzero_element(T)
    y = nonzero_element(S)
    # Is conversion invertible?
    @test convert_space(A, convert_space(B, x)) == x
    @test convert_space(B, convert_space(A, y)) == y
end

# Generic tests for embeddings
function test_embedding(A, B)
    @test embedded(A, B)
    T = eltype(A)
    S = eltype(B)
    x = nonzero_element(T)
    y = nonzero_element(S)
    # Does conversion work? We verify if it yields the right type.
    @test typeof(convert_space(B, x)) == S
    # And we test the left inverse
    @test restrict_space(A, convert_space(B, x)) == x
end


##################
# Specific tests
##################

function test_basic_spaces()
    Z = IntegerSpace{Int}
    Q = RationalSpace{Int}
    R = RealSpace{Float64}
    C = ComplexSpace{Float64}
    R1 = VectorSpace{1,Float64}
    R2 = VectorSpace{2,Float64}
    R3 = VectorSpace{3,Float64}
    @test ℤ == Z
    @test ℚ == Q
    @test ℝ == R
    @test ℝ1 == R1
    @test ℝ2 == R2
    @test ℝ3 == R3

    # The basic spaces are not isomorphic to each other
    @test !isomorphic(Z, Q)
    @test !isomorphic(Q, R)
    @test !isomorphic(R1, R2)
    @test !isomorphic(R2, R3)
    test_embedding(Z, Q)
    test_embedding(Z, R)
    test_embedding(Z, C)
    test_embedding(Q, R)
    test_embedding(Q, C)
    test_embedding(R, C)
    test_embedding(R, R1)
    test_embedding(R1, R2)
    test_embedding(R1, R3)

    # Isomorphism between T and SVector{1,T}
    test_isomorphism(Z, VectorSpace{1,Int})
    test_isomorphism(Q, VectorSpace{1,Rational{Int}})
    test_isomorphism(R, R1)

    # Isomorphism between C and R2
    test_isomorphism(C, R2)

    # Some tests with BigFloat and BigInt
    @test widen(Z) == IntegerSpace{widen(Int)}
    @test widen(Q) == RationalSpace{widen(Int)}
    @test widen(C) == ComplexSpace{widen(Float64)}
    @test widen(R) == RealSpace{widen(Float64)}
    @test widen(R3) == VectorSpace{3,widen(Float64)}
    test_embedding(Z, widen(Z))
    test_embedding(Q, widen(Q))
    test_embedding(C, widen(C))
    test_embedding(R, widen(R))
    test_embedding(R1, widen(R3))

    test_isomorphism(widen(R), widen(R1))
    test_isomorphism(widen(C), widen(R2))

    # Combination of isomorphism and widening
    test_embedding(R, widen(R1))
    test_embedding(R1, widen(R))
    test_embedding(C, widen(R2))
    test_embedding(R2, widen(C))
    @test !embedded(widen(R), R1)
    @test !embedded(widen(C), R2)
    @test !embedded(widen(R1), R2)
    @test !embedded(widen(R1), R3)

    # Some promotions
    @test typeof(promote_space(1.0)) == Tuple{typeof(1.0)}
    @test typeof(promote_space(1, 1.0)) == Tuple{typeof(1.0),typeof(1.0)}

    @test promote_space_type(RealSpace{Float64}, RealSpace{BigFloat}) == RealSpace{BigFloat}
    @test promote_space_type(VectorSpace{1,BigFloat}, VectorSpace{2,Float64}) == VectorSpace{2,BigFloat}
    @test promote_space_type(VectorSpace{1,Float64}, VectorSpace{2,BigFloat}) == VectorSpace{2,BigFloat}

    # Preference for scalars over 1D vectors
    @test promote_space_type(RealSpace{Float64}, VectorSpace{1,Float64}) == RealSpace{Float64}
    @test promote_space_type(VectorSpace{1,Float64}, RealSpace{Float64}) == RealSpace{Float64}
    # - even when one space is bigger
    @test promote_space_type(RealSpace{Float64}, VectorSpace{1,BigFloat}) == RealSpace{BigFloat}
    @test promote_space_type(RealSpace{BigFloat}, VectorSpace{1,Float64}) == RealSpace{BigFloat}

    # Preference for 2D vectors over complex numbers
    @test promote_space_type(ComplexSpace{Float64}, VectorSpace{2,Float64}) == VectorSpace{2,Float64}
    @test promote_space_type(VectorSpace{2,Float64}, ComplexSpace{Float64}) == VectorSpace{2,Float64}
    # - except when one space is bigger
    @test promote_space_type(ComplexSpace{BigFloat}, VectorSpace{2,Float64}) == VectorSpace{2,BigFloat}
    @test promote_space_type(ComplexSpace{Float64}, VectorSpace{2,BigFloat}) == VectorSpace{2,BigFloat}
end

function test_product_spaces()
    Z = IntegerSpace{Int}
    R = RealSpace{Float64}
    P = Z × R
    @test eltype(P) == Tuple{Int,Float64}

    R2 = R × R
    @test eltype(R2) == Tuple{Float64,Float64}
    test_isomorphism(R2, VectorSpace{2,Float64})

    R3 = cartesianproduct(R, R, R)
    @test eltype(R3) == Tuple{Float64,Float64,Float64}
    test_isomorphism(R3, VectorSpace{3,Float64})

    RR = R2 × R
    @test eltype(RR) == Tuple{Tuple{Float64,Float64},Float64}
    test_isomorphism(RR, VectorSpace{3,Float64})
end
