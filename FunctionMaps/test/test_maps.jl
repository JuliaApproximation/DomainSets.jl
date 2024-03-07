
maps_to_test(T) = [
    StaticIdentityMap{T}(),
    VectorIdentityMap{T}(10),
    ConstantMap{T}(one(T)),
    ConstantMap{T}(SVector{2,T}(1,2)),
    ZeroMap{T}(),
    UnityMap{T}(),
    AffineMap(T(1.2), T(2.4)), # scalar map
    AffineMap(-T(1.2), T(2.4)), # scalar map with negative A
    AffineMap(randvec(T, 2, 2), randvec(T, 2)), # static map
    AffineMap(randvec(T, 3, 2), randvec(T, 3)), # static map, rectangular
    AffineMap(rand(T, 2, 2), rand(T, 2)), # vector map
    AffineMap(rand(T, 3, 2), rand(T, 3)), # vector map, rectangular
    AffineMap(LinearAlgebra.I, one(T)/2),  # use UniformScaling object as A
    AffineMap(LinearAlgebra.I, randvec(T,2)),  # use UniformScaling object as A
    GenericAffineMap(randvec(T, 2, 2), randvec(T, 2)),
    GenericAffineMap(T(1.2), randvec(T, 2)),
    GenericAffineMap(randvec(T, 3, 2), randvec(T, 3)),
    Translation(randvec(T, 3)),
    LinearMap(randvec(T, 2, 2)),
    LinearMap(randvec(T, 2)),
    LinearMap(randvec(T, 2, 2)) ∘ AffineMap(T(1.2), randvec(T, 2)),
    AffineMap(5.0, 2.0) ∘ VectorToComplex{T}() ∘ UnitCircleMap{T}(),
    LinearMap(SMatrix{2,2}(1,2,3,T(4))) ∘ CartToPolarMap{T}() ∘ LinearMap(SMatrix{2,2}(1,2,3,T(4))),
    # Interval{Any}(0.0, 1.0)
]

randvec(T,n) = SVector{n,T}(rand(n))
randvec(T,m,n) = SMatrix{m,n,T}(rand(m,n))

suitable_point_to_map(m::Map) = suitable_point_to_map(m, domaintype(m))

suitable_point_to_map(m::Map, ::Type{SVector{N,T}}) where {N,T} = SVector{N,T}(rand(N))
suitable_point_to_map(m::Map, ::Type{T}) where {T<:Number} = rand(T)
suitable_point_to_map(m::Map, ::Type{<:AbstractVector{T}}) where {T} = rand(T, mapsize(m,2))

suitable_point_to_map(m::ProductMap) =
    map(suitable_point_to_map, components(m))
suitable_point_to_map(m::VcatMap{T,M,N}) where {T,M,N} =
    SVector{N,T}(rand(T,N))

suitable_point_to_map(::CartToPolarMap{T}) where {T} = randvec(T,2)
suitable_point_to_map(::PolarToCartMap{T}) where {T} = randvec(T,2)

widertype(T) = widen(T)
widertype(::Type{SVector{N,T}}) where {N,T} = SVector{N,widen(T)}
widertype(::Type{Vector{T}}) where {T} = Vector{widen(T)}
widertype(::Type{Tuple{A}}) where {A} = Tuple{widen(A)}
widertype(::Type{Tuple{A,B}}) where {A,B} = Tuple{widen(A),widen(B)}
widertype(::Type{Tuple{A,B,C}}) where {A,B,C} = Tuple{widen(A),widen(B),widen(C)}
widertype(::Type{NTuple{N,T}}) where {N,T} = NTuple{N,widen(T)}

issquarematrix(A) = false
issquarematrix(A::AbstractArray) = size(A,1)==size(A,2)


function test_maps()
    @testset "generic functionality" begin
        test_generic_functionality()
    end

    @testset "generic map tests" begin
        generic_map_tests(Float64)
        generic_map_tests(BigFloat)
    end

    # Test special maps
    @testset "identity map" begin
        test_identity_map(Float64)
        test_identity_map(BigFloat)
    end
    @testset "basic maps" begin
        test_basic_maps(Float64)
        test_basic_maps(BigFloat)
    end
    @testset "affine maps" begin
        test_affine_maps(Float64)
        test_affine_maps(BigFloat)
    end
    @testset "composite maps" begin
        test_composite_map(Float64)
        test_composite_map(BigFloat)
    end
    @testset "product maps" begin
        test_product_map(Float64)
        test_product_map(BigFloat)
    end
    @testset "wrapped maps" begin
        test_wrapped_maps(Float64)
        test_wrapped_maps(BigFloat)
    end
    @testset "scaling maps" begin
        test_scaling_maps(Float64)
        test_scaling_maps(BigFloat)
    end
    @testset "isomorphisms" begin
        test_isomorphisms(Float64)
        test_isomorphisms(BigFloat)
    end
    @testset "Mixed maps" begin
        test_mixed_maps()
    end
end

function test_composite_map(T)
    a = T(0)
    b = T(1)
    c = T(2)
    d = T(3)
    ma = StaticIdentityMap{T}()
    mb = interval_map(a, b, c, d)

    r = suitable_point_to_map(ma)
    m1 = ma∘mb
    test_generic_map(m1)
    @test m1(r) ≈ ma(mb(r))
    m2 = m1∘mb
    test_generic_map(m2)
    @test m2(r) ≈ m1(mb(r))
    m3 = mb∘m2
    test_generic_map(m3)
    @test m3(r) ≈ mb(m2(r))
    m = m2∘m3
    test_generic_map(m)
    @test m(r) ≈ m2(m3(r))

    m5 = ComposedMap(LinearMap(rand(T,2,2)), AffineMap(rand(T,2,2),rand(T,2)))
    test_generic_map(m5)
    @test jacobian(m5) isa ConstantMap
    @test m5[Component(1)] isa LinearMap
    @test m5[Component(2)] isa AffineMap
    @test ComposedMap(m5[Component(1:2)]...) == m5
    @test_throws BoundsError m5[Component(3)]

    m6 = multiply_map(ma,ma)
    @test m6(one(T)/2) == one(T)/4
    @test jacobian(m6) isa SumMap
    @test jacobian(m6)(one(T)) == 2
    @test jacobian(m6, one(T)) == 2

    m7 = sum_map(ma,ma)
    @test m7(one(T)) == 2
    @test jacobian(m7) isa ConstantMap
    @test jacobian(m7, one(T)) == 2

    @test mapsize(ComposedMap(LinearMap(2),LinearMap(rand(T,2)),LinearMap(rand(T,2,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(rand(T,2)),LinearMap(rand(T,2)'),LinearMap(rand(T,2)))) == (2,)
    @test mapsize(ComposedMap(LinearMap(rand(T,2,2)),LinearMap(rand(T,2)'),LinearMap(2))) == (1,2)
    @test mapsize(ComposedMap(LinearMap(one(T)),LinearMap(one(T)))) == ()

    @test composedmap() == ()
    @test composedmap(ma) == ma
    @test composedmap(ma,ma) == ma
    @test composedmap(ma,ma,ma) == ma

    @test composite_jacobian(ma) == jacobian(ma)

    @test multiply_map() == ()
    @test multiply_map(ma) == ma
    @test multiply_map(ma,ma)(2*one(T)) == 4
    @test multiply_map(ma,ma,ma)(2*one(T)) == 8

    @test sum_jacobian() == ()
    @test sum_jacobian(ma) == jacobian(ma)
end

function test_mixed_maps()
    m1 = composedmap(cos, sin)
    @test domaintype(m1) == Any
    @test m1(0.4) == sin(cos(0.4))

    m2 = composedmap(AffineMap(2.0, 3.0), cos)
    @test domaintype(m2) == Float64
    @test m2(0.4) ≈ cos(2*0.4+3)
    @inferred m2(0.4)
    @test repr(m2) == "cos ∘ (x -> 2.0 * x + 3.0)"

    m3 = composedmap(sin, AffineMap(2.0, 3.0), cos)
    @test m3(0.5) ≈ cos(2*sin(0.5)+3)
    @test repr(m3) == "cos ∘ (x -> 2.0 * x + 3.0) ∘ sin"
    @test domaintype(m3) == Any
    @inferred m3(0.5)

    m4 = productmap(sin, cos)
    @test m4 isa TupleProductMap
    @test domaintype(m4) == Tuple{Any,Any}
    @test m4(0.3,0.5) == (sin(0.3), cos(0.5))
end

test_maps()
