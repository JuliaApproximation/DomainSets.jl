function test_isomorphisms(T)
    m1 = FunctionMaps.VectorToNumber{T}()
    @test m1(SA[1.0]) == 1.0
    m1b = FunctionMaps.NumberToVector{T}()
    @test inverse(m1) == m1b
    @test inverse(m1, 0.4) == m1b(0.4)
    @test inverse(m1b) == m1
    @test inverse(m1b, SA[0.4]) == m1(SA[0.4])
    @test m1b(1.0) == SA[1.0]
    @test mapsize(m1) == (1,1)
    @test mapsize(m1b) == (1,)
    @test mapsize(m1 ∘ inverse(m1)) == ()
    @test mapsize(inverse(m1) ∘ m1) == (1,1)
    @test jacobian(m1 ∘ inverse(m1), 1) == one(T)
    @test jacobian(inverse(m1) ∘ m1, [1]) == ones(T,1,1)
    @test jacobian(m1) isa ConstantMap
    @test jacobian(m1b) isa ConstantMap

    m2 = FunctionMaps.VectorToComplex{T}()
    @test m2(SA[one(T), one(T)]) == 1 + im
    m2b = FunctionMaps.ComplexToVector{T}()
    @test inverse(m2) == m2b
    @test inverse(m2, zero(T)+0im) == m2b(zero(T)+0im)
    @test inverse(m2b) == m2
    @test inverse(m2b, SA[one(T),one(T)]) == one(T)+one(T)*im
    @test m2b(one(T)+one(T)*im) == SA[one(T),one(T)]
    @test mapsize(m2) == (1,2)
    @test mapsize(m2b) == (2,)
    @test mapsize(m2 ∘ m2b) == ()
    @test mapsize(m2b ∘ m2) == (2,2)

    m3 = FunctionMaps.VectorToTuple{2,T}()
    m3b = FunctionMaps.TupleToVector{2,T}()
    @test m3(SA[one(T), one(T)]) == (one(T),one(T))
    @test inverse(m3b, SA[one(T), one(T)]) == (one(T),one(T))
    @test inverse(m3) == m3b
    @test inverse(m3b) == m3
    @test m3b( (one(T),one(T)) ) == SA[one(T),one(T)]
    @test inverse(m3, (one(T),one(T)) ) == SA[one(T),one(T)]

    m4 = FunctionMaps.NestedToFlat{3,T,Tuple{Tuple{T,T},T},(2,1)}()
    m4b = FunctionMaps.FlatToNested{3,T,Tuple{Tuple{T,T},T},(2,1)}()
    x4 = ([T(1),T(2)],T(3))
    y4 = T[1,2,3]
    @test m4(x4) == y4
    @test m4b(y4) == x4
    @test inverse(m4) == m4b
    @test inverse(m4b) == m4
    @test inverse(m4, y4) == x4
    @test inverse(m4b, x4) == y4
end

function test_scaling_maps(T)
    test_generic_map(LinearMap(T(2)))
    test_generic_map(LinearMap(T(2), T(3)))
    test_generic_map(LinearMap(T(2), T(3), T(4)))
    test_generic_map(LinearMap(T(2), T(3), T(4), T(5)))
end

function test_identity_map(T)
    i1 = StaticIdentityMap{T}()
    i2 = StaticIdentityMap{SVector{2,T}}()
    test_generic_map(i1)
    test_generic_map(i2)
    @test i1 == i2
    @test hash(i1) == hash(i2)
    @test islinearmap(i1)
    @test isaffinemap(i1)
    @test convert(StaticIdentityMap{SVector{2,T}}, i1) === i2
    @test jacobian(i1) isa ConstantMap
    @test jacobian(i1, 1) == 1
    @test jacdet(i1, 1) == 1
    m1 = convert(FunctionMaps.AbstractAffineMap{T}, i1)
    @test m1 isa LinearMap{T}
    @test jacdet(m1, 1) == 1
    @test convert(FunctionMaps.AbstractAffineMap, i1) isa LinearMap{T}
    m2 = convert(FunctionMaps.LinearMap{T}, i1)
    @test m2 isa LinearMap{T}
    @test jacdet(m2, 1) == 1

    i3 = VectorIdentityMap{T}(10)
    test_generic_map(i3)
    r = rand(T, 10)
    @test i3(r) ≈ r
    @test hash(i3) == hash(VectorIdentityMap{Int}(10))

    @test IdentityMap() ∘ LinearMap(2) == LinearMap(2.0)
end

function test_basic_maps(T)
    @test UnityMap{SVector{2,Float64}}() == UnityMap{SVector{2,Float64},Float64}()
    @test hash(ConstantMap(2)) == hash(ConstantMap(2.0))
    @test FunctionMaps.absmap(ConstantMap(-2)) == ConstantMap(2)
end

function test_wrapped_maps(T)
    m1 = WrappedMap{T}(cos)
    m2 = WrappedMap{T}(sin)
    @test m1(one(T)) ≈ cos(one(T))
    @test m2(one(T)) ≈ sin(one(T))
    @test isequalmap(m1, cos)
    @test isequalmap(sin, m2)
    m3 = m1 ∘ m2
    @test m3(one(T)) ≈ cos(sin(one(T)))

    @test WrappedMap(cos) isa WrappedMap{Float64}
    @test convert(Map, cos) isa WrappedMap
    @test convert(Map{BigFloat}, m1) isa WrappedMap{BigFloat}

    @test convert(Map{T}, cos) isa FunctionMaps.WrappedMap{T,typeof(cos)}

    @test convert(Map, LinearAlgebra.I) isa GenericLinearMap{Vector{Any}}
    @test convert(Map{Vector{T}}, LinearAlgebra.I) isa GenericLinearMap{Vector{T}}
end
