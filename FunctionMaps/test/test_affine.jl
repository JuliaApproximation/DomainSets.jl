function test_affine_maps(T)
    A = rand(T,2,2)
    @test FunctionMaps.to_matrix(Vector{T}, A) == A
    @test FunctionMaps.to_matrix(T, 2) == 2
    @test FunctionMaps.to_matrix(SVector{2,T}, 2) == SMatrix{2,2}(2,0,0,2)
    @test FunctionMaps.to_matrix(SVector{2,T}, LinearAlgebra.I) == SMatrix{2,2}(1,0,0,1)
    @test FunctionMaps.to_matrix(Vector{T}, 2) == UniformScaling(2)
    @test FunctionMaps.to_matrix(Vector{T}, LinearAlgebra.I) == LinearAlgebra.I
    # test fallback with nonsensical call
    @test FunctionMaps.to_matrix(Tuple{Int}, 2) == 2

    @test FunctionMaps.to_matrix(T, A, 2) == A
    @test FunctionMaps.to_matrix(T, 2, 3) == 2
    @test FunctionMaps.to_matrix(T, UniformScaling(2), 3) == 2
    @test FunctionMaps.to_matrix(T, LinearAlgebra.I, zero(T)) isa T
    @test FunctionMaps.to_matrix(SVector{2,T}, 2, SVector(1,1)) == SMatrix{2,2}(2,0,0,2)
    @test FunctionMaps.to_matrix(Vector{T}, 2, [1,2]) == [2 0 ; 0 2]

    @test FunctionMaps.to_vector(T, 2) == 0
    @test FunctionMaps.to_vector(SVector{2,T}, 2) == SVector(0,0)
    @test FunctionMaps.to_vector(Vector{T}, A) == [0,0]
    @test FunctionMaps.to_vector(T, 2, 3) == 3

    if T != BigFloat    # BigFloat's make pinv fail for StaticArrays
        @test FunctionMaps.matrix_pinv(SMatrix{2,2}(rand(T),rand(T),rand(T),rand(T))) isa SMatrix{2,2}
        @test FunctionMaps.matrix_pinv(SVector(rand(T),rand(T))) isa Transpose{T,SVector{2,T}}
    end

    test_linearmap(T)
    test_translation(T)
    test_affinemap(T)
    test_abstractarraymap(T)
end



function test_linearmap(T)
    m1 = LinearMap(2one(T))
    @test m1 isa LinearMap{T}
    @test domaintype(m1) == T
    @test islinearmap(m1)
    @test isaffinemap(m1)
    @test m1(3) == 6
    @test m1(3one(T)) == 6
    @test m1(3one(widen(T))) isa widen(T)
    @test m1(SVector(1,2)) == SVector(2,4)
    @test m1([1.0,2.0]) == [2.0,4.0]
    mw1 = convert(Map{widen(T)}, m1)
    @test mw1 isa Map{widen(T)}
    @test mw1.A isa widen(T)
    @test jacobian(m1) isa ConstantMap{T}
    @test jacobian(m1, 1) == 2
    @test LinearMap(one(T)) == StaticIdentityMap{T}()

    m2 = LinearMap(2)
    @test domaintype(m2) == Int
    @test m2 === ScalarLinearMap(2)
    @test m2(one(T)) isa T
    @test m2 == convert(Map, 2)
    @test jacobian(m2, 1) == 2
    @test jacobian(m2) isa ConstantMap{Int}
    @test jacobian(m2, 1) == 2
    @test LinearMap(1) == StaticIdentityMap{T}()

    m3 = LinearMap(SMatrix{2,2}(one(T), 2one(T), 3one(T), 4one(T)))
    @test m3 isa LinearMap{SVector{2,T}}
    @test m3 === StaticLinearMap(m3.A)
    @test m3(SVector(1,2)) == SVector(7, 10)
    @test m3(SVector{2,T}(1,2)) == SVector{2,T}(7, 10)
    @test m3 ∘ m3 isa LinearMap
    @test LinearMap(SA[1 0; 0 1]) == StaticIdentityMap{domaintype(m3)}()

    A = rand(T,2,2)
    m4 = LinearMap(A)
    @test m4 isa LinearMap{Vector{T}}
    @test m4([1,2]) ==  A * [1,2]
    y = zeros(T,2)
    @test (FunctionMaps.applymap!(y, m4, [1,2]); y == A * [1,2])
    @test jacobian(m4, [1,2]) == A

    m5 = LinearMap{Vector{T}}(UniformScaling(2*one(T)))
    @test m5 isa LinearMap{Vector{T}}
    @test m5([1,2]) ==  2 * [1,2]
    @test jacobian(m5, [1,2]) == UniformScaling(2)

    m6 = LinearMap{SVector{2,T}}(2)
    @test m6 isa GenericLinearMap{SVector{2,T},T}
    @test m6(SA[one(T),one(T)]) == [2,2]
    @test leftinverse(m6) isa LinearMap{SVector{2,T}}
    @test leftinverse(m6)([2,2]) == [1,1]

    # Test construction and conversion
    @test LinearMap{T}(1) isa ScalarLinearMap{T}
    @test LinearMap{Vector{T}}(rand(Int,5,5)) isa VectorLinearMap{T}
    @test LinearMap{SVector{2,T}}(SMatrix{2,2}(1, 2, 3, 4)) isa StaticLinearMap{T,2,2,4}

    @test GenericLinearMap(one(T)) isa GenericLinearMap{T,T}
    @test GenericLinearMap(SMatrix{2,2,T}(1,1,1,1)) isa GenericLinearMap{SVector{2,T},SMatrix{2,2,T,4}}

    @test VectorLinearMap(rand(T,3,3)) isa VectorLinearMap{T}
    @test StaticLinearMap(SMatrix{2,2,T}(1,2,3,4)) isa StaticLinearMap{T,2,2,4}

    @test convert(Map{SVector{2,T}}, LinearMap(rand(T,2,2))) isa StaticLinearMap{T,2,2,4}
    @test convert(Map{T}, 1) isa LinearMap{T}
end



function test_translation(T)
    v = randvec(T,3)
    m = Translation(v)
    @test !islinearmap(m)
    @test isaffinemap(m)
    @test inverse(inverse(m)) == m
    @test jacobian(m) isa ConstantMap
    @test affinevector(jacobian(m)) == [1 0 0; 0 1 0; 0 0 1]
    @test jacdet(m) isa UnityMap{SVector{3,T},T}

    @test mapsize(Translation(one(T))) == ()

    # Test construction and conversion
    @test Translation(one(T)) isa ScalarTranslation{T}
    @test Translation{T}(1) isa ScalarTranslation{T}
    @test Translation(SVector{2,T}(1,2)) isa StaticTranslation{T,2}
    @test Translation{SVector{2,T}}(rand(T,2)) isa StaticTranslation{T,2}
    @test Translation(rand(Int,2)) isa VectorTranslation{Int}
    @test Translation{Vector{T}}(rand(Int,2)) isa VectorTranslation{T}
    @test Translation(1:10) isa GenericTranslation{Vector{Int}}
    @test Translation{Vector{T}}(1:10) isa GenericTranslation{Vector{T}}

    @test StaticTranslation(MVector(one(T),2)) isa StaticTranslation{T,2}
    @test StaticTranslation{T}(SVector(1,2)) isa StaticTranslation{T,2}

    @test VectorTranslation(1:10) isa VectorTranslation{Int}
    @test VectorTranslation{T}(1:10) isa VectorTranslation{T}

    @test GenericTranslation(1:10) isa GenericTranslation{Vector{Int}}
    @test GenericTranslation{Vector{T}}(1:10) isa GenericTranslation{Vector{T}}

    @test GenericTranslation(one(T)) isa GenericTranslation{T,T}
    @test GenericTranslation{T}(1) isa GenericTranslation{T,T}
    @test_throws InexactError GenericTranslation{Int}(one(T)/2)

    @test convert(Map{SVector{2,T}}, Translation(rand(T,2))) isa StaticTranslation{T,2}
end


function test_affinemap(T)
    m1 = AffineMap(T(2), T(3))
    @test m1 isa ScalarAffineMap{T}
    @test !islinearmap(m1)
    @test isaffinemap(m1)
    @test m1(2) == 7

    @test convert(ScalarAffineMap{BigFloat}, m1) == m1

    @test m1 ∘ m1 isa AffineMap
    @test (m1 ∘ m1)(2) == 2*(2*2+3)+3

    m2 = AffineMap(T(2), SVector{2,T}(1,2))
    @test m2 isa GenericAffineMap{SVector{2,T}}
    @test m2(SVector(1,2)) == SVector(3,6)
    @test mapsize(m2) == (2,2)

    m3 = AffineMap(UniformScaling(2*one(T)), [one(T),2*one(T)])
    @test m3 isa AffineMap{Vector{T}}
    @test mapsize(m3) == (2,2)
    @test m3([1,2]) ==  2 * [1,2] + [1,2]
    y = zeros(T,2)
    @test (FunctionMaps.applymap!(y, m3, [1,2]); y == m3([1,2]))
    @test jacobian(m3, [1,2]) == [2 0; 0 2]
    @test jacdet(m3, [1,2]) == 4

    @test AffineMap(LinearAlgebra.I, 2.0) isa ScalarAffineMap
    m4 = GenericAffineMap(LinearAlgebra.I, 2.0)
    @test m4.A isa UniformScaling{Float64}
    @test jacdet(m4, 3.0) == 1.0

    @test mapsize(AffineMap(rand(3),rand(3))) == (3,)
    @test mapsize(AffineMap(LinearAlgebra.I,2)) == ()

    # Test construction and conversion
    @test AffineMap(1, 2*one(T)) isa ScalarAffineMap{T}
    @test AffineMap{BigFloat}(1, 2.0) isa ScalarAffineMap{BigFloat}
    @test AffineMap(rand(T,2),rand(T,2)) isa GenericAffineMap{T}
    @test AffineMap(rand(T,2,2),rand(T,2)) isa VectorAffineMap{T}
    @test AffineMap{Vector{T}}(rand(Int,2,2),rand(Int,2)) isa VectorAffineMap{T}
    @test AffineMap{SVector{2,T}}(rand(T,2,2),rand(T,2)) isa StaticAffineMap{T,2,2,4}
    @test AffineMap{SVector{2,T}}(rand(Int,2,2),rand(T,2)) isa StaticAffineMap{T,2,2,4}
    @test AffineMap{SVector{2,T}}(rand(Int,2,2),rand(Int,2)) isa StaticAffineMap{T,2,2,4}
    @test AffineMap{SVector{2,T}}(1,rand(Int,2)) isa GenericAffineMap{SVector{2,T},T,SVector{2,T}}
    @test AffineMap{SVector{2,T}}(SMatrix{3,2}(3,2,1,4,5,6),SVector(1,2,3)) isa StaticAffineMap{T,2,3,6}

    @test GenericAffineMap(SMatrix{2,2,T}(1,2,3,4), [1,2]) isa GenericAffineMap{SVector{2,T},SMatrix{2,2,T,4}}
    @test GenericAffineMap(SMatrix{1,1,Int}(1),SVector{1,T}(1)) isa GenericAffineMap{SVector{1,T}}
    @test GenericAffineMap(rand(T,2,2),rand(Int,2)) isa GenericAffineMap{Vector{T},Matrix{T},Vector{T}}
    @test GenericAffineMap{Vector{T}}(1, [1,2]) isa GenericAffineMap{Vector{T},T,Vector{Int}}
    @test GenericAffineMap{Vector{T}}(rand(Int,2,2), rand(Int,2)) isa GenericAffineMap{Vector{T},Matrix{T},Vector{T}}
    @test convert(GenericAffineMap{Vector{T}}, GenericAffineMap(rand(Int,2,2), rand(Int,2))) isa GenericAffineMap{Vector{T}}
    @test GenericAffineMap(rand(T,2),rand(T,2)) isa GenericAffineMap{T}
    @test GenericAffineMap{T}(rand(Int,2),rand(Int,2)) isa GenericAffineMap{T}

    @test VectorAffineMap(rand(T,2,2),rand(T,2)) isa VectorAffineMap{T}

    @test StaticAffineMap(SMatrix{1,1,T}(1),SVector{1,T}(1)) isa StaticAffineMap{T}
    @test StaticAffineMap(SMatrix{1,1,T}(1),SVector{1,Int}(1)) isa StaticAffineMap{T}

    @test StaticAffineMap(rand(T,2,2),SVector{2,T}(1,1)) isa StaticAffineMap{T}
    @test_throws AssertionError StaticAffineMap(rand(T,3,2),SVector{2,T}(1,1))
    @test StaticAffineMap(SMatrix{2,2,T}(1,2,3,4),[1,1]) isa StaticAffineMap{T,2}
    @test StaticAffineMap(SMatrix{3,2,T}(1,2,3,4,5,6),SVector{3,T}(1,2,3)) isa StaticAffineMap{T,2,3,6}

    @test convert(Map{SVector{2,T}}, AffineMap(rand(2,2),rand(2))) isa StaticAffineMap{T,2,2,4}
end

function test_abstractarraymap(T)
    A = rand(T, 3, 3)
    @test FunctionMaps.MapStyle(A) isa FunctionMaps.IsMap
    @test FunctionMaps.checkmap(A) == A
    @test Map(A) isa LinearMap
    @test Map(A) isa VectorLinearMap
    @test Map(SA[1 2; 3 4]) isa StaticLinearMap
    @test Map(Diagonal(rand(3))) isa GenericLinearMap
    @test mapsize(A) == size(A)
    @test applymap(A, 1:3) == A*(1:3)
    @test islinearmap(A)
    @test isaffinemap(A)
    @test affinematrix(A) == A
    @test affinevector(A) == [0,0,0]
    @test jacobian(A) isa ConstantMap
    @test jacobian(A, 1:3) == A
    @test inverse(A) ≈ inv(A)
    @test inverse(A, 1:3) ≈ A \ (1:3)
    @test canonicalmap(A) isa LinearMap
    @test FunctionMaps.equalmap(A) isa LinearMap
    @test isequalmap(A, LinearMap(A))
end
