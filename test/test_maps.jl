# test_maps.jl

function test_maps()
    @testset "$(rpad("Maps",80))" begin
        test_maps(Float64)
        test_maps(BigFloat)
        test_embedding_maps()
        test_composite_map(Float64)
    end
end

function suitable_point_to_map(m, n)
    x = @SVector ones(n)
end

# Test a map m with dimensions n
function test_generic_map(T, m, n)
    # Try to map a random vector
    isreal(T(0)) && (@test isreal(m))
    r = randvec(T,n)

    x = suitable_point_to_map(m, n)
    y1 = applymap(m, x)
    y2 = m * x
    @test y1 == y2
    y3 = m(x)
    @test y1 == y3

    mi = inv(m)
    xi1 = applymap(mi, y1)
    @test xi1 ≈ x
    xi2 = mi * y1
    @test xi2 ≈ x
    xi3 = m\y1
    @test xi3 ≈ x

    # if is_linear(m)
    #     x = suitable_point_to_map(m, n)
    #     y1 = applymap(m, x)
    #     x0 = @SVector zeros(T,n)
    #     a,b = linearize(m, x0)
    #     y2 = a*x+b
    #     @test y1 ≈ y2
    # end
end

randvec(T,n) = SVector{n,T}(rand(n))
randvec(T,m,n) = SMatrix{m,n,T}(rand(m,n))

function test_maps(T)
    a = T(0)
    b = T(1)
    c = T(2)
    d = T(3)
    m = interval_map(a, b, c, d)
    @test m(a) ≈ c
    @test m(b) ≈ d

    test_generic_map(T, m, 1)

    m2 = AffineMap(randvec(T, 2, 2), randvec(T, 2))
    test_generic_map(T, m2, 2)

    m3 = LinearMap(randvec(T, 2, 2))
    test_generic_map(T, m3, 2)

    # Test an affine map with a a scalar and b a vector
    m4 = AffineMap(T(1.2), randvec(T, 2))
    test_generic_map(T, m4, 2)

    m5 = AffineMap(randvec(T, 3, 3), randvec(T, 3))
    test_generic_map(T, m5, 3)

    m6 = m3∘m4
    @test typeof(m6) <: AffineMap
    test_generic_map(T, m6, 2)

    # Test special maps
    test_generic_map(T, scaling_map(T(2)), 1)

    test_generic_map(T, scaling_map(T(2)), 2)
    test_generic_map(T, scaling_map(T(2), T(3)), 2)
    test_generic_map(T, scaling_map(T(2), T(3), T(4)), 3)

    test_generic_map(T, IdentityMap{T}(), 1)
    test_generic_map(T, IdentityMap{SVector{2,T}}(), 2)

    # Diagonal maps

end

function test_embedding_maps()
    T1 = Float64
    T2 = Complex{Float64}
    T3 = SVector{1,Float64}
    T4 = SVector{2,Float64}
    T5 = SVector{2,Complex{Float64}}

    test_embedding_map(T1, T2)
    test_embedding_map(T3, T4)
    test_isomorphism_map(T2, T4)
end

function test_embedding_map(T1, T2)
    x = nonzero_element(T1)
    m = embedding_map(T2, T1)
    @test applymap(m, x) == convert_space(spacetype(T2), x)
    m2 = restriction_map(T1, T2)
    y = nonzero_element(T2)
    @test applymap(m2, y) == restrict_space(spacetype(T1), y)
end

function test_isomorphism_map(T1, T2)
    test_embedding_map(T1, T2)
    test_embedding_map(T2, T1)
    @test inv(isomorphism_map(T1,T2)) == isomorphism_map(T2, T1)
    x = nonzero_element(T1)
    m = isomorphism_map(T2, T1)
    y = applymap(m, x)
    @test apply_inverse(m, y) == x
    m2 = isomorphism_map(T1, T2)
    y2 = nonzero_element(T2)
    x2 = applymap(m2, y2)
    @test apply_inverse(m2, x2) == y2
end

function test_composite_map(T)
  a = T(0)
  b = T(1)
  c = T(2)
  d = T(3)
  ma = IdentityMap{T}()
  mb = interval_map(a, b, c, d)

  r = randvec(T, 2)
  m1 = ma∘mb
  test_generic_map(T, m1, 2)
  @test m1(r) ≈ ma(mb(r))
  m2 = m1∘mb
  test_generic_map(T, m2, 2)
  @test m2(r) ≈ m1(mb(r))
  m3 = mb∘m2
  test_generic_map(T, m3, 2)
  @test m3(r) ≈ mb(m2(r))
  m = m2∘m3
  test_generic_map(T, m, 2)
  @test m(r) ≈ m2(m3(r))

end
