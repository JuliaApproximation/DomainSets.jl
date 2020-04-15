
function suitable_point_to_map(m)
    T = domaintype(m)
    if T <: SVector
        randvec(eltype(T),length(T))
    else
        x = T(rand())
    end
end

function suitable_point_to_map(m::DomainSets.ProductMap)
    x = ()
    for map in elements(m)
        x = (x..., suitable_point_to_map(map))
    end
    x
end

suitable_point_to_map(::CartToPolarMap{T}) where {T} = randvec(T,2)
suitable_point_to_map(::PolarToCartMap{T}) where {T} = randvec(T,2)
# Test a map m with dimensions n
function test_generic_map(T, m)
    # Try to map a random vector
    isreal(T(0)) && (@test isreal(m))
    x = suitable_point_to_map(m)
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
    xi4 = apply_inverse(m, y1)
    @test xi4 ≈ x

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

    test_generic_map(T, m)

    m2 = AffineMap(randvec(T, 2, 2), randvec(T, 2))
    test_generic_map(T, m2)

    m3 = LinearMap(randvec(T, 2, 2))
    test_generic_map(T, m3)
    r = randvec(T,2,1)
    @test jacobian(m3, r)*r ≈ m3*r

    # Test an affine map with a a scalar and b a vector
    m4 = AffineMap(T(1.2), randvec(T, 2))
    test_generic_map(T, m4)

    m5 = AffineMap(randvec(T, 3, 3), randvec(T, 3))
    test_generic_map(T, m5)

    m6 = m3∘m4
    @test typeof(m6) <: AffineMap
    test_generic_map(T, m6)

    # Test special maps
    test_composite_map(T)

    test_product_map(T)

    test_scaling_maps(T)

    test_identity_map(T)

    test_rotation_map(T)

    test_translation_map(T)

    test_cart_polar_map(T)

    # Diagonal maps

end

function test_scaling_maps(T)
    test_generic_map(T, scaling_map(T(2)))

    test_generic_map(T, scaling_map(T(2)))
    test_generic_map(T, scaling_map(T(2), T(3)))
    test_generic_map(T, scaling_map(T(2), T(3), T(4)))
    test_generic_map(T, scaling_map(T(2), T(3), T(4), T(5)))
end

function test_identity_map(T)
    test_generic_map(T, IdentityMap{T}())
    test_generic_map(T, IdentityMap{SVector{2,T}}())
    @test islinear(IdentityMap{T}())
end

function test_rotation_map(T)
    ϕ = T(pi)/4
    m = rotation_map(ϕ)
    x = [one(T), zero(T)]
    y = m*x
    @test y[1] ≈ sqrt(T(2))/2
    @test y[2] ≈ sqrt(T(2))/2

    ϕ = T(pi)/4
    m = rotation_map(ϕ, 0, 0)
    x = [zero(T), one(T), zero(T)]
    y = m*x
    @test y[1] ≈ 0
    @test y[2] ≈ sqrt(T(2))/2
    @test y[3] ≈ sqrt(T(2))/2

    # TODO: add more tests for a 3D rotation

    theta = T(rand())
    phi = T(rand())
    psi = T(rand())
    m2 = rotation_map(theta)
    test_generic_map(T, m2)
    m3 = rotation_map(phi, theta, psi)
    test_generic_map(T, m3)

    r = suitable_point_to_map(m2)
    @test norm(m2*r)≈norm(r)

    r = suitable_point_to_map(m3)
    @test norm(m3*r)≈norm(r)
    @test islinear(m3)
end

function test_translation_map(T)
    v = randvec(T,3)
    m = translation_map(v)
    test_generic_map(T, m)
    @test islinear(m)
end

function test_cart_polar_map(T)
    m1 = CartToPolarMap{T}()
    test_generic_map(T, m1)
    @test !islinear(m1)

    m2 = PolarToCartMap{T}()
    test_generic_map(T, m2)
    @test !islinear(m2)
end

function test_composite_map(T)
    a = T(0)
    b = T(1)
    c = T(2)
    d = T(3)
    ma = IdentityMap{T}()
    mb = interval_map(a, b, c, d)

    r = suitable_point_to_map(ma)
    m1 = ma∘mb
    test_generic_map(T, m1)
    @test m1(r) ≈ ma(mb(r))
    m2 = m1∘mb
    test_generic_map(T, m2)
    @test m2(r) ≈ m1(mb(r))
    m3 = mb∘m2
    test_generic_map(T, m3)
    @test m3(r) ≈ mb(m2(r))
    m = m2∘m3
    test_generic_map(T, m)
    @test m(r) ≈ m2(m3(r))

    @test !(typeof(CompositeMap(ma)) <: CompositeMap)
end

function test_product_map(T)
    a = T(0)
    b = T(1)
    c = T(2)
    d = T(3)
    ma = IdentityMap{T}()
    mb = interval_map(a, b, c, d)

    r1 = suitable_point_to_map(ma)
    r2 = suitable_point_to_map(ma)
    r3 = suitable_point_to_map(ma)
    r4 = suitable_point_to_map(ma)
    r5 = suitable_point_to_map(ma)

    m1 = tensorproduct(ma,mb)
    test_generic_map(T, m1)
    @test compare_tuple(m1((r1,r2)), (ma(r1),mb(r2)))
    m2 = tensorproduct(m1,mb)
    test_generic_map(T, m2)
    @test compare_tuple(m2((r1,r2,r3)), (ma(r1),mb(r2),mb(r3)) )
    m3 = tensorproduct(mb,m2)
    test_generic_map(T, m3)
    @test compare_tuple(m3((r1,r2,r3,r4)),(mb(r1),ma(r2),mb(r3),mb(r4)))
    m = tensorproduct(m1,m2)
    test_generic_map(T, m)
    @test compare_tuple(m((r1,r2,r3,r4,r5)),(m1((r1,r2))...,m2((r3,r4,r5))...))
end

Base.isapprox(a::NTuple{L,SVector{N,T}}, b::NTuple{L,SVector{N,T}}) where {L,N,T} = compare_tuple(a,b)
Base.isapprox(a::NTuple{L,T}, b::NTuple{L,T}) where {L,T} = compare_tuple(a,b)

# Compare two iterable sequences for element-wise equality
compare_tuple(a, b) = reduce(&, map(isapprox, a, b))


@testset "maps" begin
    test_maps(Float64)
    test_maps(BigFloat)
end
