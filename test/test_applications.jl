
function test_rotation_map(T)
    ϕ = T(pi)/4
    m = rotation_map(ϕ)
    x = [one(T), zero(T)]
    y = m(x)
    @test y[1] ≈ sqrt(T(2))/2
    @test y[2] ≈ sqrt(T(2))/2

    ϕ = T(pi)/4
    m = rotation_map(ϕ, 0, 0)
    x = [zero(T), one(T), zero(T)]
    y = m(x)
    @test y[1] ≈ 0
    @test y[2] ≈ sqrt(T(2))/2
    @test y[3] ≈ sqrt(T(2))/2

    # TODO: add more tests for a 3D rotation

    theta = T(rand())
    phi = T(rand())
    psi = T(rand())
    m2 = rotation_map(theta)
    test_generic_map(m2)
    m3 = rotation_map(phi, theta, psi)
    test_generic_map(m3)

    r = suitable_point_to_map(m2)
    @test norm(m2(r))≈norm(r)

    r = suitable_point_to_map(m3)
    @test norm(m3(r))≈norm(r)
    @test islinear(m3)
end

function test_cart_polar_map(T)
    m1 = CartToPolarMap{T}()
    test_generic_map(m1)
    @test !islinear(m1)
    @test isreal(m1)

    m2 = PolarToCartMap{T}()
    test_generic_map(m2)
    @test !islinear(m2)
    @test isreal(m2)

    @test inverse(m1) == m2
    @test inverse(m2) == m1
end

function test_rand(T)
    r = Rectangle(T[-1, 2], T[3, 4])
    @test @inferred(Random.gentype(r)) == Vector{T}
    @test typeof(rand(r)) == Random.gentype(r)
    @test @inferred(rand(r)) in r

    r = Rectangle(SA[T(-1), T(2)], SA[T(3), T(4)])
    @test @inferred(Random.gentype(r)) == SVector{2, T}
    @test typeof(rand(r)) == Random.gentype(r)
    @test @inferred(rand(r)) in r

    hybrid_product = ProductDomain(["a", "b"], 1.0..2.0)
    @test @inferred(Random.gentype(hybrid_product)) == Tuple{String, Float64}
    @test typeof(rand(hybrid_product)) == Random.gentype(hybrid_product)
    @test @inferred(rand(hybrid_product)) in hybrid_product

    # b = Ball(2.0, SA[1.0, 2.0])
    # @test @inferred(Random.gentype(r)) == SVector{2, T}
    # @test typeof(rand(r)) == Random.gentype(r)
    # @test @inferred(rand(r)) in r
end

function test_applications(T)
    test_rotation_map(T)
    test_cart_polar_map(T)
    test_rand(T)
end

@testset "applications" begin
    test_applications(Float64)
    test_applications(BigFloat)
end
