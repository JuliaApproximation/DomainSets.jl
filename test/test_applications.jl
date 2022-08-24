
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

    hybrid_product = ProductDomain(["a", "b"], T(1.0)..T(2.0))
    @test @inferred(Random.gentype(hybrid_product)) == Tuple{String, T}
    @test typeof(rand(hybrid_product)) == Random.gentype(hybrid_product)
    @test @inferred(rand(hybrid_product)) in hybrid_product

    b = Ball(2.0, SA[T(1.0), T(2.0)])
    @test @inferred(Random.gentype(b)) == SVector{2, T}
    if T == BigFloat
        @test_broken typeof(rand(b)) == Random.gentype(b) # all rand tests for BigFloat StaticVector Balls are broken due to issue #114
    else
        @test typeof(rand(b)) == Random.gentype(b)
        @test @inferred(rand(b)) in b
        @test all(p in b for p in rand(b, 100))
        test_rng_consistency(b)
    end

    b = Ball(2.0, [T(1.0), T(2.0)])
    @test @inferred(Random.gentype(b)) == Vector{T}
    @test typeof(rand(b)) == Random.gentype(b)
    @test @inferred(rand(b)) in b
    @test all(p in b for p in rand(b, 100))
    test_rng_consistency(b)

    b = Ball(2.0, T(1.0))
    @test @inferred(Random.gentype(b)) == T
    @test_broken typeof(rand(b)) == Random.gentype(b) # broken due to issue #113
    @test_broken @inferred(rand(b)) in b # broken due to issue #113
    # test_rng_consistency(b) # uncomment when issue #113 is fixed

    # Higher dimension
    b = Ball(2.0, SA[T(1.0), T(1.0), T(1.0), T(1.0)])
    @test @inferred(Random.gentype(b)) == SVector{4, T}
    if T == BigFloat
        @test_broken typeof(rand(b)) == Random.gentype(b) # all rand tests for BigFloat StaticVector Balls are broken due to issue #114
    else
        @test typeof(rand(b)) == Random.gentype(b)
        @test @inferred(rand(b)) in b
        @test all(p in b for p in rand(b, 100))
        test_rng_consistency(b)
    end

    if T == Float64
        # Test numerical accuracy - two rectangles of the same size should have the same number of points
        rng = StableRNG(1)
        n = 1_000_000
        region_1 = Rectangle(zeros(4), [0.3, -0.3, 0.3, -0.3])
        lower_corner = fill(-0.5, 4)
        region_2 = Rectangle(lower_corner, lower_corner.+0.3)
        rs = rand(rng, b, n)
        @show n_1 = sum(r in region_1 for r in rs)
        @show n_2 = sum(r in region_2 for r in rs)
        @test isapprox(n_1/n, n_2/n, atol=0.01) 
    end
end

function test_rng_consistency(set)
    rng1 = Random.MersenneTwister(1)
    rng2 = Random.MersenneTwister(1)
    @test rand(rng1, set) == rand(rng2, set)
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
