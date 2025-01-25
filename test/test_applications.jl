
# these functions are copied from the FunctionMaps.jl tests:
# --- start of copy
randvec(T,n) = SVector{n,T}(rand(n))
randvec(T,m,n) = SMatrix{m,n,T}(rand(m,n))
suitable_point_to_map(m::Map) = suitable_point_to_map(m, domaintype(m))
suitable_point_to_map(m::Map, ::Type{SVector{N,T}}) where {N,T} = SVector{N,T}(rand(N))
suitable_point_to_map(m::Map, ::Type{T}) where {T<:Number} = rand(T)
suitable_point_to_map(m::Map, ::Type{<:AbstractVector{T}}) where {T} = rand(T, mapsize(m,2))
suitable_point_to_map(m::FunctionMaps.ProductMap) =
    map(suitable_point_to_map, components(m))
suitable_point_to_map(m::FunctionMaps.VcatMap{T,M,N}) where {T,M,N} =
    SVector{N,T}(rand(T,N))
# --- end of copy

suitable_point_to_map(::CartToPolarMap{T}) where {T} = randvec(T,2)
suitable_point_to_map(::PolarToCartMap{T}) where {T} = randvec(T,2)

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
    # FunctionMapsTests.test_generic_map(m2)
    m3 = rotation_map(phi, theta, psi)
    # FunctionMapsTests.test_generic_map(m3)

    r = suitable_point_to_map(m2)
    @test norm(m2(r))≈norm(r)

    r = suitable_point_to_map(m3)
    @test norm(m3(r))≈norm(r)
    @test islinearmap(m3)
end

function test_cart_polar_map(T)
    m1 = CartToPolarMap{T}()
    # FunctionMapsTests.test_generic_map(m1)
    @test !islinearmap(m1)
    @test isrealmap(m1)

    m2 = PolarToCartMap{T}()
    # FunctionMapsTests.test_generic_map(m2)
    @test !islinearmap(m2)
    @test isrealmap(m2)

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
    @test typeof(rand(b)) == Random.gentype(b)
    @test @inferred(rand(b)) in b
    @test all(p in b for p in rand(b, 100))
    test_rng_consistency(b)

    b = Ball(2.0, [T(1.0), T(2.0)])
    @test @inferred(Random.gentype(b)) == Vector{T}
    @test typeof(rand(b)) == Random.gentype(b)
    @test @inferred(rand(b)) in b
    @test all(p in b for p in rand(b, 100))
    test_rng_consistency(b)

    b = Ball(2.0, T(1.0))
    @test @inferred(Random.gentype(b)) == T
    @test typeof(rand(b)) == Random.gentype(b)
    @test @inferred(rand(b)) in b
    test_rng_consistency(b)

    # Higher dimension
    # Only works for Float64 since there is no randn(BigFloat)
    if T == Float64
        b = Ball(1.5, SA[1.0, -1.0, 2.0, -3.0])
        @test @inferred(Random.gentype(b)) == SVector{4, T}
        @test typeof(rand(b)) == Random.gentype(b)
        @test @inferred(rand(b)) in b
        @test all(p in b for p in rand(b, 100))
        test_rng_consistency(b)

        # Test numerical accuracy - two rectangles of the same size should have the same number of points
        rng = StableRNG(1)
        n = 1_000_000
        region_1 = Rectangle([0.0, -0.3, 0.0, -0.3], [0.3, 0.0, 0.3, 0.0]) .+ center(b)
        lower_corner = radius(b)*fill(-0.5, 4)
        region_2 = Rectangle(lower_corner, lower_corner.+0.3) .+ center(b)
        rs = rand(rng, b, n)
        n_1 = sum(r in region_1 for r in rs)
        n_2 = sum(r in region_2 for r in rs)
        @test isapprox(n_1, n_2, rtol=0.1)
    end
end

function test_rng_consistency(set)
    rng1 = Random.MersenneTwister(1)
    rng2 = Random.MersenneTwister(1)
    @test rand(rng1, set) == rand(rng2, set)
end

function test_coordinates()
    p1 = DomainSets.EuclideanSpherePoint(SA[1/sqrt(3),1/sqrt(3),1/sqrt(3)])
    @test p1 ∈ UnitSphere(Val(3))
    p2 = DomainSets.EuclideanSpherePoint([1/sqrt(3),1/sqrt(3),1/sqrt(3)])
    @test p2 ∈ UnitSphere(3)
    p3 = DomainSets.SphericalCoordinate(0.4, 0.5)
    @test p3 ∈ UnitSphere(Val(3))
    @test DomainSets.SphericalCoordinate(1, 0.5) isa DomainSets.SphericalCoordinate{Float64}
    @test DomainSets.SphericalCoordinate(1, 2) isa DomainSets.SphericalCoordinate{Float64}
end

function test_applications(T)
    test_rotation_map(T)
    test_cart_polar_map(T)
    test_rand(T)
    test_coordinates()
end

@testset "applications" begin
    test_applications(Float64)
    test_applications(BigFloat)

    # Test an additional composite map
    m1 = LinearMap(SMatrix{2,2}(1,2,3,4.0))
    m2 = CartToPolarMap()
    cmap = m1 ∘ m2 ∘ m1
    @test cmap isa Map
    # FunctionMapsTests.test_generic_map(cmap)
end
