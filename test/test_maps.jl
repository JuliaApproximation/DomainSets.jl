# test_maps.jl

function test_maps()
    @testset "$(rpad("Maps",80))" begin
        test_maps(Float64)
        test_maps(BigFloat)
    end
end

function suitable_point_to_map(m, n)
    x = @SVector ones(n)
end

# Test a map m with dimensions n
function test_generic_map(T, m, n)
    # Try to map a random vector
    x = suitable_point_to_map(m, n)
    y1 = applymap(m, x)
    y2 = m * x
    @test y1 == y2

    mi = inv(m)
    xi1 = applymap(mi, y1)
    @test xi1 ≈ x
    xi2 = mi * y1
    @test xi2 ≈ x

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
