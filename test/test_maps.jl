
randvec(T,n) = SVector{n,T}(rand(n))
randvec(T,m,n) = SMatrix{m,n,T}(rand(m,n))

suitable_point_to_map(m::Map) = suitable_point_to_map(m, domaintype(m))

suitable_point_to_map(m::Map, ::Type{SVector{N,T}}) where {N,T} = SVector{N,T}(rand(N))
suitable_point_to_map(m::Map, ::Type{T}) where {T<:Number} = rand(T)
suitable_point_to_map(m::Map, ::Type{<:AbstractVector{T}}) where {T} = rand(T, dimension(m))

suitable_point_to_map(m::DomainSets.ProductMap) =
    map(suitable_point_to_map, elements(m))

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

# Test a map m with dimensions n
function test_generic_map(T, m)
    # Try to map a random vector
    isreal(zero(T)) && (@test isreal(m))

    @test convert(Map{domaintype(m)}, m) == m

    x = suitable_point_to_map(m)
    y1 = applymap(m, x)
    y2 = m(x)
    @test y1 == y2

    x = suitable_point_to_map(m)
    S = domaintype(m)
    U = codomaintype(m)
    @test x isa S
    @test m(x) isa U

    try # try because the inverse may not be defined
        x = suitable_point_to_map(m)
        y = applymap(m, x)
        mi = inv(m)
        xi1 = applymap(mi, y)
        @test xi1 ≈ x
        xi2 = mi(y)
        @test xi2 ≈ x
        xi3 = m\y
        @test xi3 ≈ x
    catch
    end

    try # try because jacobian may not be implemented
        jac = jacobian(m)
        x = suitable_point_to_map(m)
        @test jac(x) == jacobian(m, x)
        if issquarematrix(jac(x))
            @test jacdet(m, x) == det(jacobian(m, x))
        end
        δ = sqrt(eps(T))
        x2 = x .+ δ
        if !(m isa ProductMap)
            @test norm(m(x2) .- (m(x)+jac(x)*(x-x2))) < 100δ
        end
    catch
    end

    if domaintype(m) == Float64
        @test convert(Map{BigFloat}, m) isa Map{BigFloat}
        @test convert(Map{BigFloat}, m) == m
    end
    if eltype(T) == Float64
        U = widertype(domaintype(m))
        @test convert(Map{U}, m) isa Map{U}
        @test convert(Map{U}, m) == m
    end

    if isaffine(m)
        A = matrix(m)
        b = vector(m)
        @test m(x) == A*x+b
    end

    # leftinv and rightinv tests currently fail on pinv for BigFloat
    if !(T == BigFloat)
        try # try because the inverse may not be defined
            mli = leftinv(m)
            @test mli(m(x)) ≈ x
            mri = rightinv(m)
            @test m(mri(x)) ≈ x
        catch
        end
    end
end

function test_maps(T)
    generic_tests(T)

    # Test special maps
    test_affine_maps(T)
    test_composite_map(T)
    test_product_map(T)
    test_wrapped_maps(T)
    test_scaling_maps(T)
    test_identity_map(T)
    test_rotation_map(T)
    test_cart_polar_map(T)
end

function generic_tests(T)
    maps = [
        IdentityMap{T}(),
        VectorIdentityMap{T}(10),
        ConstantMap{T}(one(T)),
        ConstantMap{T}(SVector{2,T}(1,2)),
        ZeroMap{T}(),
        UnityMap{T}(),
        interval_map(T(0), T(1), T(2), T(3)),
        AffineMap(randvec(T, 2, 2), randvec(T, 2)),
        LinearMap(randvec(T, 2, 2)),
        AffineMap(T(1.2), randvec(T, 2)),
        AffineMap(randvec(T, 3, 3), randvec(T, 3)),
        LinearMap(randvec(T, 2, 2)) ∘ AffineMap(T(1.2), randvec(T, 2))
    ]
    for map in maps
        test_generic_map(T, map)
    end
end

function test_affine_maps(T)
    test_linearmap(T)
    test_translation(T)
    test_affinemap(T)
end

function test_linearmap(T)
    m1 = LinearMap(2one(T))
    @test m1 isa LinearMap{T}
    @test domaintype(m1) == T
    @test islinear(m1)
    @test isaffine(m1)
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

    m2 = LinearMap(2)
    @test domaintype(m2) == Int
    @test m2(one(T)) isa T
    @test jacobian(m2, 1) == 2
    @test jacobian(m2) isa ConstantMap{Int}
    @test jacobian(m2, 1) == 2

    m3 = LinearMap(SMatrix{2,2}(one(T), 2one(T), 3one(T), 4one(T)))
    @test m3 isa LinearMap{SVector{2,T}}
    @test m3(SVector(1,2)) == SVector(7, 10)
    @test m3(SVector{2,T}(1,2)) == SVector{2,T}(7, 10)

    A = rand(T,2,2)
    m4 = LinearMap(A)
    @test m4 isa LinearMap{Vector{T}}
    @test m4([1,2]) ==  A * [1,2]
    @test jacobian(m4, [1,2]) == A

    # Test conversions
    @test LinearMap{T}(1) isa LinearMap{T,T}
    @test LinearMap{Vector{T}}(rand(Int,5,5)) isa LinearMap{Vector{T},Matrix{T}}
    @test LinearMap{SVector{2,T}}(rand(Int,2,2)) isa LinearMap{SVector{2,T},<:SMatrix{2,2,T}}
    @test LinearMap{SVector{2,T}}(2) isa LinearMap{SVector{2,T},T}

    @test convert(Map{SVector{2,T}}, m4) isa LinearMap{SVector{2,T},<:SMatrix{2,2,T}}
end

function test_translation(T)
    v = randvec(T,3)
    m = Translation(v)
    test_generic_map(T, m)
    @test !islinear(m)
    @test isaffine(m)
    @test inv(inv(m)) == m
    @test jacobian(m) isa ConstantMap
    @test vector(jacobian(m)) == [1 0 0; 0 1 0; 0 0 1]

    # Test conversions
    @test Translation{T}(1) isa Translation{T,T}
    @test Translation{SVector{2,T}}(rand(Int,2)) isa Translation{SVector{2,T},SVector{2,T}}
    @test Translation{Vector{T}}(rand(Int,2)) isa Translation{Vector{T},Vector{T}}
end

function test_affinemap(T)
    m1 = AffineMap(T(2), T(3))
    @test m1 isa AffineMap{T}
    @test !islinear(m1)
    @test isaffine(m1)
    @test m1(2) == 7

    m2 = AffineMap(T(2), SVector{2,T}(1,2))
    @test m2 isa AffineMap{SVector{2,T}}
    @test m2(SVector(1,2)) == SVector(3,6)

    @test m1 ∘ m1 isa AffineMap

    # Test conversions
    @test AffineMap(1, 2*one(T)) isa AffineMap{T,T,T}
    @test AffineMap{BigFloat}(1, 2.0) isa AffineMap{BigFloat,BigFloat,BigFloat}
    @test AffineMap(rand(T,2,2),rand(T,2)) isa AffineMap{Vector{T},Matrix{T},Vector{T}}
    @test AffineMap{Vector{T}}(rand(Int,2,2),rand(Int,2)) isa AffineMap{Vector{T},Matrix{T},Vector{T}}
    @test AffineMap{SVector{2,T}}(rand(T,2,2),rand(T,2)) isa AffineMap{SVector{2,T},<:SMatrix{2,2,T},SVector{2,T}}
    @test AffineMap{SVector{2,T}}(rand(Int,2,2),rand(T,2)) isa AffineMap{SVector{2,T},<:SMatrix{2,2,T},SVector{2,T}}
    @test AffineMap{SVector{2,T}}(rand(Int,2,2),rand(Int,2)) isa AffineMap{SVector{2,T},<:SMatrix{2,2,T},SVector{2,T}}
    @test AffineMap{SVector{2,T}}(1,rand(Int,2)) isa AffineMap{SVector{2,T},T,SVector{2,T}}
    @test AffineMap{SVector{2,T}}(rand(Int,3,2),SVector(1,2,3)) isa AffineMap{SVector{2,T},<:SMatrix{3,2,T},SVector{3,T}}

    @test convert(Map{SVector{2,T}}, AffineMap(rand(2,2),rand(2))) isa Map{SVector{2,T}}
end

function test_scaling_maps(T)
    test_generic_map(T, scaling_map(T(2)))
    test_generic_map(T, scaling_map(T(2), T(3)))
    test_generic_map(T, scaling_map(T(2), T(3), T(4)))
    test_generic_map(T, scaling_map(T(2), T(3), T(4), T(5)))
end

function test_identity_map(T)
    i1 = IdentityMap{T}()
    i2 = IdentityMap{SVector{2,T}}()
    test_generic_map(T, i1)
    test_generic_map(T, i2)
    @test i1 == i2
    @test islinear(i1)
    @test isaffine(i1)
    @test convert(IdentityMap{SVector{2,T}}, i1) === i2
    @test jacobian(i1) isa ConstantMap
    @test jacobian(i1, 1) == 1
    @test jacdet(i1, 1) == 1
    m1 = convert(DomainSets.AbstractAffineMap{T}, i1)
    @test m1 isa LinearMap{T}
    @test jacdet(m1, 1) == 1
    m2 = convert(DomainSets.LinearMap{T}, i1)
    @test m2 isa LinearMap{T}
    @test jacdet(m2, 1) == 1

    i3 = VectorIdentityMap{T}(10)
    test_generic_map(T, i3)
    r = rand(T, 10)
    @test i3(r) ≈ r
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

function test_wrapped_maps(T)
    m1 = WrappedMap{T}(cos)
    m2 = WrappedMap{T}(sin)
    @test m1(one(T)) ≈ cos(one(T))
    @test m2(one(T)) ≈ sin(one(T))
    m3 = m1 ∘ m2
    @test m3(one(T)) ≈ cos(sin(one(T)))
end

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
    test_generic_map(T, m2)
    m3 = rotation_map(phi, theta, psi)
    test_generic_map(T, m3)

    r = suitable_point_to_map(m2)
    @test norm(m2(r))≈norm(r)

    r = suitable_point_to_map(m3)
    @test norm(m3(r))≈norm(r)
    @test islinear(m3)
end

function test_cart_polar_map(T)
    m1 = CartToPolarMap{T}()
    test_generic_map(T, m1)
    @test !islinear(m1)

    m2 = PolarToCartMap{T}()
    test_generic_map(T, m2)
    @test !islinear(m2)
end


Base.isapprox(a::NTuple{L,SVector{N,T}}, b::NTuple{L,SVector{N,T}}) where {L,N,T} = compare_tuple(a,b)
Base.isapprox(a::NTuple{L,T}, b::NTuple{L,T}) where {L,T} = compare_tuple(a,b)

# Compare two iterable sequences for element-wise equality
compare_tuple(a, b) = reduce(&, map(isapprox, a, b))


@testset "maps" begin
    test_maps(Float64)
    test_maps(BigFloat)
end
