function generic_map_tests(T)
    for map in maps_to_test(T)
        @test prectype(map) == T
        test_generic_map(map)
    end
end

# generic map test suite
function test_generic_map(m)
    @test convert(Map{domaintype(m)}, m) == m

    x = suitable_point_to_map(m)
    @test applymap(m, x) == m(x)
    y = m(x)
    @test isreal(m) == (isreal(x) && isreal(y))

    S = domaintype(m)
    U = codomaintype(m)
    @test x isa S
    @test y isa U

    if isaffinemap(m) && !isconstantmap(m) && !(prectype(m) == BigFloat)
        test_generic_inverse(m)
    else
        try
            # The map may not support an inverse, let's try
            test_generic_inverse(m)
        catch e
            # make sure it failed because the inverse is not defined
            @test e isa MethodError || prectype(m) == BigFloat
        end
    end

    if isaffinemap(m)
        M = affinematrix(m)
        @test size(M) == mapsize(m)
        test_generic_jacobian(m)
    else
        try # jacobian may not be implemented
            test_generic_jacobian(m)
        catch e
        end
    end

    if domaintype(m) == Float64
        @test convert(Map{BigFloat}, m) isa Map{BigFloat}
        @test convert(Map{BigFloat}, m) == m
    end
    if prectype(m) == Float64
        U = widertype(domaintype(m))
        @test convert(Map{U}, m) isa Map{U}
        @test convert(Map{U}, m) == m
    end

    if isaffinemap(m)
        A = affinematrix(m)
        b = affinevector(m)
        @test m(x) == A*x+b
    end

    @test hash(m) == map_hash(m)
end

function test_generic_inverse(m)
    M = mapsize(m,1)
    N = mapsize(m,2)
    x = suitable_point_to_map(m)
    y = m(x)

    if M==N
        inverse(m, y)      # trigger exception outside test if not implemented

        minv = inverse(m)
        @test minv(y) ≈ x
        @test inverse(m)(y) ≈ x
        @test inverse(m, y) ≈ x
        @test m\y ≈ x
        @test LazyInverse(m)(y) ≈ inverse(m, y)
    end
    if M >= N && numtype(m)!=BigFloat
        leftinverse(m, y)      # trigger exception outside test if not implemented

        mli = leftinverse(m)
        @test mli(y) ≈ x
        @test leftinverse(m, y) ≈ x
    end
    if M <= N && numtype(m)!=BigFloat
        rightinverse(m, y)      # trigger exception outside test if not implemented

        mri = rightinverse(m)
        @test m(mri(y)) ≈ y
        @test m(rightinverse(m, y)) ≈ y
    end
end

function test_generic_jacobian(m)
    x = suitable_point_to_map(m)
    # Trigger exception outside of test if jacobian(m, x) is not implemented.
    # test_generic_jacobian(m) should be called within a try/catch block
    jacobian(m, x)

    δ = sqrt(eps(prectype(m)))
    x2 = x .+ δ
    # we intentionally test jacobian(m, x) before testing jacobian(m)
    if !(m isa ProductMap)
        @test norm(m(x2) - m(x) + jacobian(m, x)*(x-x2)) < 100δ
    end
    jac = jacobian(m)
    @test jac(x) ≈ jacobian(m, x)
    j = jacobian(m, x)
    @test size(j) == mapsize(m)
    if j isa AbstractArray{Float64}
        y = similar(j)
        jacobian!(y, m, x)
        @test y ≈ jac(x)
    end
    if issquarematrix(jac(x))
        @test jacdet(m, x) ≈ det(jacobian(m, x))
    end
    if mapsize(m) == ()
        d = diffvolume(m, x)
        s = sqrt(det(transpose(jacobian(m,x))*jacobian(m,x)))
        @test d ≈ s || d ≈ -s
    else
        @test abs(diffvolume(m, x) - sqrt(det(transpose(jacobian(m,x))*jacobian(m,x)))) < 1e-7
    end
    if isaffinemap(m)
        @test diffvolume(m) isa ConstantMap
    end
end

function test_generic_functionality()
    m = LinearMap(2)
    @test convert_domaintype(Float64, m) isa ScalarLinearMap{Float64}
    @test convert_domaintype(Complex{Float64}, m) isa ScalarLinearMap{Complex{Float64}}
    @test convert_domaintype(SVector{2,BigFloat}, m) isa GenericLinearMap{SVector{2, BigFloat}, BigFloat}
    @test convert_domaintype(Any, m) === m
    @test convert_codomaintype(Float64, m) isa ScalarLinearMap{Float64}
    @test convert_codomaintype(Complex{Float64}, m) isa ScalarLinearMap{Complex{Float64}}
    @test_throws ArgumentError convert_codomaintype(SVector{2,BigFloat}, m)

    m2 = LinearMap{SVector{2,Float64}}(2)
    @test convert_domaintype(SVector{2,BigFloat}, m2) isa GenericLinearMap{SVector{2, BigFloat}, BigFloat}
    @test convert_codomaintype(SVector{2,BigFloat}, m2) isa GenericLinearMap{SVector{2, BigFloat}, BigFloat}
end
