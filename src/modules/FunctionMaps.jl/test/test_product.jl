function test_product_map(T)
    ma = StaticIdentityMap{T}()
    mb = DomainSets.interval_map(T(0), T(1), T(2), T(3))

    r1 = suitable_point_to_map(ma)
    r2 = suitable_point_to_map(ma)
    r3 = suitable_point_to_map(ma)
    r4 = suitable_point_to_map(ma)
    r5 = suitable_point_to_map(ma)

    m1 = productmap(ma,mb)
    test_generic_map(m1)
    @test m1(SVector(r1,r2)) ≈ SVector(ma(r1),mb(r2))
    m2 = productmap(m1,mb)
    test_generic_map(m2)
    @test m2(SVector(r1,r2,r3)) ≈ SVector(ma(r1),mb(r2),mb(r3))
    m3 = productmap(mb,m2)
    test_generic_map(m3)
    @test m3(SVector(r1,r2,r3,r4)) ≈ SVector(mb(r1),ma(r2),mb(r3),mb(r4))
    m4 = productmap(m1,m2)
    test_generic_map(m4)
    @test m4(SVector(r1,r2,r3,r4,r5)) ≈ SVector(m1(SVector(r1,r2))...,m2(SVector(r3,r4,r5))...)

    @test ProductMap(ma, mb) == ProductMap([ma,mb])
    @test hash(ProductMap(ma, mb)) == hash(ProductMap([ma,mb]))

    m5 = productmap(AffineMap(SMatrix{2,2,T}(1.0,2,3,4), SVector{2,T}(1,3)), LinearMap{T}(2.0))
    @test domaintype(component(m5,1)) == SVector{2,T}
    @test domaintype(component(m5,2)) == T
    test_generic_map(m5)
    x = SVector{3,T}(rand(T,3))
    @test m5(x) ≈ SVector(component(m5,1)(SVector(x[1],x[2]))...,component(m5,2)(x[3]))

    m6 = ProductMap([ma,mb])
    @test m6 isa DomainSets.VectorProductMap
    @test mapsize(m6) == (2,2)
    @test convert(Map{SVector{2,T}}, m6) isa DomainSets.VcatMap
    test_generic_map(m6)

    d1 = ProductDomain(T(1)..T(2), T(1)..T(2))
    d2 = ProductDomain(T(2)..T(4), T(2)..T(4))
    @test mapto(d1, d2) isa DomainSets.VcatMap
    @test jacobian(mapto(d1,d2)) isa ConstantMap
    @test jacdet(mapto(d1,d2)) == ConstantMap{SVector{2,T}}(4)
    d1v = ProductDomain([T(1)..T(2), T(1)..T(2)])
    d2v = ProductDomain([T(2)..T(4), T(2)..T(4)])
    @test mapto(d1v, d2v) isa DomainSets.VectorProductMap
    @test jacobian(mapto(d1v,d2v)) isa ConstantMap
    @test jacdet(mapto(d1v,d2v)) == ConstantMap{Vector{T}}(4)

    @test components(productmap(LinearMap(1), LinearMap(one(T)))) isa Tuple{<:Map{T},<:Map{T}}

    # test a non-rectangular product map
    m_rect = AffineMap(SA[one(T) 2; 3 4; 5 6], SA[one(T),one(T),zero(T)])
    pm = productmap(m_rect, m_rect)
    @test pm isa DomainSets.VcatMap{T,6,4,(3,3),(2,2)}
    @test jacobian(pm, SA[one(T),0,0,0]) isa SMatrix{6,4,T}
    @test diffvolume(pm, SA[one(T),0,0,0]) == 24
end
