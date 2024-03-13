@testset "map interval" begin
        @test interval_map(1.0, Inf, 2.0, Inf) == AffineMap(1.0, 1.0)
        @test interval_map(1.0, Inf, Inf, 2.0) == AffineMap(-1.0, 3.0)
        @test interval_map(Inf, 1.0, 2.0, Inf) == AffineMap(-1.0, 3.0)
        @test interval_map(Inf, 1.0, Inf, 2.0) == AffineMap(1.0, 1.0)
        @test interval_map(-Inf, Inf, -Inf, Inf) == IdentityMap()
        @test interval_map(-Inf, Inf, Inf, -Inf) == LinearMap(-1)
        @test interval_map(Inf, -Inf, -Inf, Inf) == LinearMap(-1)
        @test interval_map(Inf, -Inf, Inf, -Inf) == IdentityMap()
        @test interval_map(Inf, Inf, Inf, Inf) == IdentityMap()
        @test interval_map(-Inf, -Inf, -Inf, -Inf) == IdentityMap()
        @test interval_map(Inf, Inf, -Inf, -Inf) == LinearMap(-1)
        @test interval_map(-Inf, -Inf, Inf, Inf) == LinearMap(-1)
        @test_throws ArgumentError interval_map(-Inf, Inf, -Inf, -Inf)
        @test_throws ArgumentError interval_map(-1, 1, -1, Inf)
        @test_throws ArgumentError interval_map(-1, 1, -1, Inf)

        a = 1.0; b = 2.0; c = 3.0; d = 5;
        @test interval_map(a, b, c, d)(a) == c
        @test interval_map(a, b, c, d)(b) == d
        @test interval_map(a, b, c, d) == bounded_interval_map(a, b, c, d)
end
