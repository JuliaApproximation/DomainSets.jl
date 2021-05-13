@testset "canonical domains" begin
    @test canonicaldomain(2..3) isa ChebyshevInterval
    @test canonicaldomain(2..3) == ChebyshevInterval{Float64}()
    @test mapto(UnitInterval(), UnitInterval()) isa IdentityMap
    @test mapto(UnitInterval(), ChebyshevInterval()) isa AffineMap
end
