module DomainSetsMakieTests
using DomainSets, StaticArrays, Test
import Makie
import Makie: plot, Poly, Lines
using DomainSets: Sphere

@testset "Plotting" begin
    @testset "2D" begin
        @test plot((0..1) × (1..2)).plot isa Poly
        @test plot(UnitCircle()).plot isa Lines
        @test plot(UnitDisk()).plot isa Poly
        @test plot(Sphere(2.0, SVector(1.0, 0.5))).plot isa Lines
        @test plot(Sphere(3.0, [1.0, 0.5])).plot isa Lines
        Point(SVector(0.1,0.2))
    end
end
end # module