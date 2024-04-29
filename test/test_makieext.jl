module DomainSetsMakieTests
using DomainSets, StaticArrays, Test
import Makie
import Makie: plot, Poly, Lines, Scatter
using DomainSets: Sphere

@testset "Plotting" begin
    @testset "2D" begin
        @test plot((0..1) × (1..2)).plot isa Poly
        @test plot(UnitCircle()).plot isa Lines
        @test plot(UnitDisk()).plot isa Poly
        @test plot(Sphere(2.0, SVector(1.0, 0.5))).plot isa Lines
        @test plot(Sphere(3.0, [1.0, 0.5])).plot isa Lines
        @test plot(Point(SVector(0.1,0.2))).plot isa Scatter
    end

    @testset "3D" begin
        @test plot((0..1) × (1..2) × (3..4)).plot isa Poly
        @test plot(UnitBall()).plot isa Poly
        @test plot(Sphere(2.0, SVector(1.0, 0.5,0.5))).plot isa Lines
        @test plot(Sphere(3.0, [1.0, 0.5,0.5])).plot isa Lines
    end
end
end # module