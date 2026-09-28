module DomainSetsMakieTests
using DomainSets, StaticArrays, Test
import Makie
import Makie: plot, Poly, Lines, Scatter
using DomainSets: Sphere, ×

@testset "Plotting" begin
    @testset "2D" begin
        @test plot(cartesianproduct(0..1, 1..2)).plot isa Poly
        @test plot(UnitCircle()).plot isa Lines
        @test plot(UnitDisk()).plot isa Poly
        @test plot(Sphere(2.0, SVector(1.0, 0.5))).plot isa Lines
        @test plot(Sphere(3.0, [1.0, 0.5])).plot isa Lines
        @test plot(Point(SVector(0.1,0.2))).plot isa Scatter
    end

    @testset "3D" begin
        @test plot(cartesianproduct(0..1, 1..2, 3..4)).plot isa Poly
        @test plot(UnitBall()).plot isa Poly
        @test plot(Sphere(2.0, SVector(1.0, 0.5,0.5))).plot isa Lines
        @test plot(Sphere(3.0, [1.0, 0.5,0.5])).plot isa Lines
    end

    @testset "unit cubes" begin
        for d in (UnitSquare(), UnitCube())
            p = plot(d).plot
            @test p isa Poly
        end
        @test convert(Makie.HyperRectangle, UnitSquare()) == Makie.Rect(Makie.Vec(0.0, 0.0), Makie.Vec(1.0, 1.0))
        @test convert(Makie.HyperRectangle, UnitCube()) == Makie.Rect(Makie.Vec(0.0, 0.0, 0.0), Makie.Vec(1.0, 1.0, 1.0))
    end

    @testset "simplices" begin
        p = plot(UnitSimplex(Val(2))).plot
        @test p isa Poly
        @test p[1][] == Makie.Point2f.([(0,0), (1,0), (0,1)])
        @test plot(UnitSimplex(2)).plot isa Poly
        p = plot(UnitSimplex(Val(3))).plot
        @test p isa Makie.Mesh
        @test Makie.GeometryBasics.coordinates(p[1][]) == Makie.Point3f.([(0,0,0), (1,0,0), (0,1,0), (0,0,1)])
        @test length(Makie.GeometryBasics.faces(p[1][])) == 4
        @test_throws ArgumentError plot(UnitSimplex(Val(4)))
    end

    @testset "unions" begin
        color(c) = Makie.to_color(c[])
        u = UnionDomain(UnitDisk(), Ball(0.5, SVector(2.5, 0.0)))
        fig, ax, p = plot(u)
        @test length(p.plots) == 2
        @test all(c -> c isa Poly, p.plots)
        @test color(p.plots[1].color) == color(p.plots[2].color) == color(p.color)
        # colors cycle between unions (Makie < 0.24 only cycles between unions of the same type)
        @test color(Makie.plot!(ax, UnionDomain(UnitDisk(), Ball(0.5, SVector(2.5, 2.5)))).color) ≠ color(p.color)
        # the components may be different kinds of domains
        q = Makie.plot!(ax, UnionDomain(Ball(0.5, SVector(0.0, 2.5)), Rectangle(SVector(1.5,2.0), SVector(2.5,3.0))))
        @test all(c -> c isa Poly, q.plots)
        @test color(q.plots[1].color) == color(q.plots[2].color)
        r = Makie.plot!(ax, u; color=:red)
        @test all(c -> color(c.color) == Makie.to_color(:red), r.plots)

        # intervals are plotted by IntervalSets
        if Base.get_extension(DomainSets.IntervalSets, :IntervalSetsMakieExt) !== nothing
            @test length(plot(UnionDomain(0..1, 2..3)).plot.plots) == 2
        end
    end
end
end # module