module DomainSetsMakieExt
using DomainSets
using DomainSets.StaticArrays
import Makie
using DomainSets: leftendpoint, rightendpoint, Rectangle, HyperRectangle, Point, DomainPoint, pointval, point
import Makie: Point2f, Point3f, Rect, Circle, Poly, Lines, Mesh, convert_arguments, HyperSphere, Vec, Scatter, PointBased
import Base: convert

# infimum and supremum are defined for every HyperRectangle, e.g. UnitSquare()
function convert(::Type{Makie.HyperRectangle}, r::HyperRectangle)
    N = dimension(r)
    l, u = infimum(r), supremum(r)
    Rect(convert(Vec{N}, l), convert(Vec{N}, u .- l))
end

convert(::Type{<:HyperSphere}, r::Union{Sphere{SVector{N,T}},Ball{SVector{N,T}}}) where {N,T} = HyperSphere{N,T}(Point(center(r)), radius(r))

function convert(::Type{<:HyperSphere}, r::Union{Sphere{<:AbstractVector{T}}, Ball{<:AbstractVector{T}}}) where T
    N = length(center(r))
    HyperSphere{N,T}(Makie.Point{N}(center(r)), radius(r))
end

convert(::Type{<:Makie.Point}, r::Point) = Makie.Point(pointval(r))
convert(::Type{<:Makie.Point}, r::DomainPoint) = convert(Makie.Point, point(r))

convert_arguments(::Type{Plt}, r::HyperRectangle; kwds...) where Plt <: Poly = convert_arguments(Plt, convert(Makie.HyperRectangle, r); kwds...)
convert_arguments(::Type{Plt}, r::Ball; kwds...) where Plt <: Poly = convert_arguments(Plt, convert(HyperSphere, r); kwds...)
convert_arguments(::Type{Plt}, r::Sphere; kwds...) where Plt <: Lines = convert_arguments(Plt, convert(HyperSphere, r); kwds...)
convert_arguments(::Type{Plt}, r::Point; kwds...) where Plt <: Scatter = convert_arguments(Plt, convert(Makie.Point, r); kwds...)

# the unit simplex is a triangle in 2D and a tetrahedron in 3D
function convert_arguments(::Type{Plt}, s::UnitSimplex; kwds...) where Plt <: Poly
    dimension(s) == 2 || throw(ArgumentError("Only 2D simplices are plotted as polygons, use a mesh in 3D"))
    convert_arguments(Plt, [Point2f(0, 0), Point2f(1, 0), Point2f(0, 1)]; kwds...)
end
function convert_arguments(::Type{Plt}, s::UnitSimplex; kwds...) where Plt <: Mesh
    dimension(s) == 3 || throw(ArgumentError("Only 3D simplices are plotted as meshes"))
    convert_arguments(Plt, [Point3f(0, 0, 0), Point3f(1, 0, 0), Point3f(0, 1, 0), Point3f(0, 0, 1)], [1 2 3; 1 2 4; 1 3 4; 2 3 4]; kwds...)
end

# a union is plotted by plotting each of its components in the same color
Makie.@recipe DomainUnionPlot (domain,) begin
    "Color of every component of the union."
    color = @inherit linecolor
    cycle = [:color]
end

function Makie.plot!(p::DomainUnionPlot)
    for d in components(p[:domain][])
        Makie.plot!(p, d; color = p[:color])
    end
    p
end



Makie.plottype(a::HyperRectangle) = Poly
Makie.plottype(a::Sphere) = Lines
Makie.plottype(a::Ball) = Poly
Makie.plottype(a::Union{DomainPoint,Point}) = Scatter
Makie.plottype(s::UnitSimplex) = dimension(s) == 3 ? Mesh : Poly
Makie.plottype(a::UnionDomain) = DomainUnionPlot

end # module