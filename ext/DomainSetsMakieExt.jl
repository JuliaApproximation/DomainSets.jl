module DomainSetsMakieExt
using DomainSets
using DomainSets.StaticArrays
import Makie
using DomainSets: leftendpoint, rightendpoint, Rectangle
using Makie: Point2f, Rect, Circle, Poly, Point, Lines
import Base: convert

function convert(::Type{Vector{Point2f}}, r::Rectangle)
    (a,c) = leftendpoint(r)
    (b,d) = rightendpoint(r)
    Point2f[(a,c), (b,c), (b,d), (a,d)]
end

function convert(::Type{Rect}, r::Rectangle)
    (a,c) = leftendpoint(r)
    (b,d) = rightendpoint(r)
    Rect(a, c, b-a, d-c)
end

convert(::Type{Circle}, r::Sphere{<:SVector{2}}) = Circle(Point(center(r)), radius(r))
convert(::Type{Circle}, r::Ball{<:SVector{2}}) = Circle(Point(center(r)), radius(r))

Makie.convert_arguments(::Type{<:Poly}, r::Rectangle) = (convert(Rect, r),)
Makie.convert_arguments(::Type{<:Poly}, r::Ball{<:SVector{2}}) = (convert(Circle, r),)
Makie.convert_arguments(::Type{<:Lines}, r::Sphere{<:SVector{2}}) = (convert(Circle, r),)


Makie.plottype(a::Rectangle) = Poly
Makie.plottype(a::Sphere{<:SVector{2}}) = Lines
Makie.plottype(a::Ball{<:SVector{2}}) = Poly

end # module