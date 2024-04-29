module DomainSetsMakieExt
using DomainSets
using DomainSets.StaticArrays
import Makie
using DomainSets: leftendpoint, rightendpoint, Rectangle, HyperRectangle, Point, DomainPoint, pointval, point
import Makie: Point2f, Rect, Circle, Poly, Lines, convert_arguments, HyperSphere, Vec, Scatter, PointBased
import Base: convert

function convert(::Type{Makie.HyperRectangle}, r::Rectangle{<:SVector{N}}) where N
    l = leftendpoint(r)
    r = rightendpoint(r)
    Rect(convert(Vec{N}, l), convert(Vec{N}, r .- l))
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



Makie.plottype(a::HyperRectangle) = Poly
Makie.plottype(a::Sphere) = Lines
Makie.plottype(a::Ball) = Poly
Makie.plottype(a::Union{DomainPoint,Point}) = Scatter

end # module