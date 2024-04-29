module DomainSetsMakieExt
using DomainSets
using DomainSets.StaticArrays
import Makie
using DomainSets: leftendpoint, rightendpoint, Rectangle, HyperRectangle, DomainPoint, pointval, point
import Makie: Point2f, Rect, Circle, Poly, Point, Lines, convert_arguments, HyperSphere, Vec, Scatter
import Base: convert

function convert(::Type{Vector{Point2f}}, r::Rectangle)
    (a,c) = leftendpoint(r)
    (b,d) = rightendpoint(r)
    Point2f[(a,c), (b,c), (b,d), (a,d)]
end

function convert(::Type{Makie.HyperRectangle}, r::Rectangle{<:SVector{N}}) where N
    l = leftendpoint(r)
    r = rightendpoint(r)
    Rect(convert(Vec{N}, l), convert(Vec{N}, r .- l))
end

convert(::Type{<:HyperSphere}, r::Union{Sphere{SVector{N,T}},Ball{SVector{N,T}}}) where {N,T} = HyperSphere{N,T}(Point(center(r)), radius(r))

function convert(::Type{<:HyperSphere}, r::Union{Sphere{<:AbstractVector{T}}, Ball{<:AbstractVector{T}}}) where T
    N = length(center(r))
    HyperSphere{N,T}(Point{N}(center(r)), radius(r))
end

convert(::Type{<:Makie.Point}, r::Point) = Makie.Point(pointval(r))
convert(::Type{<:Makie.Point}, r::DomainPoint) = convert(Makie.Point, point(r))

convert_arguments(::Type{Plt}, r::HyperRectangle) where Plt <: Poly = convert_arguments(Plt, convert(Makie.HyperRectangle, r))
convert_arguments(::Type{Plt}, r::Ball) where Plt <: Poly = convert_arguments(Plt, convert(HyperSphere, r))
convert_arguments(::Type{Plt}, r::Sphere) where Plt <: Lines = convert_arguments(Plt, convert(HyperSphere, r))
convert_arguments(::Type{Plt}, r::Union{DomainPoint,Point}) where Plt <: Scatter = convert_arguments(Plt, convert(Makie.Point, r))



Makie.plottype(a::HyperRectangle) = Poly
Makie.plottype(a::Sphere) = Lines
Makie.plottype(a::Ball) = Poly
Makie.plottype(a::Union{DomainPoint,Point}) = Scatter

end # module