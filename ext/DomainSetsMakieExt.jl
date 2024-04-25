module DomainSetsMakieExt
using DomainSets, Makie
using DomainSets: leftendpoint, rightendpoint, Rectangle
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

Makie.convert_arguments(::Type{<:Poly}, r::Rectangle) = (convert(Rect, r),)

Makie.plottype(a::Rectangle) = Poly

end # module