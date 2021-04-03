
"""
    Point(x)

represents a single point at `x`.
"""
struct Point{T} <: Domain{T}
    x::T
end

similardomain(d::Point, ::Type{T}) where {T} = Point{T}(d.x)

convert(::Type{Number}, d::Point{<:Number}) = d.x
convert(::Type{N}, d::Point{<:Number}) where N<:Number = convert(N, convert(Number, d.x))
Number(d::Point) = convert(Number, d)

convert(::Type{Domain}, c::Number) = Point(c)
convert(::Type{Domain{T}}, c::Number) where T = Point{T}(c)

==(a::Point,b::Point) = a.x == b.x
indomain(x, d::Point) = x == d.x
isempty(::Point) = false

approx_indomain(x, d::Point, tolerance) = norm(x-d.x) <= tolerance

isopenset(d::Point) = false
isclosedset(d::Point) = true

boundary(d::Point) = d

interior(d::Point{T}) where {T} = EmptySpace{T}()
closure(d::Point) = d

point_in_domain(d::Point) = d.x

dimension(d::Point{Vector{T}}) where {T} = length(d.x)

function map_domain(map, p::Point)
    x = applymap(map, p.x)
    Point(x)
end

mapped_domain(invmap, p::Point) = map_domain(inv(invmap), p)

for op in (:+,:-)
    @eval $op(a::Point, b::Point) = Point($op(a.x,b.x))
end

# Interval minus a point:
setdiffdomain(d::Interval, x::Number) = setdiffdomain(d, Point(x))
setdiffdomain(d::Interval, p::Point) = setdiffdomain(promote_domains((d,p))...)
function setdiffdomain(d::Interval{L,R,T}, p::Point{T}) where {L,R,T}
    a = leftendpoint(d)
    b = rightendpoint(d)
    x = p.x

    a == x && return Interval{:open,R,T}(a,b)
    a < x < b && return UnionDomain(Interval{L,:open,T}(a,p.x), Interval{:open,R,T}(p.x,b))
    b == x && return Interval{L,:open,T}(a,b)
    return d
end

issubset(p::Point, d::Domain) = p.x ∈ d

setdiffdomain1(p::Point, d2) = issubset(p, d2) ? EmptySpace{eltype(p)}() : p
