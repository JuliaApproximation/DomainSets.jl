
"""
    Point(x)

represents a single point at `x`.
"""
struct Point{T} <: Domain{T}
    x::T
end

similardomain(d::Point, ::Type{T}) where {T} = Point{T}(d.x)

convert(::Type{Number}, d::Point{<:Number}) = d.x
convert(::Type{N}, d::Point{<:Number}) where N<:Number =
    convert(N, convert(Number, d.x))
Number(d::Point) = convert(Number, d)

convert(::Type{Domain}, c::Number) = Point(c)
convert(::Type{Domain{T}}, c::Number) where T = Point{T}(c)

pointval(d::Point) = d.x

isequaldomain(d1::Point, d2::Point) = pointval(d1) == pointval(d2)
hash(d::Point, h::UInt) = hashrec("Point", pointval(d), h)

indomain(x, d::Point) = x == pointval(d)
isempty(::Point) = false

approx_indomain(x, d::Point, tolerance) = norm(x-pointval(d)) <= tolerance

dimension(d::Point{<:AbstractVector}) = length(pointval(d))

canonicaldomain(d::Point{T}) where {T<:StaticTypes} = Point(zero(T))
canonicaldomain(d::Point{T}) where {T<:AbstractVector} =
    Point(zeros(eltype(T),dimension(d)))

mapfrom_canonical(d::Point) = Translation(pointval(d))

isopenset(d::Point) = false
isclosedset(d::Point) = true

boundary(d::Point) = d
boundingbox(d::Point) = pointval(d)..pointval(d)

infimum(d::Point) = pointval(d)
supremum(d::Point) = pointval(d)

interior(d::Point) = emptyspace(d)
closure(d::Point) = d

choice(d::Point) = pointval(d)

distance_to(d::Point, x) = norm(x-pointval(d))

mapped_domain(invmap, p::Point) = Point(inverse(invmap, pointval(p)))
map_domain(map, p::Point) = Point(applymap(map, pointval(p)))
parametric_domain(map, p::Point) = Point(applymap(map, pointval(p)))

for op in (:+,:-)
    @eval $op(a::Point, b::Point) = Point($op(pointval(a),pointval(b)))
end


intersectdomain1(d1::Point, d2) = pointval(d1) ∈ d2 ? d1 : emptyspace(d1)
intersectdomain2(d1, d2::Point) = pointval(d2) ∈ d1 ? d2 : emptyspace(d2)

# Interval minus a point:
setdiffdomain(d::Interval, x::Number) = setdiffdomain(d, Point(x))
setdiffdomain(d::Interval, p::Point) = setdiffdomain(promote_domains((d,p))...)
function setdiffdomain(d::Interval{L,R,T}, p::Point{T}) where {L,R,T}
    a = leftendpoint(d)
    b = rightendpoint(d)
    x = pointval(p)

    a == x && return Interval{:open,R,T}(a,b)
    a < x < b && return UnionDomain(Interval{L,:open,T}(a,x), Interval{:open,R,T}(x,b))
    b == x && return Interval{L,:open,T}(a,b)
    return d
end

issubset1(d1::Point, d2) = pointval(d1) ∈ d2

setdiffdomain1(p::Point, d2) = issubset(p, d2) ? emptyspace(p) : p

intersectdomain(d1::Point, d2::Point) = pointval(d1) ∈ d2 ? d1 : emptyspace(d1)

show(io::IO,d::Point) = print(io,"Point(", pointval(d), ")")
