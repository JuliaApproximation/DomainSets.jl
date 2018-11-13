# simple.jl
# A collection of simple domains.

###########################
# The unit ball and sphere
###########################

"""
The unit ball (of radius 1) in `N` dimensions.
"""
struct UnitHyperBall{N,T} <: EuclideanDomain{N,T} end

UnitHyperBall{N}() where N = UnitHyperBall{N,Float64}()

const UnitDisk{T} = UnitHyperBall{2,T}
const UnitBall{T} = UnitHyperBall{3,T}

indomain(x, ::UnitHyperBall) = norm(x) <= 1

approx_indomain(x, ::UnitHyperBall, tolerance) = norm(x) <= 1+tolerance

isempty(::UnitHyperBall) = false

show(io::IO, d::UnitHyperBall{N}) where {N} = print(io, "the $(N)-dimensional unit ball")

# We choose the origin here
point_in_domain(d::UnitHyperBall) = zero(eltype(d))


###########################
# An n-dimensional simplex
###########################

struct UnitSimplex{N,T} <: EuclideanDomain{N,T} end
UnitSimplex{N}() where N = UnitSimplex{N,Float64}()

indomain(x, ::UnitSimplex) = mapreduce( t-> t >= 0, &, x) && norm(x,1) <= 1

approx_indomain(x, ::UnitSimplex, tolerance) = mapreduce( t-> t >= -tolerance, &, x) && norm(x,1) <= 1+tolerance

isempty(::UnitSimplex) = false

# We pick the origin, because it belongs to the domain regardless of what T is
point_in_domain(d::UnitSimplex) = zero(eltype(d))



##############
# A cylinder
##############

cylinder(::Type{T} = Float64) where {T} = UnitDisk{T}() × UnitInterval{T}()

cylinder(radius::T, length::T) where {T} = radius * UnitDisk{T}() × (0 .. length)

"""
    Point(x)

represents a single point at `x`.
"""
struct Point{T} <: Domain{T}
    x::T
end

convert(::Type{Number}, d::Point{<:Number}) = d.x
convert(::Type{N}, d::Point{<:Number}) where N<:Number = convert(N, convert(Number, d.x))
Number(d::Point) = convert(Number, d)

convert(::Type{Domain{T}}, d::Point{T}) where T = d
convert(::Type{Domain{T}}, d::Point) where T = Point(T(d.x))
convert(::Type{Domain}, c::Number) = Point(c)
convert(::Type{Domain{T}}, c::Number) where T = Point(convert(T,c))
convert(::Type{Domain}, s::Set) = UnionDomain(map(Domain,collect(s)))
Domain(d) = convert(Domain, d)

==(a::Point,b::Point) = a.x == b.x
indomain(x, d::Point) = x == d.x
isempty(::Point) = false

approx_indomain(x, d::Point, tolerance) = norm(x-d.x) <= tolerance

isopen(d::Point) = false
isclosed(d::Point) = true

point_in_domain(d::Point) = d.x

for op in (:*,:+,:-)
    @eval begin
        $op(c::Number, d::Point)  = Point($op(c,d.x))
        $op(d::Point,  c::Number) = Point($op(d.x,c))
    end
end


/(d::Point,  c::Number) = Point(d.x/c)
\(c::Number, d::Point)  = Point(c\d.x)

for op in (:+,:-)
    @eval $op(a::Point, b::Point) = Point($op(a.x,b.x))
end


for op in (:*,:+)
    @eval begin
        $op(a::Point, v::AbstractVector) = map(y->$op(a,y),v)
        $op(v::AbstractVector, a::Point) = map(y->$op(y,a),v)
    end
end


function setdiff(d::Interval{L,R,T}, p::Point{T}) where {L,R,T}
    a = leftendpoint(d)
    b = rightendpoint(d)

    a == p.x && return Interval{:open,R}(a,b)
    a < p.x < b && return UnionDomain(Interval{L,:open}(a,p.x)) ∪ UnionDomain(Interval{:open,R}(p.x,b))
    b == p.x && return Interval{L,:open}(a,b)

    return d
end

issubset(p::Point, d::Domain) = p.x ∈ d
