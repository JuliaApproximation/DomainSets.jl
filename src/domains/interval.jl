
iscompact(d::TypedEndpointsInterval{:closed,:closed}) = true
iscompact(d::TypedEndpointsInterval) = false


##################
### An interval
##################

approx_indomain(x, d::AbstractInterval, tolerance) =
    (x >= leftendpoint(d)-tolerance) && (x <= rightendpoint(d)+tolerance)

approx_indomain(x, d::TypedEndpointsInterval{:closed,:closed}, tolerance) =
    (x >= leftendpoint(d)-tolerance) && (x <= rightendpoint(d)+tolerance)
approx_indomain(x, d::TypedEndpointsInterval{:open,:closed}, tolerance) =
    (x > leftendpoint(d)-tolerance) && (x <= rightendpoint(d)+tolerance)
approx_indomain(x, d::TypedEndpointsInterval{:closed,:open}, tolerance) =
    (x >= leftendpoint(d)-tolerance) && (x < rightendpoint(d)+tolerance)
approx_indomain(x, d::TypedEndpointsInterval{:open,:open}, tolerance) =
    (x > leftendpoint(d)-tolerance) && (x < rightendpoint(d)+tolerance)


function point_in_domain(d::AbstractInterval)
    isempty(d) && throw(BoundsError())
    mean(d)
end

isapprox(d1::AbstractInterval, d2::AbstractInterval; kwds...) =
    isapprox(leftendpoint(d1), leftendpoint(d2); kwds...) &&
    isapprox(rightendpoint(d1), rightendpoint(d2); kwds...)


boundary(d::AbstractInterval) = Point(leftendpoint(d)) ∪ Point(rightendpoint(d))

similar_interval(d::AbstractInterval, a, b) = similar_interval(d, promote(a, b)...)

## Some special intervals follow, e.g.:
# - the unit interval [0,1]
# - the 'Chebyshev' interval [-1,1]
# - ...
# Unlike a generic interval, these specific intervals have no data, and only one
# type parameter (T).

"""
The abstract type `FixedInterval` is the supertype of intervals with endpoints
determined by the type, rather than field values. Examples include `UnitInterval`
and `ChebyshevInterval`.
"""
abstract type FixedInterval{L,R,T} <: TypedEndpointsInterval{L,R,T} end
const ClosedFixedInterval{T} = FixedInterval{:closed,:closed,T}

closure(d::TypedEndpointsInterval) = ClosedInterval(endpoints(d)...)
closure(d::ClosedFixedInterval) = d

"""
Return an interval that is similar to the given interval, but with endpoints
`a` and `b` instead.
"""
similar_interval(d::ClosedFixedInterval{T}, a::S, b::S) where {S,T} =
    ClosedInterval{promote_type(float(T),S)}(a, b)

"The closed unit interval [0,1]."
struct UnitInterval{T} <: ClosedFixedInterval{T} end

UnitInterval() = UnitInterval{Float64}()

endpoints(d::UnitInterval{T}) where {T} = (zero(T), one(T))

similardomain(::UnitInterval, ::Type{T}) where {T} = UnitInterval{T}()

canonicaldomain(d::AbstractInterval) = isinf(width(d)) ? d : UnitInterval{eltype(d)}()

tocanonical(d::UnitInterval{T}) where {T} = IdentityMap{T}()
fromcanonical(d::UnitInterval{T}) where {T} = IdentityMap{T}()

tocanonical(d::AbstractInterval{T}) where {T} =
    isinf(width(d)) ? IdentityMap{T}() : interval_map(endpoints(d)..., 0, 1)
fromcanonical(d::AbstractInterval{T}) where {T} =
    isinf(width(d)) ? IdentityMap{T}() : interval_map(0, 1, endpoints(d)...)


"The closed interval [-1,1]."
struct ChebyshevInterval{T} <: ClosedFixedInterval{T}
end

ChebyshevInterval() = ChebyshevInterval{Float64}()

endpoints(d::ChebyshevInterval{T}) where {T} = (-one(T),one(T))

similardomain(::ChebyshevInterval, ::Type{T}) where {T} = ChebyshevInterval{T}()

-(d::ChebyshevInterval) = d


"The half-open positive halfline `[0,∞)`."
struct HalfLine{T} <: FixedInterval{:closed,:open,T} end
HalfLine() = HalfLine{Float64}()

endpoints(d::HalfLine{T}) where {T} = (zero(T), T(Inf))

boundary(d::HalfLine) = Point(leftendpoint(d))

similardomain(::HalfLine, ::Type{T}) where {T} = HalfLine{T}()

# intercept and simplify the definition of IntervalSets
in(x, d::HalfLine) = x >= 0

approx_indomain(x, d::HalfLine, tolerance) = x >= -tolerance

function similar_interval(d::HalfLine{T}, a::S, b::S) where {T,S}
    @assert a == 0
    @assert isinf(b) && b > 0
    HalfLine{promote_type(float(T),S)}()
end

point_in_domain(d::HalfLine) = zero(eltype(d))


"The open negative halfline `(-∞,0)`."
struct NegativeHalfLine{T} <: FixedInterval{:open,:open,T} end
NegativeHalfLine() = NegativeHalfLine{Float64}()

similardomain(::NegativeHalfLine, ::Type{T}) where {T} = NegativeHalfLine{T}()

endpoints(d::NegativeHalfLine{T}) where {T} = (-T(Inf), zero(T))

boundary(d::NegativeHalfLine) = Point(rightendpoint(d))

# intercept and simplify the definition of IntervalSets
in(x, d::NegativeHalfLine) = x < 0

approx_indomain(x, d::NegativeHalfLine, tolerance) = x < tolerance

function similar_interval(d::NegativeHalfLine{T}, a::S, b::S) where {S,T}
    @assert isinf(a) && a < 0
    @assert b == 0
    NegativeHalfLine{promote_type(S,float(T))}()
end

point_in_domain(d::NegativeHalfLine) = -one(eltype(d))

similar_interval(d::Interval{L,R,T}, a::S, b::S) where {L,R,T,S} =
    Interval{L,R,promote_type(float(T),S)}(a, b)


#########################################
# A few set operations with known result
#########################################

# We define an exhaustive list of combinations of the four fixed intervals
# combined above in the routines 'intersect', 'union' and 'setdiff' where
# the output is known explicitly.

# Since fixed intervals are fully determined by their type,
# the result of intersect, union or setdiff is always known for two
# domains of the same type.
intersect(d1::D, d2::D) where {D <: FixedInterval} = d1
union(d1::D, d2::D) where {D <: FixedInterval} = d1
setdiff(d1::D, d2::D) where {D <: FixedInterval} = EmptySpace{eltype(D)}()

# [0,1] ∩ [-1,1] = [0,1]
intersect(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = UnitInterval{T}()
intersect(d1::ChebyshevInterval{T}, d2::UnitInterval{T}) where {T} = UnitInterval{T}()
# [0,1] ∩ [0,∞) = [0,1]
intersect(d1::UnitInterval{T}, d2::HalfLine{T}) where {T} = UnitInterval{T}()
intersect(d1::HalfLine{T}, d2::UnitInterval{T}) where {T} = UnitInterval{T}()
# [0,1] ∩ (-∞,0) = {}
intersect(d1::UnitInterval{T}, d2::NegativeHalfLine{T}) where {T} = EmptySpace{T}()
intersect(d1::NegativeHalfLine{T}, d2::UnitInterval{T}) where {T} = EmptySpace{T}()
# [-1,1] ∩ [0,∞) = [0,1]
intersect(d1::ChebyshevInterval{T}, d2::HalfLine{T}) where {T} = UnitInterval{T}()
intersect(d1::HalfLine{T}, d2::ChebyshevInterval{T}) where {T} = UnitInterval{T}()
# [0,∞) ∩ (-∞,0) = {}
intersect(d1::HalfLine{T}, d2::NegativeHalfLine{T}) where {T} = EmptySpace{T}()
intersect(d1::NegativeHalfLine{T}, d2::HalfLine{T}) where {T} = EmptySpace{T}()

# [0,1] ∪ [-1,1] = [-1,1]
union(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = ChebyshevInterval{T}()
union(d1::ChebyshevInterval{T}, d2::UnitInterval{T}) where {T} = ChebyshevInterval{T}()
# [0,1] ∪ [0,∞) = [0,∞)
union(d1::UnitInterval{T}, d2::HalfLine{T}) where {T} = HalfLine{T}()
union(d1::HalfLine{T}, d2::UnitInterval{T}) where {T} = HalfLine{T}()

# (-∞,0) ∪ [0,∞) = (-∞,∞)
# Note: T<:real to ensure that FullSpace{T} is not larger than intended.
union(d1::NegativeHalfLine{T}, d2::HalfLine{T}) where {T<:Real} = FullSpace{T}()
union(d1::HalfLine{T}, d2::NegativeHalfLine{T}) where {T<:Real} = FullSpace{T}()


# [0,1] ∖ [-1,1] = {}
setdiff(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = EmptySpace{T}()
# [0,1] ∖ [0,∞) = {}
setdiff(d1::UnitInterval{T}, d2::HalfLine{T}) where {T} = EmptySpace{T}()
# [0,1] ∖ (-∞,0) = [0,1]
setdiff(d1::UnitInterval{T}, d2::NegativeHalfLine{T}) where {T} = UnitInterval{T}()
# [-1,1] ∖ (-∞,0) = [0,1]
setdiff(d1::ChebyshevInterval{T}, d2::NegativeHalfLine{T}) where {T} = UnitInterval{T}()
# [0,∞) ∖ (-∞,0) = [0,∞)
setdiff(d1::HalfLine{T}, d2::NegativeHalfLine{T}) where {T} = HalfLine{T}()
# (-∞,0) ∖ [0,1] = (-∞,0)
setdiff(d1::NegativeHalfLine{T}, d2::UnitInterval{T}) where {T} = NegativeHalfLine{T}()
# (-∞,0) ∖ [0,∞) = (-∞,0)
setdiff(d1::NegativeHalfLine{T}, d2::HalfLine{T}) where {T} = NegativeHalfLine{T}()


#################################
# Conversions between intervals
#################################
# Avoid depcrecated warning: Warning: Constructors no longer fall back to `convert`.
# example: A constructor `AbstractInterval{Float64}(::IntervalSets.ClosedInterval{Float64})` should be defined instead.

function convert(::Type{UnitInterval{T}}, d::AbstractInterval) where {T}
    endpoints(d) == (0,1) || throw(InexactError(:convert,UnitInterval,d))
    UnitInterval{T}()
end

function convert(::Type{ChebyshevInterval{T}}, d::AbstractInterval) where {T}
    endpoints(d) == (-1,1)|| throw(InexactError(:convert,ChebyshevInterval,d))
    ChebyshevInterval{T}()
end

UnitInterval{T}(d::AbstractInterval) where T = convert(UnitInterval{T}, d)
ChebyshevInterval{T}(d::AbstractInterval) where T = convert(ChebyshevInterval{T}, d)


function show(io::IO, d::ChebyshevInterval)
    print(io, Interval(d))
    print(io, " (Chebyshev)")
end
function show(io::IO, d::UnitInterval)
    print(io, Interval(d))
    print(io, " (Unit)")
end
function show(io::IO, d::HalfLine)
    print(io, Interval(d))
    print(io, " (HalfLine)")
end


########################
# Arithmetic operations
########################

function setdiff(d1::AbstractInterval{T}, d2::AbstractInterval{T}) where T
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    isempty(d1) && return d1
    isempty(d2) && return d1
    b1 < a2 && return d1
    a1 < a2 ≤ b1 ≤ b2 && return (a1 .. a2)
    a1 < a2 ≤ b2 < b1 && return UnionDomain(a1 .. a2) ∪ UnionDomain(b2 .. b1)
    a2 ≤ a1 < b2 < b1 && return (b2 .. b1)
    a2 ≤ a1 ≤ b1 ≤ b2 && return EmptySpace{T}()

    @assert b2 ≤ a1
    d1
end

setdiff(d1::AbstractInterval, d2::AbstractInterval) = setdiff(promote(d1,d2)...)

switch_open_closed(d::AbstractInterval) = d
switch_open_closed(d::Interval{L,R,T}) where {L,R,T} =
    Interval{R,L,T}(leftendpoint(d),rightendpoint(d))

function map_domain(map::AbstractAffineMap, domain::AbstractInterval)
    le = map(leftendpoint(domain))
    re = map(rightendpoint(domain))
    if le<re
        similar_interval(domain,le,re)
    else
        similar_interval(switch_open_closed(domain),re,le)
    end
end

mapped_domain(invmap::AbstractAffineMap, domain::AbstractInterval) =
    map_domain(inv(invmap), domain)

# Preserve maps when acting on the half line and negative halfline, because
# the map may not always preserve the domain type
map_domain(map::AbstractAffineMap, domain::HalfLine) = MappedDomain(domain, inv(map))
map_domain(map::AbstractAffineMap, domain::NegativeHalfLine) = MappedDomain(domain, inv(map))
