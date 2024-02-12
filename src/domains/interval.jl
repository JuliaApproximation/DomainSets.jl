
iscompact(d::TypedEndpointsInterval{:closed,:closed}) = true
iscompact(d::TypedEndpointsInterval) = false

isinterval(d) = false
isinterval(d::AbstractInterval) = true

hash(d::AbstractInterval, h::UInt) =
    hashrec(isleftopen(d), isrightopen(d), leftendpoint(d), rightendpoint(d), h)

Display.object_parentheses(::AbstractInterval) = true

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


function choice(d::AbstractInterval)
    isempty(d) && throw(BoundsError())
    IntervalSets.mean(d)
end

center(d::AbstractInterval) = IntervalSets.mean(d)

# For an interval of integers, try to find an integer point
function choice(d::AbstractInterval{T}) where {T<:Integer}
    isempty(d) && throw(BoundsError())
    x = round(T, IntervalSets.mean(d))
    if x ∈ d
        x
    else
        # the mean is not inside the interval when rounded, this means
        # we have to choose an endpoint
        if isleftclosed(d)
            leftendpoint(d)
        elseif isrightclosed(d)
            rightendpoint(d)
        else
            # the interval is open and contains no integers in its interior
            throw(BoundsError())
        end
    end
end


boundary(d::AbstractInterval) = Point(leftendpoint(d)) ∪ Point(rightendpoint(d))
corners(d::AbstractInterval) = [leftendpoint(d), rightendpoint(d)]

normal(d::AbstractInterval, x) = (abs(minimum(d)-x) < abs(maximum(d)-x)) ? -one(domaineltype(d)) : one(domaineltype(d))

distance_to(d::AbstractInterval, x) = x ∈ d ? zero(domaineltype(d)) : min(abs(x-supremum(d)), abs(x-infimum(d)))

boundingbox(d::AbstractInterval) = d

volume(d::AbstractInterval) = width(d)

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

convert(::Type{AbstractInterval{T}}, d::FixedInterval{L,R,T}) where {L,R,T} = d
convert(::Type{AbstractInterval{T}}, d::FixedInterval{L,R,S}) where {L,R,S,T} =
    similardomain(d, T)

closure(d::AbstractInterval) = ClosedInterval(endpoints(d)...)
closure(d::ClosedInterval) = d

interior(d::AbstractInterval) = OpenInterval(endpoints(d)...)

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


"The closed interval [-1,1]."
struct ChebyshevInterval{T} <: ClosedFixedInterval{T}
end

ChebyshevInterval() = ChebyshevInterval{Float64}()

endpoints(d::ChebyshevInterval{T}) where {T} = (-one(T),one(T))

similardomain(::ChebyshevInterval, ::Type{T}) where {T} = ChebyshevInterval{T}()

Base.:-(d::ChebyshevInterval) = d


interval_map(a, b, c, d) = interval_map(promote(a,b,c,d)...)

"""
Map the interval `[a,b]` to the interval `[c,d]`.

This function deals with infinite intervals, and the type of the
map returned may depend on the value (finiteness) of the given endpoints.
"""
function interval_map(a::T, b::T, c::T, d::T) where {T}
    FT = float(T)
    if isfinite(a) && isfinite(b) && isfinite(c) && isfinite(d)
        bounded_interval_map(a, b, c, d)
    elseif isfinite(a) && !isfinite(b) && isfinite(c) && !isfinite(d)
        # (a,Inf) to (c,Inf)
        AffineMap(one(FT), c-a)
    elseif isfinite(a) && !isfinite(b) && !isfinite(c) && isfinite(d)
        # (a,Inf) to (Inf,d)
        AffineMap(-one(FT), d+a)
    elseif !isfinite(a) && isfinite(b) && isfinite(c) && !isfinite(d)
        # (Inf,b) to (c,Inf)
        AffineMap(-one(FT), c+b)
    elseif !isfinite(a) && isfinite(b) && !isfinite(c) && isfinite(d)
        # (Inf,b) to (Inf,d)
        AffineMap(one(FT), d-b)
    elseif !isfinite(a) && !isfinite(b) && !isfinite(c) && !isfinite(d)
        if (a < 0) && (b > 0) && (c < 0) && (d > 0)
            # (-Inf,Inf) to (-Inf,Inf)
            StaticIdentityMap{FT}()
        elseif (a < 0) && (b > 0) && (c > 0) && (d < 0)
            # (-Inf,Inf) to (Inf,-Inf)
            LinearMap(-one(FT))
        elseif (a > 0) && (b < 0) && (c < 0) && (d > 0)
            # (Inf,-Inf) to (-Inf,Inf)
            LinearMap(-one(FT))
        elseif (a > 0) && (b < 0) && (c > 0) && (d < 0)
            # (Inf,-Inf) to (Inf,-Inf)
            StaticIdentityMap{FT}()
        elseif (a > 0) && (b > 0) && (c > 0) && (d > 0)
            # (Inf,Inf) to (Inf,Inf)
            StaticIdentityMap{FT}()
        elseif (a < 0) && (b < 0) && (c < 0) && (d < 0)
            # (-Inf,-Inf) to (-Inf,-Inf)
            StaticIdentityMap{FT}()
        else
            throw(ArgumentError("Requested affine map is unbounded"))
        end
    else
        throw(ArgumentError("Requested affine map is unbounded"))
    end
end

"Like interval_map, but guaranteed to return a scalar affine map."
bounded_interval_map(a, b, c, d) = bounded_interval_map(promote(a,b,c,d)...)
bounded_interval_map(a::T, b::T, c::T, d::T) where {T} =
    AffineMap((d-c)/(b-a), c - a*(d-c)/(b-a))

mapto(d1::D, d2::D) where {D <: FixedInterval} = identitymap(d1)
mapto(d1::AbstractInterval, d2::AbstractInterval) =
    interval_map(leftendpoint(d1), rightendpoint(d1), leftendpoint(d2), rightendpoint(d2))


canonicaldomain(d::ClosedInterval{T}) where {T} = ChebyshevInterval{float(T)}()
canonicaldomain(d::FixedInterval) = d

function canonicaldomain(d::Interval{:open,:open,T}) where {T}
    FT = float(T)
    if isfinite(leftendpoint(d)) && isfinite(rightendpoint(d))
        Interval{:open,:open,FT}(-1, 1)
    elseif isfinite(leftendpoint(d)) || isfinite(rightendpoint(d))
        PositiveRealLine{FT}()
    elseif leftendpoint(d) < 0 && rightendpoint(d) > 0
        RealLine{FT}()
    elseif leftendpoint(d) > 0 && rightendpoint(d) < 0
        RealLine{FT}()
    else
        Point(zero(T))
    end
end
function canonicaldomain(d::Interval{:open,:closed,T}) where {T}
    FT = float(T)
    if isfinite(leftendpoint(d)) && isfinite(rightendpoint(d))
        Interval{:open,:closed,FT}(-1, 1)
    elseif isfinite(rightendpoint(d))
        NonnegativeRealLine{FT}()
    else
        throw(ArgumentError("Canonical domains can not be closed at infinity."))
    end
end
function canonicaldomain(d::Interval{:closed,:open,T}) where {T}
    FT = float(T)
    if isfinite(leftendpoint(d)) && isfinite(rightendpoint(d))
        Interval{:closed,:open,FT}(-1, 1)
    elseif isfinite(leftendpoint(d))
        NonnegativeRealLine{FT}()
    else
        throw(ArgumentError("Canonical domains can not be closed at infinity."))
    end
end

mapfrom_canonical(d::ClosedInterval) = bounded_interval_map(-1, 1, endpoints(d)...)
mapfrom_canonical(d::Interval) = mapto(canonicaldomain(d), d)


"The positive halfline `[0,∞)` or `(0,∞)`, left-closed or left-open."
struct HalfLine{T,C} <: FixedInterval{C,:open,T} end
HalfLine() = HalfLine{Float64}()
HalfLine{T}() where {T} = HalfLine{Float64,:closed}()

const NonnegativeRealLine{T} = HalfLine{T,:closed}
const PositiveRealLine{T} = HalfLine{T,:open}
NonnegativeRealLine() = NonnegativeRealLine{Float64}()
PositiveRealLine() = PositiveRealLine{Float64}()

endpoints(d::HalfLine{T}) where {T} = (zero(T), T(Inf))
boundary(d::HalfLine) = Point(leftendpoint(d))
interior(d::HalfLine{T}) where {T} = PositiveRealLine{T}()
closure(d::HalfLine{T}) where {T} = NonnegativeRealLine{T}()

similardomain(::HalfLine{S,C}, ::Type{T}) where {S,T,C} = HalfLine{T,C}()

# intercept and simplify the definition of IntervalSets
in(x, d::NonnegativeRealLine) = x >= 0
in(x, d::PositiveRealLine) = x > 0

approx_indomain(x, d::HalfLine, tolerance) = x >= -tolerance

function similar_interval(d::HalfLine{T,C}, a::S, b::S) where {T,S,C}
    @assert a == 0
    @assert isinf(b) && b > 0
    HalfLine{promote_type(float(T),S),C}()
end

choice(d::NonnegativeRealLine) = zero(domaineltype(d))
choice(d::PositiveRealLine) = one(domaineltype(d))


"The negative halfline `(-∞,0]` or `(-∞,0)`, right-closed or right-open."
struct NegativeHalfLine{T,C} <: FixedInterval{:open,C,T} end
NegativeHalfLine() = NegativeHalfLine{Float64}()
NegativeHalfLine{T}() where {T} = NegativeHalfLine{T,:open}()

const NonpositiveRealLine{T} = NegativeHalfLine{T,:closed}
const NegativeRealLine{T} = NegativeHalfLine{T,:open}
NonpositiveRealLine() = NonpositiveRealLine{Float64}()
NegativeRealLine() = NegativeRealLine{Float64}()

similardomain(::NegativeHalfLine{S,C}, ::Type{T}) where {S,T,C} =
    NegativeHalfLine{T,C}()

endpoints(d::NegativeHalfLine{T}) where {T} = (-T(Inf), zero(T))
boundary(d::NegativeHalfLine) = Point(rightendpoint(d))
interior(d::NegativeHalfLine{T}) where {T} = NegativeRealLine{T}()
closure(d::NegativeHalfLine{T}) where {T} = NonpositiveRealLine{T}()

# intercept and simplify the definition of IntervalSets
in(x, d::NonpositiveRealLine) = x <= 0
in(x, d::NegativeRealLine) = x < 0

approx_indomain(x, d::NegativeHalfLine, tolerance) = x < tolerance

function similar_interval(d::NegativeHalfLine{T,C}, a::S, b::S) where {S,T,C}
    @assert isinf(a) && a < 0
    @assert b == 0
    NegativeHalfLine{promote_type(S,float(T)),C}()
end

choice(d::NegativeRealLine) = -one(domaineltype(d))
choice(d::NonpositiveRealLine) = zero(domaineltype(d))


"The real line `(-∞,∞)`."
struct RealLine{T} <: FixedInterval{:open,:open,T} end
RealLine() = RealLine{Float64}()

choice(d::RealLine) = zero(eltype(d))

similardomain(::RealLine, ::Type{T}) where {T} = RealLine{T}()

function similar_interval(d::RealLine{T}, a::S, b::S) where {S,T}
    @assert isinf(a) && a < 0
    @assert isinf(b) && b > 0
    RealLine{promote_type(S,float(T))}()
end

endpoints(d::RealLine{T}) where {T} = (-T(Inf), T(Inf))
boundary(d::RealLine) = emptyspace(d)
interior(d::RealLine) = d

isfullspace(d::RealLine) = true

similar_interval(d::Interval{L,R,T}, a::S, b::S) where {L,R,T,S} =
    Interval{L,R,promote_type(float(T),S)}(a, b)


#########################################
# A few set operations with known result
#########################################

# We define an exhaustive list of combinations of the four fixed intervals
# combined above in the routines 'intersect', 'union' and 'setdiff' where
# the output is known explicitly.

# Override the definition of Intervals.jl for FixedInterval's defined here
Base.union(d1::FixedInterval, d2::FixedInterval) = uniondomain(d1, d2)
Base.intersect(d1::FixedInterval, d2::FixedInterval) = intersectdomain(d1, d2)

# Promotion to joint type T
uniondomain(d1::TypedEndpointsInterval, d2::TypedEndpointsInterval) =
    uniondomain(promote_domains(d1,d2)...)
intersectdomain(d1::TypedEndpointsInterval, d2::TypedEndpointsInterval) =
    intersectdomain(promote_domains(d1,d2)...)
setdiffdomain(d1::AbstractInterval, d2::AbstractInterval) =
    setdiffdomain(promote_domains(d1,d2)...)



# type-unstable union of intervals. This function differs from `union` in
# IntervalSets.jl because that one throws an error if the intervals do not overlap
function uniondomain(d1::TypedEndpointsInterval{L1,R1,T}, d2::TypedEndpointsInterval{L2,R2,T}) where {L1,R1,L2,R2,T}
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    # are they empty?
    isempty(d1) && return d2
    isempty(d2) && return d1
    # are they equal?
    d1 == d2 && return d1
    # does one lie within the other?
    a2 > a1 && b2 < b1 && return d1
    a1 > a2 && b1 < b2 && return d2
    # are they disjoint?
    b2 < a1 && return UnionDomain(d1, d2)   # return a UnionDomain for disjoint intervals
    b1 < a2 && return UnionDomain(d1, d2)
    if a1 ≤ b1 == a2 ≤ b2
        if R1 == :closed || L2 == :closed
            return Interval{L1,R2,T}(a1,b2)
        else
            return UnionDomain(Interval{L1,:open,T}(a1,b1), Interval{:open,R2,T}(a2,b2))
        end
    end
    # at this stage they must overlap
    if a1 < a2
        a = a1
        L = L1
    elseif a1 == a2
        a = a1
        L = (L1==L2==:open) ? :open : :closed
    else
        a = a2
        L = L2
    end
    if b1 > b2
        b = b1
        R = R1
    elseif b1 == b2
        b = b1
        R = (R1==R2==:open) ? :open : :closed
    else
        b = b2
        R = R2
    end
    Interval{L,R,T}(a, b)
end

function intersectdomain(d1::TypedEndpointsInterval{L1,R1,T}, d2::TypedEndpointsInterval{L2,R2,T}) where {L1,R1,L2,R2,T}
    # go back to the definition of IntervalSets.jl
    # to that end, convert to Interval in order to avoid stack overflow
    d = intersect(Interval(d1), Interval(d2))
    if isempty(d)
        # but avoid returning an interval like 2..1
        emptyspace(d)
    elseif leftendpoint(d)==rightendpoint(d)
        # and avoid an interval like 1..1
        Point(leftendpoint(d))
    else
        d
    end
end

isequaldomain(d1::TypedEndpointsInterval, d2::Point) =
    isclosedset(d1) && (leftendpoint(d1)==rightendpoint(d1)==pointval(d2))

# Since fixed intervals are fully determined by their type,
# the result of intersect, union or setdiff is always known for two
# domains of the same type.
intersectdomain(d1::D, d2::D) where {D <: FixedInterval} = d1
uniondomain(d1::D, d2::D) where {D <: FixedInterval} = d1
setdiffdomain(d1::D, d2::D) where {D <: FixedInterval} = emptyspace(d1)

# [0,1] ∩ [-1,1] = [0,1]
intersectdomain(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = UnitInterval{T}()
intersectdomain(d1::ChebyshevInterval{T}, d2::UnitInterval{T}) where {T} = UnitInterval{T}()
# [0,1] ∩ [0,∞) = [0,1]
intersectdomain(d1::UnitInterval{T}, d2::NonnegativeRealLine{T}) where {T} = UnitInterval{T}()
intersectdomain(d1::NonnegativeRealLine{T}, d2::UnitInterval{T}) where {T} = UnitInterval{T}()
# [0,1] ∩ (-∞,0) = {}
intersectdomain(d1::UnitInterval{T}, d2::NegativeRealLine{T}) where {T} = emptyspace(T)
intersectdomain(d1::NegativeRealLine{T}, d2::UnitInterval{T}) where {T} = emptyspace(T)
# [-1,1] ∩ [0,∞) = [0,1]
intersectdomain(d1::ChebyshevInterval{T}, d2::NonnegativeRealLine{T}) where {T} = UnitInterval{T}()
intersectdomain(d1::NonnegativeRealLine{T}, d2::ChebyshevInterval{T}) where {T} = UnitInterval{T}()
# open and closed halfline
intersectdomain(d1::HalfLine{T}, d2::NegativeRealLine{T}) where {T} = emptyspace(T)
intersectdomain(d1::NegativeRealLine{T}, d2::HalfLine{T}) where {T} = emptyspace(T)
intersectdomain(d1::NonnegativeRealLine{T}, d2::NonpositiveRealLine{T}) where {T} = Point(zero(T))
intersectdomain(d1::NonpositiveRealLine{T}, d2::NonnegativeRealLine{T}) where {T} = Point(zero(T))
intersectdomain(d1::NonnegativeRealLine{T}, d2::PositiveRealLine{T}) where {T} = d2
intersectdomain(d1::PositiveRealLine{T}, d2::NonnegativeRealLine{T}) where {T} = d1
# [a,b] ∩ (-∞,∞) = [a,b]
intersectdomain(d1::AbstractInterval{T}, d2::RealLine{T}) where {T} = d1
intersectdomain(d1::TypedEndpointsInterval{L,R,T}, d2::RealLine{T}) where {L,R,T} = d1
intersectdomain(d1::RealLine{T}, d2::AbstractInterval{T}) where {T} = d2
intersectdomain(d1::RealLine{T}, d2::TypedEndpointsInterval{L,R,T}) where {L,R,T} = d2


# [0,1] ∪ [-1,1] = [-1,1]
uniondomain(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = ChebyshevInterval{T}()
uniondomain(d1::ChebyshevInterval{T}, d2::UnitInterval{T}) where {T} = ChebyshevInterval{T}()
# [0,1] ∪ [0,∞) = [0,∞)
uniondomain(d1::UnitInterval{T}, d2::NonnegativeRealLine{T}) where {T} = d2
uniondomain(d1::NonnegativeRealLine{T}, d2::UnitInterval{T}) where {T} = d1
# open and closed halflines
uniondomain(d1::NonnegativeRealLine{T}, d2::NonpositiveRealLine{T}) where {T} = RealLine{T}()
uniondomain(d1::NonpositiveRealLine{T}, d2::NonnegativeRealLine{T}) where {T} = RealLine{T}()
uniondomain(d1::PositiveRealLine{T}, d2::NonpositiveRealLine{T}) where {T} = RealLine{T}()
uniondomain(d1::NonpositiveRealLine{T}, d2::PositiveRealLine{T}) where {T} = RealLine{T}()
uniondomain(d1::NonnegativeRealLine{T}, d2::NegativeRealLine{T}) where {T} = RealLine{T}()
uniondomain(d1::NegativeRealLine{T}, d2::NonnegativeRealLine{T}) where {T} = RealLine{T}()
uniondomain(d1::NonnegativeRealLine{T}, d2::PositiveRealLine{T}) where {T} = d1
uniondomain(d1::PositiveRealLine{T}, d2::NonnegativeRealLine{T}) where {T} = d2
uniondomain(d1::NonpositiveRealLine{T}, d2::NegativeHalfLine{T}) where {T} = d1
uniondomain(d1::NegativeHalfLine{T}, d2::NonpositiveRealLine{T}) where {T} = d2
# [a,b] ∪ (-∞,∞) = (-∞,∞)
uniondomain(d1::AbstractInterval{T}, d2::RealLine{T}) where {T} = d2
uniondomain(d1::TypedEndpointsInterval{L,R,T}, d2::RealLine{T}) where {L,R,T} = d2
uniondomain(d1::RealLine{T}, d2::AbstractInterval{T}) where {T} = d1
uniondomain(d1::RealLine{T}, d2::TypedEndpointsInterval{L,R,T}) where {L,R,T} = d1


# [0,1] ∖ [-1,1] = {}
setdiffdomain(d1::UnitInterval{T}, d2::ChebyshevInterval{T}) where {T} = emptyspace(T)
# [0,1] ∖ [0,∞) = {}
setdiffdomain(d1::UnitInterval{T}, d2::NonnegativeRealLine{T}) where {T} = emptyspace(T)
# [0,1] ∖ (-∞,0) = [0,1]
setdiffdomain(d1::UnitInterval{T}, d2::NegativeRealLine{T}) where {T} = d1
# [-1,1] ∖ (-∞,0) = [0,1]
setdiffdomain(d1::ChebyshevInterval{T}, d2::NegativeRealLine{T}) where {T} = UnitInterval{T}()
# [0,∞) ∖ (-∞,0) = [0,∞)
setdiffdomain(d1::HalfLine{T}, d2::NegativeRealLine{T}) where {T} = d1
setdiffdomain(d1::HalfLine{T}, d2::NonpositiveRealLine{T}) where {T} = HalfLine{T,:open}()
# (-∞,0) ∖ [0,1] = (-∞,0)
setdiffdomain(d1::NegativeRealLine{T}, d2::UnitInterval{T}) where {T} = d1
# (-∞,0) ∖ [0,∞) = (-∞,0)
setdiffdomain(d1::NegativeHalfLine{T}, d2::PositiveRealLine{T}) where {T} = d1
setdiffdomain(d1::NegativeHalfLine{T}, d2::NonnegativeRealLine{T}) where {T} = NegativeRealLine{T}()
# (-∞,∞) ∖ [0,∞) = (-∞,0)
setdiffdomain(d1::RealLine{T}, d2::NonnegativeRealLine{T}) where {T} = NegativeRealLine{T}()
setdiffdomain(d1::RealLine{T}, d2::PositiveRealLine{T}) where {T} = NonpositiveRealLine{T}()
setdiffdomain(d1::RealLine{T}, d2::NonpositiveRealLine{T}) where {T} = PositiveRealLine{T}()
setdiffdomain(d1::RealLine{T}, d2::NegativeRealLine{T}) where {T} = NonnegativeRealLine{T}()


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
function show(io::IO, d::RealLine)
    print(io, Interval(d))
    print(io, " (RealLine)")
end


########################
# Arithmetic operations
########################

complement(L) = L == :open ? :closed : :open

function setdiffdomain(d1::TypedEndpointsInterval{L1,R1,T}, d2::TypedEndpointsInterval{L2,R2,T}) where {L1,R1,L2,R2,T}
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    isempty(d1) && return d1
    isempty(d2) && return d1
    # intervals aren't empty: we now know that a1 ≤ b1 and a2 ≤ b2
    d1 == d2 && return emptyspace(T)
    # Order: a1 b1 a2 b2
    b1 < a2 && return d1
    if a1 < a2 == b1 ≤ b2
        # if a2==b1==b2 then [a2,b2] is closed, because it is non-empty
        (L2 == :open && R1 == :closed) ? R = :closed : R = :open
        return Interval{L1,R,T}(a1, a2)
    end
    # Order: a1 a2 b1 b2
    if a1 < a2 < b1 ≤ b2
        (L2 == :open) ? R = :closed : R = :open
        return Interval{L1,R,T}(a1, a2)
    end
    # Order: a1 a2 b2 b1
    if a1 < a2 < b2 < b1
        return uniondomain(Interval{L1,complement(L1),T}(a1,a2), Interval{complement(R2),R1,T}(b2,b1))
    end
    if a1 < a2 == b2 < b1
        # since a2==b2 and the interval isn't empty, [a2,b2] is closed
        return uniondomain(Interval{L1,:open,T}(a1,a2), Interval{:open,R1,T}(b2,b1))
    end
    # Order: a2 a1 b2 b1
    if a2 ≤ a1 < b2 < b1
        return Interval{complement(R2),R1,T}(b2,b1)
    end
    # Order: a2 a1 b1 b2
    if a2 < a1 ≤ b1 < b2
        return emptyspace(T)
    end
    if a2 ≤ a1 ≤ b1 ≤ b2
        return (L2 == :open && R2 == :open) ? d1 : emptyspace(T)
    end
    # Order: a2 b2 a1 b1
    if b2 == a1 ≤ b1
        return (R2 == :open) ? d1 : Interval{:open,R1}(a1,b1)
    end
    if b2 < a1
        return d1
    end
    error("Can't reach this line: please file an issue.")
end


switch_open_closed_if_applicable(d::AbstractInterval) = d
switch_open_closed_if_applicable(d::Interval{L,R,T}) where {L,R,T} =
    Interval{R,L,T}(leftendpoint(d),rightendpoint(d))

function map_domain(map::AbstractAffineMap{<:Number}, domain::AbstractInterval)
    le = map(leftendpoint(domain))
    re = map(rightendpoint(domain))
    if le<re
        similar_interval(domain,le,re)
    else
        similar_interval(switch_open_closed_if_applicable(domain),re,le)
    end
end

mapped_domain(invmap::AbstractAffineMap{<:Number}, domain::AbstractInterval) =
    map_domain(inverse(invmap), domain)

# Preserve maps when acting on the half line and negative halfline, because
# the map may not always preserve the domain type
map_domain(map::AbstractAffineMap{<:Number}, domain::HalfLine) = MappedDomain(inverse(map), domain)
map_domain(map::AbstractAffineMap{<:Number}, domain::NegativeHalfLine) = MappedDomain(inverse(map), domain)
