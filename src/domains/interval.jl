# interval.jl

##################
### An interval
##################

abstract type AbstractInterval{T} <: Domain{T}
end

ndims(::Type{D}) where {D <: AbstractInterval} = 1

"The left endpoint of the interval."
leftendpoint(d::AbstractInterval) = d.a

"The right endpoint of the interval."
rightendpoint(d::AbstractInterval) = d.b


## Some special intervals
# - the unit interval [0,1]
# - the 'Chebyshev' interval [-1,1]

"""
The abstract type `FixedInterval` is the supertype of intervals with endpoints
determined by the type, rather than field values. Examples include `UnitInterval`
and `ChebyshevInterval`.
"""
abstract type FixedInterval{T} <: AbstractInterval{T}
end

# We assume by default that fixed intervals are closed. Override if they aren't.
isclosed(d::FixedInterval) = true
isopen(d::FixedInterval) = false

# We also assume that the domain is compact. Override if it is not.
iscompact(d::FixedInterval) = true

# We assume a closed domain for membership.
indomain(x, d::FixedInterval) = leftendpoint(d) <= x <= rightendpoint(d)

"""
Return an interval that is similar to the given interval, but with endpoints
`a` and `b` instead.
"""
# Assume a closed interval by default
similar_interval(d::FixedInterval{T}, a, b) where {T} = ClosedInterval{T}(a, b)


"The closed unit interval [0,1]."
struct UnitInterval{T} <: FixedInterval{T}
end

# By default we create a Float64 interval
UnitInterval() = UnitInterval{Float64}()

leftendpoint(d::UnitInterval{T}) where {T} = zero(T)
rightendpoint(d::UnitInterval{T}) where {T} = one(T)


"The closed interval [-1,1]."
struct ChebyshevInterval{T} <: FixedInterval{T}
end

leftendpoint(d::ChebyshevInterval{T}) where {T} = -one(T)
rightendpoint(d::ChebyshevInterval{T}) where {T} = one(T)


"The real line `(-∞,∞)`."
struct RealLine{T} <: FixedInterval{T}
end

leftendpoint(d::RealLine{T}) where {T} = zero(T)
rightendpoint(d::RealLine{T}) where {T} = T(Inf)

iscompact(d::RealLine) = false

indomain(x::T, d::RealLine{T}) where {T} = true

function similar_interval(d::RealLine, a, b)
    @assert isinf(a) && a < 0
    @assert isinf(b) && b > 0
    d
end


"The half-open positive halfline `[0,∞)`."
struct HalfLine{T} <: FixedInterval{T}
end

leftendpoint(d::HalfLine{T}) where {T} = zero(T)
rightendpoint(d::HalfLine{T}) where {T} = T(Inf)

# A half-open domain is neither open nor closed
isclosed(d::HalfLine) = false
isopen(d::HalfLine) = false

iscompact(d::HalfLine) = false

indomain(x, d::HalfLine) = x >= 0

function similar_interval(d::HalfLine, a, b)
    @assert a == 0
    @assert isinf(b) && b > 0
    d
end


"The open negative halfline `(-∞,0)`."
struct NegativeHalfLine{T} <: FixedInterval{T}
end

leftendpoint(d::NegativeHalfLine{T}) where {T} = -T(Inf)
rightendpoint(d::NegativeHalfLine{T}) where {T} = zero(T)

isclosed(d::NegativeHalfLine) = false
isopen(d::NegativeHalfLine) = true

iscompact(d::NegativeHalfLine) = false

indomain(x, d::NegativeHalfLine) = x < 0

function similar_interval(d::NegativeHalfLine, a, b)
    @assert isinf(a) && a < 0
    @assert b == 0
    d
end


"""
A general interval with endpoints `a` and `b`. The interval can be open or
closed at each endpoint. This is determined by the `L` and `R` type parameters,
which may be `:open` or `:closed`.
"""
struct Interval{L,R,T} <: AbstractInterval{T}
    a     ::  T
    b     ::  T

    # The interval defaults to the unit interval.
    function Interval{L,R,T}(a = 0, b = 1) where {L,R,T}
        # We only allow finite values for a and b
        @assert isfinite(a)
        @assert isfinite(b)
        new{L,R,T}(a,b)
    end
end

"A closed interval `[a,b]`."
const ClosedInterval{T} = Interval{:closed,:closed,T}

"An open interval `(a,b)`."
const OpenInterval{T} = Interval{:open,:open,T}

"A half-open interval `(a,b]`."
const HalfOpenLeftInterval{T} = Interval{:open,:closed,T}

"A half-open interval `[a,b)`."
const HalfOpenRightInterval{T} = Interval{:closed,:open,T}

# By default we create a closed interval
interval(args...) = closed_interval(args...)

closed_interval(args...) = ClosedInterval(args...)
open_interval(args...) = OpenInterval(args...)

# By default we use a Float64 type
Interval{L,R}() where {L,R} = Interval{L,R,Float64}()

Interval{L,R}(a::T, b::S) where {L,R,T,S} = Interval{L,R}(promote(a,b)...)

Interval{L,R}(a::T, b::T) where {L,R,T} = Interval{L,R,T}(a, b)

leftendpoint(d::Interval) = d.a
rightendpoint(d::Interval) = d.b

# The interval is closed if it is closed at both endpoints, and open if it
# is open at both endpoints. In all other cases, it is neither open nor closed.
isclosed(d::ClosedInterval) = true
isopen(d::OpenInterval) = true
isclosed(d::Interval) = false
isopen(d::Interval) = false

iscompact(d::Interval) = true

indomain(x, d::OpenInterval) = d.a < x < d.b
indomain(x, d::ClosedInterval) = d.a <= x <= d.b
indomain(x, d::HalfOpenLeftInterval) = d.a < x <= d.b
indomain(x, d::HalfOpenRightInterval) = d.a <= x < d.b


similar_interval(d::Interval{L,R,T}, a, b) where {L,R,T} =
    Interval{L,R,T}(a, b)


#################################
# Conversions between intervals
#################################

"Return a default interval (the unit interval `[0,1]`)."
interval() = UnitInterval()


convert(::Type{Interval{L,R,T}}, d::AbstractInterval{S}) where {L,R,T,S} =
    Interval{L,R,T}(leftendpoint(d), rightendpoint(d))

function convert(::Type{UnitInterval{T}}, d::AbstractInterval{S}) where {T,S}
    @assert leftendpoint(d) == 0
    @assert rightendpoint(d) == 1
    UnitInterval{T}()
end

function convert(::Type{ChebyshevInterval{T}}, d::AbstractInterval{S}) where {T,S}
    @assert leftendpoint(d) == -1
    @assert rightendpoint(d) == 1
    UnitInterval{T}()
end



########################
# Arithmetic operations
########################

# Some computations with intervals simplify without having to use a mapped domain.
# This is only the case for Interval{L,R,T}, and not for any of the FixedIntervals
# because the endpoints of the latter are, well, fixed.
(+)(d::Interval, x::Number) = similar_interval(d, leftendpoint(d)+x, rightendpoint(d)+x)
(*)(a::Number, d::Interval) = similar_interval(d, a*leftendpoint(d), a*rightendpoint(d))
(/)(d::Interval, a::Number) = similar_interval(d, leftendpoint(d)/a, rightendpoint(d)/a)


show(io::IO, d::AbstractInterval) = print(io, "the interval [", d.a, ", ", d.b, "]")

function union(d1::Interval{L1,R1,T}, d2::Interval{L2,R2,T}) where {L1,R1,L2,R2,T}
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    if (b1 < a2) || (a1 > b2)
        UnionDomain(d1, d2)
    else
        # TODO: add some logic to determine open and closed nature of endpoints of new interval
        Interval(min(a1, a2), max(b1, b2))
    end
end

function intersect(d1::Interval{L1,R1,T}, d2::Interval{L2,R2,T}) where {L1,R1,L2,R2,T}
    a1 = leftendpoint(d1)
    b1 = rightendpoint(d1)
    a2 = leftendpoint(d2)
    b2 = rightendpoint(d2)

    if (b1 < a2) || (a1 > b2)
        EmptyDomain(Val{1}())
    else
        # TODO: add some logic to determine open and closed nature of endpoints of new interval
        Interval(max(a1, a2), min(b1, b2))
    end
end
