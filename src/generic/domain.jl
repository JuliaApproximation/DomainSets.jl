# Definition of the abstract Domain type and its interface

eltype(::Type{<:Domain{T}}) where {T} = T
prectype(::Type{<:Domain{T}}) where {T} = prectype(T)
numtype(::Type{<:Domain{T}}) where {T} = numtype(T)

convert_numtype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{convert_numtype(T,U)}, d)
convert_prectype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{convert_prectype(T,U)}, d)

Domain(d) = convert(Domain, d)

convert(::Type{Domain{T}}, d::Domain{T}) where {T} = d
convert(::Type{Domain{T}}, d::Domain{S}) where {S,T} = similardomain(d, T)

"A `EuclideanDomain` is any domain whose eltype is `SVector{N,T}`."
const EuclideanDomain{N,T} = Domain{SVector{N,T}}

"What is the Euclidean dimension of the domain?"
dimension(::Domain{<:Number}) = 1
dimension(::EuclideanDomain{N}) where {N} = N
dimension(::Domain{<:NTuple{N,Any}}) where {N} = N

"""
A `VectorDomain` is any domain whose eltype is `Vector{T}`. In this case
the dimension of the domain is not included in its type.
"""
const VectorDomain{T} = Domain{Vector{T}}

const AbstractVectorDomain{T} = Domain{<:AbstractVector{T}}

# At the level of Domain we attempt to promote the arguments to compatible
# types, then we invoke indomain. Concrete subtypes should implement
# indomain. They may also implement `in` in order to accept more types.
in(x, d::Domain) = _in(x, d)
_in(x, d::Domain{Any}) = indomain(x, d)
_in(x::S, d::Domain{T}) where {S,T} = _in(x, d, promote_type(S,T))
_in(x::T, domain::Domain{T}, ::Type{T}) where {T} = indomain(x, domain)
_in(x::S, domain::Domain{T}, ::Type{U}) where {S,T,U} = indomain(convert(U, x), convert(Domain{U}, domain))

_in(x::Tuple, d::Domain{<:Tuple}) = indomain(x, d)

function _in(x::S, domain::Domain{T}, ::Type{Any}) where {S,T}
    @warn "in(x,domain): incompatible types $(S) and $(T). Returning false."
    return false
end

# Treatment of abstract vectors for Euclidean domains and vector domains:

# - We allow any abstract vector for Vector domains
_in(x::AbstractVector, domain::AbstractVectorDomain) = indomain(x, domain)

# - for Euclidean domains, we convert x to SVector
function _in(x::AbstractVector{S}, domain::EuclideanDomain{N,T}) where {N,S,T}
    # Avoid an error in an attempt to convert to SVector{N,T}
    if length(x) != N
        @warn "in(x,domain): incompatible dimension $(length(x)) of x and $(N) of the domain. Returning false."
        false
    else
        U = promote_type(S,T)
        indomain(SVector{N,U}(x), convert(Domain{SVector{N,U}}, domain))
    end
end


function indomain(x, domain::Domain)
    @warn "Please implement `indomain` for domain $(domain) for points of type $(typeof(x)). Returning false."
    return false
end




"""
Return a suitable tolerance to use for verifying whether a point is close to
a domain. Typically, the tolerance is close to the precision limit of the numeric
type associated with the domain.
"""
default_tolerance(d::Domain) = default_tolerance(prectype(d))
default_tolerance(::Type{T}) where {T <: AbstractFloat} = 100eps(T)


"""
`approx_in(x, domain[, tolerance])`

Verify whether a point lies in the given domain with a certain tolerance.

The tolerance has to be positive. The meaning of the tolerance, in relation
to the possible distance of the point to the domain, is domain-dependent.
Usually, if the outcome is true, it means that the distance of the point to
the domain is smaller than a constant times the tolerance. That constant may
depend on the domain.

Up to inexact computations due to floating point numbers, it should also be
the case that `approx_in(x, d, 0) == in(x,d)`. This implies that `approx_in`
reflects whether a domain is open or closed.
"""
approx_in(x, domain::Domain) = approx_in(x, domain, default_tolerance(domain))

# We mimick the promotion rules of `in` above.
function approx_in(x, d::Domain, tolerance)
    tolerance >= 0 || error("Tolerance has to be positive in `approx_in`.")
    _approx_in(x, d, tolerance)
end
_approx_in(x::S, d::Domain{T}, tolerance) where {S,T} =
    _approx_in(x, d, tolerance, promote_type(S,T))
_approx_in(x::S, domain::Domain{T}, tolerance, ::Type{U}) where {S,T,U} =
    approx_indomain(convert(U, x), convert(Domain{U}, domain), tolerance)
_approx_in(x::T, domain::Domain{T}, tolerance, ::Type{T}) where {T} =
    approx_indomain(x, domain, tolerance)

_approx_in(x::Tuple, d::Domain{<:Tuple}, tolerance) = approx_indomain(x, d)

function _approx_in(x::S, domain::Domain{T}, tolerance, ::Type{Any}) where {S,T}
    @warn "in: incompatible types $(S) and $(T). Returning false."
    return false
end

# Special cases analogous to the ones above for `in` follow:
_approx_in(x::AbstractVector, domain::AbstractVectorDomain, tolerance) =
    approx_indomain(x, domain, tolerance)

function _approx_in(x::AbstractVector{S}, domain::EuclideanDomain{N,T}, tolerance) where {N,S,T}
    # Avoid an error in an attempt to convert to SVector{N,T}
    if length(x) != N
        @warn "approx_in(x,domain): incompatible dimension $(length(x)) of x and $(N) of the domain. Returning false."
        false
    else
        U = promote_type(S,T)
        approx_indomain(SVector{N,U}(x), convert(Domain{SVector{N,U}}, domain), tolerance)
    end
end


function approx_indomain(x, domain::Domain, tolerance)
    @warn "Please consider implementing `approx_in` for domain $(domain) for points of type $(typeof(x))."
    # Fall back to in without tolerance
    return in(x, domain)
end

isapprox(d1::Domain, d2::Domain; kwds...) = d1 == d2

isreal(d::Domain) = isreal(eltype(d))

infimum(d::Domain) = minimum(d)  # if the minimum exists, then it is also the infimum
supremum(d::Domain) = maximum(d)  # if the maximum exists, then it is also the supremum


# override minimum and maximum for closed sets
function boundary end
const ∂ = boundary
