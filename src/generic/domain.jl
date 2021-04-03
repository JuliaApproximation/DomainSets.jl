# Definition of the abstract Domain type and its interface

# The type Domain{T} is defined in IntervalSets.jl

eltype(::Type{<:Domain{T}}) where {T} = T
prectype(::Type{<:Domain{T}}) where {T} = prectype(T)
numtype(::Type{<:Domain{T}}) where {T} = numtype(T)

convert_numtype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{convert_numtype(T,U)}, d)
convert_prectype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{convert_prectype(T,U)}, d)

Domain(d) = convert(Domain, d)

# Concrete types can implement similardomain(d, ::Type{T}) where {T}
# to support convert(Domain{T}, d) functionality.
convert(::Type{Domain{T}}, d::Domain{T}) where {T} = d
convert(::Type{Domain{T}}, d::Domain{S}) where {S,T} = similardomain(d, T)

"Promote the given domains to have a common element type."
promote_domains() = ()
promote_domains(domains) = convert_eltype.(mapreduce(eltype, promote_type, domains), domains)

convert_eltype(::Type{T}, d::Domain) where {T} = convert(Domain{T}, d)
convert_eltype(::Type{T}, d) where {T} = _convert_eltype(T, d, eltype(d))
_convert_eltype(::Type{T}, d, ::Type{T}) where {T} = d
_convert_eltype(::Type{T}, d, ::Type{S}) where {S,T} =
    error("Don't know how to convert the `eltype` of $(d).")
# Some standard cases
convert_eltype(::Type{T}, d::AbstractArray) where {T} = convert(AbstractArray{T}, d)
convert_eltype(::Type{T}, d::Set) where {T} = convert(Set{T}, d)

compatible_eltype(d1, d2) = isconcretetype(promote_type(eltype(d1),eltype(d2)))



"A `EuclideanDomain` is any domain whose eltype is `SVector{N,T}`."
const EuclideanDomain{N,T} = Domain{SVector{N,T}}

"A `VectorDomain` is any domain whose eltype is `Vector{T}`."
const VectorDomain{T} = Domain{Vector{T}}

const AbstractVectorDomain{T} = Domain{<:AbstractVector{T}}

"What is the Euclidean dimension of the domain?"
dimension(::Domain{<:Number}) = 1
dimension(::EuclideanDomain{N}) where {N} = N
dimension(::Domain{<:NTuple{N,Any}}) where {N} = N
# We can't say anything generically about VectorDomain's here


"Is the given combination of point and domain compatible?"
iscompatible(x, d::Domain) = _iscompatible(x, d, promote_type(typeof(x),eltype(d)))
_iscompatible(x, d, ::Type{T}) where {T} = true
_iscompatible(x, d, ::Type{Any}) = false
_iscompatible(x, d::Domain{Any}, ::Type{Any}) = true

# Some generic cases where we can be sure:
iscompatible(x::SVector{N}, ::EuclideanDomain{N}) where {N} = true
iscompatible(x::SVector{N}, ::EuclideanDomain{M}) where {N,M} = false
iscompatible(x::AbstractVector, ::EuclideanDomain{N}) where {N} = length(x)==N

compatible_or_false(x, domain) =
    iscompatible(x, domain) ? true : (@warn "`in`: incompatible combination of point: $(typeof(x)) and domain eltype: $(eltype(domain)). Returning false."; false)

compatible_or_false(x::AbstractVector, domain::AbstractVectorDomain) =
    iscompatible(x, domain) ? true : (@warn "`in`: incompatible combination of vector with length $(length(x)) and domain '$(domain)' with dimension $(dimension(domain)). Returning false."; false)


"Promote point and domain to compatible types."
promote_pair(x, d) = _promote_pair(x, d, promote_type(typeof(x),eltype(d)))
_promote_pair(x, d, ::Type{T}) where {T} = convert(T, x), convert(Domain{T}, d)
_promote_pair(x, d, ::Type{Any}) = x, d
# Some exceptions:
# - matching types: avoid promotion just in case it is expensive
promote_pair(x::T, d::Domain{T}) where {T} = x, d
# - `Any` domain
promote_pair(x, d::Domain{Any}) = x, d
# - tuples: these are typically composite domains and the elements may be promoted later on
promote_pair(x::Tuple, d::Domain{<:Tuple}) = x, d
# - abstract vectors: promotion may be expensive
promote_pair(x::AbstractVector, d::AbstractVectorDomain) = x, d
# - SVector: promotion is likely cheap
promote_pair(x::AbstractVector{S}, d::EuclideanDomain{N,T}) where {N,S,T} =
    _promote_pair(x, d, SVector{N,promote_type(S,T)})


# At the level of Domain we attempt to promote the arguments to compatible
# types, then we invoke indomain. Concrete subtypes should implement
# indomain. They may also implement `in` in order to accept more types.
#
# Note that if the type of x and the element type of d don't match, then
# both x and the domain may be promoted (using convert(Domain{T}, d)) syntax).
in(x, d::Domain) = compatible_or_false(x, d) && indomain(promote_pair(x, d)...)


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
approx_in(x, d::Domain) = approx_in(x, d, default_tolerance(d))

function compatible_or_false(x, d, tol)
    tol >= 0 || error("Tolerance has to be positive in `approx_in`.")
    compatible_or_false(x, d)
end

approx_in(x, d, tol) = compatible_or_false(x, d, tol) && approx_indomain(promote_pair(x, d)..., tol)

# Fallback to `in`
approx_indomain(x, d, tol) = in(x, d)


isapprox(d1::Domain, d2::Domain; kwds...) = d1 == d2

isreal(d::Domain) = isreal(eltype(d))

infimum(d::Domain) = minimum(d)  # if the minimum exists, then it is also the infimum
supremum(d::Domain) = maximum(d)  # if the maximum exists, then it is also the supremum


# override minimum and maximum for closed sets
function boundary end
const ∂ = boundary

## Mappings to and from a canonical domain

"""
Return an associated canonical domain, if any, of the given domain.

Optionally, additional arguments can be used to specify one of several
canonical domains.
"""
canonicaldomain(d::Domain, args...) = d

"Return a map from the domain to its canonical domain."
tocanonical(d, args...) = IdentityMap{eltype(d)}()

"Return a map to a domain from its canonical domain."
fromcanonical(d, args...) = IdentityMap{eltype(d)}()

"Return a bijective map between domains `d1` and `d2`."
bijection(d1, d2) = bijection1(d1, d2)

# simplify the first argument
bijection1(d1, d2) = _bijection1(d1, d2, canonicaldomain(d1))
function _bijection1(d1, d2, cd)
    if d1 == cd
        bijection2(d1, d2)
    else
        bijection(cd, d2) ∘ tocanonical(d1)
    end
end
# simplify the second argument
bijection2(d1, d2) = _bijection2(d1, d2, canonicaldomain(d2))
function _bijection2(d1, d2, cd)
    if d2 == cd
        no_known_bijection(d1, d2)
    else
        fromcanonical(d2) ∘ bijection(d1, cd)
    end
end

no_known_bijection(d1, d2) = d1 == d2 ? IdentityMap{eltype(d1)}() : error("No bijection known between $(d1) and $(d2).")

"Return a parameterization of the given domain."
parameterization(d::Domain) = fromcanonical(d)
