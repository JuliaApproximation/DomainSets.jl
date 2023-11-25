# Definition of the abstract Domain type and its interface

Domain(d) = convert(Domain, d)

prectype(::Type{<:Domain{T}}) where T = prectype(T)
numtype(::Type{<:Domain{T}}) where T = numtype(T)

# Domain-specific prectype and numtype default to using domaineltype
domain_prectype(d) = prectype(domaineltype(d))
domain_numtype(d) = numtype(domaineltype(d))
prectype(d::DomainRef) = prectype(domaineltype(d))
numtype(d::DomainRef) = numtype(domaineltype(d))

# Concrete types can implement similardomain(d, ::Type{T}) where {T}
# to support convert(Domain{T}, d) functionality.
convert(::Type{Domain{T}}, d::Domain{T}) where {T} = d
convert(::Type{Domain{T}}, d::Domain{S}) where {S,T} = similardomain(d, T)

"Can the domains be promoted without throwing an error?"
promotable_domains(domains...) = promotable_eltypes(map(domaineltype, domains)...)
promotable_eltypes(types...) = isconcretetype(promote_type(types...))
promotable_eltypes(::Type{S}, ::Type{T}) where {S<:AbstractVector,T<:AbstractVector} =
    promotable_eltypes(eltype(S), eltype(T))

"Promote the given domains to have a common element type."
promote_domains() = ()
promote_domains(domains...) = promote_domains(domains)
promote_domains(domains) = convert_eltype.(mapreduce(domaineltype, promote_type, domains), domains)

convert_eltype(::Type{T}, d::Domain) where {T} = convert(Domain{T}, d)
convert_eltype(::Type{T}, d) where {T} = _convert_eltype(T, d, domaineltype(d), DomainStyle(d))
_convert_eltype(::Type{T}, d, ::Type{T}, ::IsDomain) where {T} = d
_convert_eltype(::Type{T}, d, ::Type{S}, ::IsDomain) where {S,T} =
    error("Don't know how to convert the `eltype` of $(d).")
_convert_eltype(::Type{T}, d, ::Type{T}, ::NotDomain) where {T} = d
_convert_eltype(::Type{T}, d, ::Type{S}, ::NotDomain) where {S,T} =
    error("Convert eltype: argument given is not a domain.")


"A `EuclideanDomain` is any domain whose eltype is `<:StaticVector{N,T}`."
const EuclideanDomain{N,T} = Domain{<:StaticVector{N,T}}

"A `VectorDomain` is any domain whose eltype is `Vector{T}`."
const VectorDomain{T} = Domain{Vector{T}}

"An `AbstractVectorDomain` is any domain whose eltype is `<:AbstractVector{T}`."
const AbstractVectorDomain{T} = Domain{<:AbstractVector{T}}

CompositeTypes.Display.displaysymbol(d::Domain) = 'D'

"What is the Euclidean dimension of the domain?"
dimension(d) = euclideandimension(domaineltype(d))

"Is the given combination of point and domain compatible?"
iscompatiblepair(x, d) = _iscompatiblepair(x, d, typeof(x), domaineltype(d))
_iscompatiblepair(x, d, ::Type{S}, ::Type{T}) where {S,T} =
    _iscompatiblepair(x, d, S, T, promote_type(S,T))
_iscompatiblepair(x, d, ::Type{S}, ::Type{T}, ::Type{U}) where {S,T,U} = true
_iscompatiblepair(x, d, ::Type{S}, ::Type{T}, ::Type{Any}) where {S,T} = false
_iscompatiblepair(x, d, ::Type{S}, ::Type{Any}, ::Type{Any}) where {S} = true

# Some generic cases where we can be sure:
iscompatiblepair(x::SVector{N}, ::EuclideanDomain{N}) where {N} = true
iscompatiblepair(x::SVector{N}, ::EuclideanDomain{M}) where {N,M} = false
iscompatiblepair(x::AbstractVector, ::EuclideanDomain{N}) where {N} = length(x)==N

# Note: there are cases where this warning reveals a bug, and cases where it is
# annoying. In cases where it is annoying, the domain may want to specialize `in`.
compatible_or_false(x, domain) =
    iscompatiblepair(x, domain) ? true : (@warn "`in`: incompatible combination of point: $(typeof(x)) and domain eltype: $(domaineltype(domain)). Returning false."; false)

compatible_or_false(x::AbstractVector, domain::AbstractVectorDomain) =
    iscompatiblepair(x, domain) ? true : (@warn "`in`: incompatible combination of vector with length $(length(x)) and domain '$(domain)' with dimension $(dimension(domain)). Returning false."; false)


"Promote point and domain to compatible types."
promote_pair(x, d) = _promote_pair(x, d, promote_type(typeof(x), domaineltype(d)))
_promote_pair(x, d, ::Type{T}) where {T} = convert(T, x), convert_eltype(T, d)
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
domain_tolerance(d) = domain_tolerance(domain_prectype(d))
domain_tolerance(::Type{T}) where {T <: AbstractFloat} = 100eps(T)

# a version with a tolerance, for use in approx_in
function compatible_or_false(x, d, tol)
    tol >= 0 || error("Tolerance has to be positive in `approx_in`.")
    compatible_or_false(x, d)
end

"""
`approx_in(x, domain::Domain [, tolerance])`

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
approx_in(x, d) = approx_in(x, d, domain_tolerance(d))
approx_in(x, d::Domain, tol) =
    compatible_or_false(x, d, tol) && approx_indomain(promote_pair(x, d)..., tol)
approx_in(x, d::DomainRef, tol) =
    compatible_or_false(x, domain(d), tol) && approx_indomain(promote_pair(x, domain(d))..., tol)

# Fallback to `in`
approx_indomain(x, d, tol) = in(x, d)

isapprox(d1::Domain, d2::Domain; kwds...) = d1 == d2

isreal(d::Domain) = isreal(eltype(d))

"Return a point from the given domain."
function choice(d) end

choice(d::AbstractSet) = first(d)
choice(d::AbstractArray) = first(d)
