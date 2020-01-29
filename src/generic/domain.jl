# Definition of the abstract Domain type and its interface


spaceof(::Domain{T}) where {T} = spacetype(T)

eltype(::Type{<:Domain{T}}) where {T} = T

# TODO: this should go elsewhere, perhaps in common.jl
dimension(::Type{<:Domain{T}}) where {T} = dimension_type(T)
dimension(d::Domain) = dimension(typeof(d))
dimension_type(::Type{SVector{N,T}}) where {N,T} = N
dimension_type(::Type{T}) where {T <: Number} = 1


"""
If the type `T` is a container type, the elements of `T` may have a different
`subeltype`. If `T` is not a container, `subeltype` simply evaluates to `T`.
"""
subeltype(d::Domain) = subeltype(spaceof(d))
subeltype(::Type{T}) where {T} = subeltype(GSpace{T})

"A `EuclideanDomain` is any domain whose eltype is `SVector{N,T}`."
const EuclideanDomain{N,T} = Domain{SVector{N,T}}


"NoDomain is a placeholder for the absence of a suitable domain."
struct NoDomain{T} <: Domain{T} end

# At the level of Domain we attempt to promote the arguments to compatible
# types, then we invoke indomain. Concrete subtypes should implement
# indomain, and they can assume the types of the point and the domain compatible.
in(x, d::Domain) = indomain(point_domain_promote(x, d)...)
# Types of x and domain are incompatible: we return false
indomain(x, domain::NoDomain) = false

"""
Promote a point and a domain to compatible types by promoting `x` and the
element type of the domain to a joined supertype.
"""
point_domain_promote(x::T, domain::Domain{T}) where {T} = x, domain
point_domain_promote(x, domain::Domain) =
   _point_domain_promote(x, domain, promote_type(typeof(x), eltype(domain)))

_point_domain_promote(x, domain::Domain, ::Type{T}) where {T} =
   convert(T, x), convert(Domain{T}, domain)
_point_domain_promote(x, domain::Domain, ::Type{Any}) = x, NoDomain{eltype(domain)}()


"""
Return a suitable tolerance to use for verifying whether a point is close to
a domain. Typically, the tolerance is close to the precision limit of the numeric
type associated with the domain.
"""
default_tolerance(d::Domain) = default_tolerance(d, subeltype(d))
default_tolerance(d::Domain, ::Type{T}) where {T <: Real} = 100eps(T)
default_tolerance(d::Domain, ::Type{Complex{T}}) where {T <: Real} = 100eps(T)
# Default tolerance for integers and for any other type is zero
default_tolerance(d::Domain, ::Type{T}) where {T <: Integer} = zero(T)
default_tolerance(d::Domain, ::Type{T}) where {T} = zero(T)

approx_in(x, domain::Domain, tolerance = default_tolerance(domain)) =
   approx_indomain(point_domain_promote(x, domain)..., tolerance)
approx_indomain(x, domain::NoDomain, tolerance) = false


isapprox(d1::Domain, d2::Domain; kwds...) = d1 == d2

isreal(d::Domain) = isreal(spaceof(d))

infimum(d::Domain) = minimum(d)  # if the minimum exists, then it is also the infimum
supremum(d::Domain) = maximum(d)  # if the maximum exists, then it is also the supremum


# override minimum and maximum for closed sets
function boundary end
const ∂ = boundary
