# Definition of the abstract Domain type and its interface


spaceof(::Domain{T}) where {T} = spacetype(T)

eltype(::Type{<:Domain{T}}) where {T} = T

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

"""
A `VectorDomain` is any domain whose eltype is `Vector{T}`. In this case
the dimension of the domain is not included in its type.
"""
const VectorDomain{T} = Domain{Vector{T}}

# At the level of Domain we attempt to promote the arguments to compatible
# types, then we invoke indomain. Concrete subtypes should implement
# indomain. They may also implement `in` in order to accept more types.
in(x::S, d::Domain{T}) where {S,T} = _in(x, d, promote_type(S,T))
_in(x::T, domain::Domain{T}, ::Type{T}) where {T} = indomain(x, domain)
_in(x::S, domain::Domain{T}, ::Type{U}) where {S,T,U} = indomain(convert(U, x), convert(Domain{U}, domain))
_in(x::S, domain::Domain{T}, ::Type{S}) where {S,T} = indomain(x, convert(Domain{S}, domain))
_in(x::S, domain::Domain{T}, ::Type{T}) where {S,T} = indomain(convert(T, x), domain)

function _in(x::S, domain::Domain{T}, ::Type{Any}) where {S,T}
    @warn "in: incompatible types $(S) and $(T). Returning false."
    return false
end

# Special cases:
# - We allow any vector for Euclidean domains, with conversion to a static vector
in(x::AbstractVector{T}, domain::EuclideanDomain{N,T}) where {N,T} =
    indomain(convert(SVector{N,T}, x), domain)
in(x::AbstractVector{S}, domain::EuclideanDomain{N,T}) where {N,S,T} =
    in(convert(SVector{N,promote_type(T,S)}, x), convert(EuclideanDomain{N,promote_type(S,T)}, domain))

# - we allow any abstract vector for Vector domains
in(x::AbstractVector{T}, domain::VectorDomain{T}) where {T} = indomain(x, domain)
in(x::AbstractVector{S}, domain::VectorDomain{T}) where {S,T} =
   in(convert(AbstractVector{promote_type(S,T)}, x), convert(VectorDomain{promote_type(S,T)}, domain))

# - we allow tuples for Euclidean domains if they have the right length
in(x::NTuple{N,T}, domain::EuclideanDomain{N,T}) where {N,T} = indomain(x, domain)
in(x::NTuple{N,S}, domain::EuclideanDomain{N,T}) where {N,S,T} =
    in(convert(NTuple{N,promote_type(S,T)}, x), convert(EuclideanDomain{N,promote_type(S,T)}, domain))

function indomain(x, domain::Domain)
    @warn "Please implement `indomain` for domain $(domain) for points of type $(typeof(x)). Returning false."
    return false
end




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


approx_in(x, domain::Domain) = approx_in(x, domain, default_tolerance(domain))

# We mimick the promotion rules of `in` above.
approx_in(x::S, d::Domain{T}, tolerance) where {S,T} = _approx_in(x, d, tolerance, promote_type(S,T))
_approx_in(x::S, domain::Domain{T}, tolerance, ::Type{U}) where {S,T,U} =
    approx_indomain(convert(U, x), convert(Domain{U}, domain), tolerance)
_approx_in(x::S, domain::Domain{T}, tolerance, ::Type{S}) where {S,T} =
    approx_indomain(x, convert(Domain{S}, domain), tolerance)
_approx_in(x::S, domain::Domain{T}, tolerance, ::Type{T}) where {S,T} =
    approx_indomain(convert(T, x), domain, tolerance)

function _approx_in(x::S, domain::Domain{T}, tolerance, ::Type{Any}) where {S,T}
    @warn "in: incompatible types $(S) and $(T). Returning false."
    return false
end

# Special cases analogous to the ones above for `in` follow:
approx_in(x::AbstractVector{T}, domain::EuclideanDomain{N,T}, tolerance) where {N,T} =
    approx_indomain(convert(SVector{N,T}, x), domain, tolerance)
approx_in(x::AbstractVector{S}, domain::EuclideanDomain{N,T}, tolerance) where {N,S,T} =
    approx_in(convert(SVector{N,promote_type(T,S)}, x), convert(EuclideanDomain{N,promote_type(S,T)}, domain), tolerance)

approx_in(x::AbstractVector{T}, domain::VectorDomain{T}, tolerance) where {T} =
    approx_indomain(x, domain, tolerance)
approx_in(x::AbstractVector{S}, domain::VectorDomain{T}, tolerance) where {S,T} =
   approx_in(convert(AbstractVector{promote_type(S,T)}, x), convert(VectorDomain{promote_type(S,T)}, domain), tolerance)

approx_in(x::NTuple{N,T}, domain::EuclideanDomain{N,T}, tolerance) where {N,T} =
    approx_indomain(x, domain, tolerance)
approx_in(x::NTuple{N,S}, domain::EuclideanDomain{N,T}, tolerance) where {N,S,T} =
    approx_in(convert(NTuple{N,promote_type(S,T)}, x), convert(EuclideanDomain{N,promote_type(S,T)}, domain), tolerance)

function approx_indomain(x, domain::Domain, tolerance)
    @warn "Please consider implementing `approx_in` for domain $(domain) for points of type $(typeof(x))."
    # Fall back to in without tolerance
    return in(x, domain)
end



isapprox(d1::Domain, d2::Domain; kwds...) = d1 == d2

isreal(d::Domain) = isreal(spaceof(d))

infimum(d::Domain) = minimum(d)  # if the minimum exists, then it is also the infimum
supremum(d::Domain) = maximum(d)  # if the maximum exists, then it is also the supremum


# override minimum and maximum for closed sets
function boundary end
const ∂ = boundary
