# Definition of the abstract Domain type and its interface

eltype(::Type{<:Domain{T}}) where {T} = T
prectype(::Type{<:Domain{T}}) where {T} = prectype(T)
numtype(::Type{<:Domain{T}}) where {T} = numtype(T)

convert_numtype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{convert_numtype(T,U)}, d)
convert_prectype(d::Domain{T}, ::Type{U}) where {T,U} = convert(Domain{convert_prectype(T,U)}, d)

Domain(d) = convert(Domain, d)

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

# - We allow any vector for Euclidean domains
function _in(x::AbstractVector{T}, domain::EuclideanDomain{N,T}) where {N,T}
    # Avoid an error in an attempt to convert to SVector{N,T}
    if length(x) != N
        @warn "in(x,domain): incompatible dimension $(length(x)) of x and $(N) of the domain. Returning false."
        false
    else
        indomain(SVector{N,T}(x), domain)
    end
end
_in(x::AbstractVector{S}, domain::EuclideanDomain{N,T}) where {N,S,T} =
    _in(convert(AbstractVector{promote_type(S,T)}, x), convert(Domain{SVector{N,promote_type(S,T)}}, domain))

# - We allow any abstract vector for Vector domains
_in(x::AbstractVector{T}, domain::VectorDomain{T}) where {T} = indomain(x, domain)
_in(x::AbstractVector{S}, domain::VectorDomain{T}) where {S,T} =
    in(convert(AbstractVector{promote_type(S,T)}, x), convert(VectorDomain{promote_type(S,T)}, domain))

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


approx_in(x, domain::Domain) = approx_in(x, domain, default_tolerance(domain))

# We mimick the promotion rules of `in` above.
approx_in(x, d::Domain, tolerance) = _approx_in(x, d, tolerance)
_approx_in(x::S, d::Domain{T}, tolerance) where {S,T} = _approx_in(x, d, tolerance, promote_type(S,T))
_approx_in(x::S, domain::Domain{T}, tolerance, ::Type{U}) where {S,T,U} =
    approx_indomain(convert(U, x), convert(Domain{U}, domain), tolerance)
_approx_in(x::T, domain::Domain{T}, tolerance, ::Type{T}) where {T} =
    approx_indomain(x, domain, tolerance)

function _approx_in(x::S, domain::Domain{T}, tolerance, ::Type{Any}) where {S,T}
    @warn "in: incompatible types $(S) and $(T). Returning false."
    return false
end

# Special cases analogous to the ones above for `in` follow:
_approx_in(x::AbstractVector, domain::EuclideanDomain, tolerance) =
    approx_in(SVector(x...), domain, tolerance)
_approx_in(x::SVector{N,T}, domain::EuclideanDomain{N,T}, tolerance) where {N,T} =
    approx_indomain(x, domain, tolerance)
_approx_in(x::SVector{N,S}, domain::EuclideanDomain{N,T}, tolerance) where {N,S,T} =
    approx_indomain(convert(SVector{N,promote_type(S,T)}, x), convert(EuclideanDomain{N,promote_type(S,T)}, domain), tolerance)
function _approx_in(x::SVector{N,T}, domain::EuclideanDomain{M,S}, tolerance) where {N,M,S,T}
    @warn "in(x,domain): incompatible dimension $(N) of x and $(M) of the domain. Returning false."
    false
end

_approx_in(x::AbstractVector{T}, domain::VectorDomain{T}, tolerance) where {T} =
    approx_indomain(x, domain, tolerance)
_approx_in(x::AbstractVector{S}, domain::VectorDomain{T}, tolerance) where {S,T} =
   approx_in(convert(AbstractVector{promote_type(S,T)}, x), convert(VectorDomain{promote_type(S,T)}, domain), tolerance)

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
