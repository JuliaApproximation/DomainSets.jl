# domain.jl
# Definition of the abstract Domain type and its interface

"""
`Domain{T}` is the abstract supertype of all subsets of the geometric space
`GeometricSpace{T}`.

A domain is defined mathematically by its indicator or characteristic function.
This is implemented via `in`: one can write `x ∈ D` or `in(x,D)`, and this
evaluates to true or false.
"""
abstract type Domain{T}
end

spaceof(::Domain{T}) where {T} = spacetype(T)

eltype(::Type{Domain{T}}) where {T} = T
eltype(::Type{D}) where {D <: Domain} = eltype(supertype(D))

dimension(::Type{Domain{T}}) where {T} = dimension_type(T)
dimension(::Type{D}) where {D <: Domain} = dimension(supertype(D))
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


# Convenient aliases
const Domain1d{T <: Number} = Domain{T}
const Domain2d{T} = EuclideanDomain{2,T}
const Domain3d{T} = EuclideanDomain{3,T}
const Domain4d{T} = EuclideanDomain{4,T}


# We implement the indicator function by overriding `in`.
# The implementation of `in` at the level of Domain converts the point into
# the element type of the domain using promotion, and calls `indomain`.
# Concrete subtypes should implement `indomain`, rather than `in`.
in(x::T, d::Domain{T}) where {T} = indomain(x, d)

function in(x, d::Domain)
   T = promote_type(typeof(x), eltype(d))
   T == Any && return false
   in(convert(T, x), convert(Domain{T}, d))
end

# The user may supply a vector. We attempt to convert it to the right space.
in(x::Vector{T}, d::Domain) where {T} = in(convert(eltype(d), x), d)

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

approx_in(x::T, d::Domain{T}, tolerance = default_tolerance(d)) where {T} =
   approx_indomain(x, d, tolerance)

function approx_in(x, d::Domain, tolerance = default_tolerance(d))
   T = promote_type(typeof(x), eltype(d))
   T == Any && return false
   approx_in(convert(T, x), convert(Domain{T}, d), tolerance)
end


isreal(d::Domain) = isreal(spaceof(d))

infimum(d::Domain) = minimum(d)  # if the minimum exists, then it is also the infimum
supremum(d::Domain) = maximum(d)  # if the maximum exists, then it is also the supremum

# override minimum and maximum for closed sets

boundary(d::Domain) = ∂(d)
