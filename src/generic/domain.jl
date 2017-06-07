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

"We use `Point{N,T}` as a synonym for `SVector{N,T}`."
const Point{N,T} = SVector{N,T}

"A `EuclideanDomain` is any domain whose eltype is `Point{N,T}`."
const EuclideanDomain{N,T} = Domain{Point{N,T}}

# I don't like this definition of ndims. Perhaps there should be no ndims at all.
ndims(::Type{Domain{T}}) where {T} = ndims_type(T)
ndims(::Type{D}) where {D <: Domain} = ndims(supertype(D))
ndims(d::Domain) = ndims(typeof(d))
ndims_type(::Type{SVector{N,T}}) where {N,T} = N
ndims_type(::Type{T}) where {T <: Number} = 1


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
in(x::S, d::Domain{T}) where {T,S} = in(convert(T, x), d)

# Be forgiving for a 1D case. Good idea or not? This is really an isomorphism.
in(x::SVector{1,T}, d::Domain{T}) where {T <: Number} = in(x[1], d)

# The user may supply a vector. We attempt to convert it to the right space.
in(x::Vector{T}, d::Domain) where {T} = in(convert(SVector{ndims(typeof(d)),T}, x), d)
