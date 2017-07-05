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

"""
If the type `T` is a container type, the elements of `T` may have a different
`subeltype`. If `T` is not a container, `subeltype` simply evaluates to `T`.
"""
subeltype(d) = subeltype(spaceof(d))

"We use `Point{N,T}` as a synonym for `SVector{N,T}`."
const Point{N,T} = SVector{N,T}

"A `EuclideanDomain` is any domain whose eltype is `Point{N,T}`."
const EuclideanDomain{N,T} = Domain{Point{N,T}}


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

# TODO should we genererate an error or just False if we ask whether 1.5 lays in FullSpace(Int)
# If we want error, remove try catch.
function in(x::S, d::Domain{T}) where {T,S}
  try
    in(convert(T, x), d)
  catch
    false
  end
end

# Be forgiving for a 1D case. Good idea or not? This is really an isomorphism.
in(x::SVector{1,T}, d::Domain{T}) where {T <: Number}  = in(x[1], d)

# The user may supply a vector. We attempt to convert it to the right space.
in(x::Vector{T}, d::Domain) where {T}  = in(convert(SVector{length(x),T}, x), d)

isreal(d::Domain) = isreal(spaceof(d))
