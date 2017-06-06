# domain.jl
# Contains the definition of the abstract type Domain and its interface.


"""
`Domain{T}` is the abstract supertype of all subsets of the geometric space
`GeometricSpace{T}`.

Memberschip of a subset can be determined in several ways, depending on the
type of the domain. For many domains it is implemented via the indicator or
characteristic function, which overrides `in`: one can write `x ∈ D` or
`in(x,D)`, and this evaluates to true or false.
"""
abstract type Domain{T}
end

spaceof(::Domain{T}) where {T} = spacetype(T)

eltype(::Type{Domain{T}}) where {T} = T
eltype(::Type{D}) where {D <: Domain} = eltype(supertype(D))

"We use `Point{N,T}` as a synonym for `SVector{N,T}`."
const Point{N,T} = SVector{N,T}

const EuclideanDomain{N,T} = Domain{Point{N,T}}

ndims(::Type{Domain{T}}) where {T} = ndims_type(T)
ndims(::Type{D}) where {D <: Domain} = ndims(supertype(D))
ndims(d::Domain) = ndims(typeof(d))
ndims_type(::Type{SVector{N,T}}) where {N,T} = N
ndims_type(::Type{T}) where {T <: Number} = 1



# Convenient aliases
const Domain1d{T} = EuclideanDomain{1,T}
const Domain2d{T} = EuclideanDomain{2,T}
const Domain3d{T} = EuclideanDomain{3,T}
const Domain4d{T} = EuclideanDomain{4,T}

# left and right of domains falls back to bounding box domains
# TODO: these should go away
left(d::Domain) = left(boundingbox(d))
right(d::Domain) = right(boundingbox(d))

left(d::Domain, i::Int) = left(boundingbox(d), i)
right(d::Domain, i::Int) = right(boundingbox(d), i)

# We implement the indicator function by overriding `in`.
# The implementation of `in` at the level of Domain converts the point into
# the element type of the domain using promotion, and calls `indomain`.
# Concrete subtypes should implement `indomain`, rather than `in`.
in(x::T, d::Domain{T}) where {T} = indomain(x, d)
in(x::S, d::Domain{T}) where {T,S} = in(convert(T, x), d)

# The user may supply a vector. We attempt to convert it to the right space.
in(x::Vector{T}, d::Domain) where {T} = in(convert(SVector{ndims(typeof(d)),T}, x), d)

# Check whether a value is in an interval, up to 10 times machine precision
in(x::Number, a::T, b::T) where {T <: AbstractFloat} = (a-10eps(T) <= x <= b+10eps(T))
in(x::Number, a::T, b::T) where {T <: Number} = a <= x <= b

# Intercept a broadcasted call to indomain. We assume that the user wants evaluation
# in a set of points (which we call a grid), rather than in a single point.
broadcast(::typeof(in), grid, d::Domain) = indomain_broadcast(grid, d)

# # Default methods for evaluation on a grid: the default is to call eval on the domain with
# # points as arguments. Domains that have faster grid evaluation routines may define their own version.
indomain_broadcast(grid, d::Domain) = indomain_broadcast!(zeros(Bool, size(grid)), grid, d)
# TODO: use BitArray here

# Note that indomain_broadcast! only updates the result - it should be initialized to all false!
# The idea is that you can chain different calls to indomain_broadcast (as used in DomainCollection)
# TODO: that was a bad idea, change
function indomain_broadcast!(result, grid, domain::Domain)
    for (i,x) in enumerate(grid)
        result[i] |= indomain(x, domain)
    end
    result
end
