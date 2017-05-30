# domain.jl
# Contains the definition of the abstract type Domain and its interface.

"""
A domain is a subset of N-dimensional Euclidean space.

The most important functionality of a domain is its indicator function: the
function that decides whether or not a point belongs to the domain.
This is implemented by overriding `in`: e.g. one can write `x ∈ D` and this evaluates
to true or false.
"""
abstract type Domain{N} end

ndims{N}(::Type{Domain{N}}) = N
ndims{D <: Domain}(::Type{D}) = ndims(supertype(D))
ndims{N}(::Domain{N}) = N

# Convenient aliases
Domain1d = Domain{1}
Domain2d = Domain{2}
Domain3d = Domain{3}
Domain4d = Domain{4}

# left and right of domains falls back to bounding box domains
# TODO: these should go away
left(d::Domain) = left(boundingbox(d))
right(d::Domain) = right(boundingbox(d))

left(d::Domain, i::Int) = left(boundingbox(d), i)
right(d::Domain, i::Int) = right(boundingbox(d), i)

# We implement the indicator function by overriding `in`.
# The implementation of `in` at the level of Domain converts the point into
# an SVector, and calls `indomain`. The latter can be implemented by concrete
# domains, with fewer chance of ambiguities than when directly implementing `in`.
in{N}(x::SVector{N}, d::Domain{N}) = indomain(x, d)
in(x::Number, d::Domain1d) = indomain(x, d)
in(x::SVector{1}, d::Domain1d) = indomain(x[1], d)

# Convert a point given as any other vector into a SVector. This will throw an
# error if the dimension of x was wrong.
in{N}(x::AbstractVector, d::Domain{N}) = in(SVector{N}(x), d)

# Check whether a value is in an interval, up to 10 times machine precision
in{T <: AbstractFloat}(x::Number, a::T, b::T) = (a-10eps(T) <= x <= b+10eps(T))
in{T <: Number}(x::Number, a::T, b::T) = a <= x <= b

# Intercept a broadcasted call to indomain. We assume that the user wants evaluation
# in a set of points (a grid), rather than in a single point.
broadcast(::typeof(in), grid, d::Domain) = indomain_grid(grid, d)

# # Default methods for evaluation on a grid: the default is to call eval on the domain with
# # points as arguments. Domains that have faster grid evaluation routines may define their own version.
indomain_grid(grid, d::Domain) = indomain_grid!(zeros(Bool, size(grid)), grid, d)

# Note that indomain_grid! only updates the result - it should be initialized to all false!
# The idea is that you can chain different calls to indomain_grid (as used in DomainCollection)
# TODO: that was a bad idea, change
function indomain_grid!(result, grid, domain::Domain)
    for (i,x) in enumerate(grid)
        result[i] |= indomain(x, domain)
    end
    result
end
