# space.jl

"""
A domain is always a subset of a geometric space. GeometricSpace is the abstract
supertype of all possible spaces.

A geometric space has an element type, specified as type parameter `T`.
"""
abstract type GeometricSpace{T} end

eltype{T}(::GeometricSpace{T}) = T
eltype{T}(::Type{GeometricSpace{T}}) = T
eltype{S <: GeometricSpace}(::Type{S}) = eltype(supertype(S))

isreal(space::GeometricSpace) = isreal(eltype(space))

# Return the zero element
zero(space::GeometricSpace) = zero(eltype(space))

# Default definitions of element membership, based only on type
# - override if there are additional restrictions on x
in{T}(x::T, space::GeometricSpace{T}) = true

# - another type S can be in the space if S promotes to T
function in{S,T}(x::S, space::GeometricSpace{T})
    a,b = promote(x, zero(space))
    (typeof(a) == T) && (a ∈ space)
end
