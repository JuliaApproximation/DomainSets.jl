# space.jl

"""
A domain is always a subset of a geometric space. A geometric space is completely
characterized by the type `T` of its elements. There are no restrictions on
concrete values: each instance of type `T` belongs to `GeometricSpace{T}`.
Conversely, each numeric type `T` can give rise to only one geometric space.

The goal of `GeometricSpace` is to provide a lightweight framework for converting
elements between spaces. Two features are embeddings and isomorphisms. It is
lightweight because the space is characterized by the type and has no state, hence
it can be inferred from any variable based only on its type. A restriction that
results from this choice is that there can only be one type of embedding from one
space to another. More general embeddings could be implemented using domains or
specific maps.

Embeddings:
A space `A` is embedded in a space `B` if each element of `A` corresponds to an
element of `B`. In that case an element of `A` can be promoted (using `promote_space`)
to an element of type `B`. This is not necessarily the same as the promotion of
the type of `A` to the type of `B`. For example, a scalar can be embedded into a
two-dimensional space (say with the second component equal to zero). It would be
undesirable to implement that using the standard Julia promotion system, since
it would apply to all Julia code used concurrently.

Isomorphisms:
If `A` is embedded in `B` and `B` is embedded in `A` than they `A` and `B` are
isomorphic. In that case, an element of `A` can be converted to an element of
`B` and vice-versa. The conversion is invertible. One example is the isomorphism
between ℝ^2 and ℂ.

New spaces can be defined by defining a new numeric type.
"""
struct GeometricSpace{T}
end

eltype{T}(::GeometricSpace{T}) = T
eltype{T}(::Type{GeometricSpace{T}}) = T
eltype{S <: GeometricSpace}(::Type{S}) = eltype(supertype(S))

isreal{T}(space::GeometricSpace{T}) = isreal(T)

# Return the zero element
zero{T}(space::GeometricSpace{T}) = zero(T)

# Definition of element membership is based only on type:
# - x is in the space of the type of x equals that of the space
in{T}(x::T, space::GeometricSpace{T}) = true
# - or if its type can be promoted to that of the space
in{S,T}(x::S, space::GeometricSpace{T}) = promote_type(S,T) == T

# Make the space bigger by widening the numeric type
widen{T}(::GeometricSpace{T}) = GeometricSpace{widen(T)}()

"Return the space of all elements with the same type as `x`."
space(x) = space(typeof(x))
space(::Type{T}) where {T} = GeometricSpace{T}()
