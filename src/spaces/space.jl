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

New spaces can be defined by defining a new numeric type. Embeddings and
isomorphisms are defined by defining new `embedding_reduction` and
`isomorphism_reduction` rules, along with conversiong using `convert_space`
and `restrict_space`. See the documentation of these functions for
information on how to use them.
"""
struct GeometricSpace{T}
end

# An internal shorthand in order to avoid excessively long lines when there
# are several arguments of type GeometricSpace :-)
# We export a longer name in order to avoid name clashes with other Space types
const GSpace = GeometricSpace

"AnySpace is the superset of all possible geometric spaces."
const AnySpace = GeometricSpace{Any}

eltype(::GeometricSpace{T}) where {T} = T
eltype(::Type{GeometricSpace{T}}) where {T} = T
eltype(::Type{S}) where {S <: GeometricSpace} = eltype(supertype(S))

"""
Some geometric spaces have an eltype that is composed in terms of an underlying
eltype. For example, `Complex{Float64}` and `SVector{2,Float64}` are based on
`Float64`, which is their `subeltype`.
"""
subeltype(::Type{GeometricSpace{T}}) where {T} = T

"""
Create a space that is similar to the given space, but with a different `subeltype`.

For example, `similar_space(ComplexSpace{Float64}, BigFloat)` yields a
`ComplexSpace{BigFloat}`.
"""
similar_space(::Type{GeometricSpace{T}}, ::Type{S}) where {T,S} = GeometricSpace{S}


isreal(space::GeometricSpace{T}) where {T} = isreal(T)
isreal(::Type{GeometricSpace{T}}) where {T} = isreal(T)

# Return the zero element
zero(space::GeometricSpace{T}) where {T} = zero(T)
zero(::Type{GeometricSpace{T}}) where {T} = zero(T)

"The origin of a space is its zero element."
origin(space::GeometricSpace) = zero(space)

# Definition of element membership is strictly based on type:
in(x, ::GeometricSpace) = false
in(x, ::Type{GeometricSpace{T}}) where {T} = false
in(x::T, ::Type{GeometricSpace{T}}) where {T} = true


# Make the space bigger by widening the numeric type
widen(A::Type{GeometricSpace{T}}) where {T} = _widen(A, eltype(A), subeltype(A))
# eltype and subeltype are the same: widen those
_widen(A, ::Type{T}, ::Type{T}) where {T} = GeometricSpace{widen(T)}
# eltype and subeltype are different: widen the subeltype
_widen(A, ::Type{T}, ::Type{S}) where {T,S} = similar_space(A, widen(S))

"Return the geometric space type with eltype `T`."
spacetype(::Type{T}) where {T} = GeometricSpace{T}

"Return the geometric space of all elements with the same type as `x`."
spaceof(x::T) where {T} = spacetype(T)
# For convenience
spaceof(x::GeometricSpace{T}) where {T} = typeof(x)


"""
Return the superspace of the given geometric space `A`. The superspace of `A`
should be larger than `A`, and `A` should be embedded in it.

Superspaces are used to automatically discover embeddings and promotion rules.
Its role is analogous to that of `supertype` in the Julia type system. Indeed,
`AnySpace` is a superspace of all spaces, much like `Any` is a supertype of
all types.
"""
superspace(::Type{GeometricSpace{T}}) where {T} = AnySpace

superspaceof(x) = superspace(spaceof(x))

"""
A space `A` is a subspace of space `B` if `B` is a supertype of `A`.
"""  # Any space is a subspace of itself
issubspace(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{T}}) where {T} = true
# So is AnySpace. This is covered by the line above, but we make it explicit for clarity
issubspace(A::Type{AnySpace}, B::Type{AnySpace}) = true
# AnySpace is not a subspace of any other space
issubspace(A::Type{AnySpace}, B::Type{GeometricSpace{T}}) where {T} = false
# In all other cases, we replace A by its superspace and recurse
issubspace(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} = issubspace(superspace(A), B)
# TODO add issubspace to specific spaces
