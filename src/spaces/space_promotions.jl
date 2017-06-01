# space_promotions.jl

# Conversion

"""
Convert the variable `x` to an element of the given space.
"""
convert_space(::Type{GeometricSpace{T}}, x::T) where {T} = x

# We supply a default conversion using the Julia conversion between types
convert_space(::Type{GeometricSpace{T}}, x) where {T} = convert(T, x)

const True = Val{true}
const False = Val{false}

(&)(::Type{True}, ::Type{True}) = True
(&)(::Type{True}, ::Type{False}) = False
(&)(::Type{False}, ::Type{True}) = False
(&)(::Type{False}, ::Type{False}) = False

result(::Type{True}) = true
result(::Type{False}) = false

#############
# Embedding
#############

"""
The function `embedding_rule` describes in the type domain whether a space with
type `T` is embedded into a space with type `B`. If so, it returns True,
otherwise it returns False.

Note that True == Val{true} and False == Val{false} are types, not boolean values.
"""
embedding_rule(::Type{GeometricSpace{T}}, ::Type{GeometricSpace{T}}) where {T} = True

# Other spaces are embedded by default only if the type S promotes to T.
embedding_rule(::Type{GeometricSpace{S}}, ::Type{GeometricSpace{T}}) where {S,T} = _embedding_via_promotion(T,promote_type(S,T))

# - ok, S promotes to T
_embedding_via_promotion(::Type{T}, ::Type{T}) where {T} = True
# - not ok, S and T promote to something else: return False
_embedding_via_promotion(::Type{T}, ::Type{U}) where {T,U} = False

"""
Returns true if `A` is embedded in `B`.
"""
embedded(a::GeometricSpace, b::GeometricSpace) = embedded(typeof(a), typeof(b))
embedded(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} = result(embedding_rule(A, B))

"The symbol ↪ (\hookrightarrow) is a synonym for `embedded`."
↪ = embedded



##############
# Ismorphisms
##############

isomorphism_rule(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} =
    embedding_rule(A,B) & embedding_rule(B,A)

"""
Returns true if `A` is embedded in `B`.
"""
isomorphic(a::GeometricSpace, b::GeometricSpace) = isomorphic(typeof(a), typeof(b))
isomorphic(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} = result(isomorphism_rule(A, B))

"The symbol ≅ (\cong) is a synonym for `isomorphic`."
≅ = isomorphic
