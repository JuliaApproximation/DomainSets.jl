# space_promotions.jl


################
# Preliminaries
################

const True = Val{true}
const False = Val{false}

(&)(::Type{True}, ::Type{True}) = True
(&)(::Type{True}, ::Type{False}) = False
(&)(::Type{False}, ::Type{True}) = False
(&)(::Type{False}, ::Type{False}) = False

result(::Type{True}) = true
result(::Type{False}) = false

"""
Return the superspace of the given geometric space. If a space has a superspace,
it can be embedded into that superspace. This is used to determine whether two
spaces can be promoted to a joined superspace.

The largest superspace is `AnySpace` and this is the default for any `T`.
"""
superspace(::Type{GeometricSpace{T}}) where {T} = AnySpace


#############
# Conversion
#############

"""
Convert the variable `x` to an element of the given space.
"""
convert_space(::Type{GeometricSpace{T}}, x::T) where {T} = x

# We supply a default conversion using the Julia conversion between types, after
# converting x to a superspace if needed.
convert_space(A::Type{GeometricSpace{T}}, x::S) where {T,S} = _convert_space(A, x, promote_type(T,S), superspace(spaceof(x)))

# - type of x promotes to T: simply convert to T
function _convert_space(A::Type{GeometricSpace{T}}, x::S, ::Type{T}, B::Type{GeometricSpace{V}}) where {T,S,V}
    convert(T,x)
end
# - type of x does not promote to T: convert x to its superspace
function _convert_space(A::Type{GeometricSpace{T}}, x::S, ::Type{U}, B::Type{GeometricSpace{V}}) where {T,S,U,V}
    convert_space(A, convert_space(B, x))
end


#############
# Embeddings
#############


"""
The function `embedding` describes in the type domain whether a space with
type `T` is embedded into a space with type `B`. If so, it returns True,
otherwise it returns False.

Note that True == Val{true} and False == Val{false} are types, not boolean values.
"""
embedding(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{T}}) where {T} = True

# Other spaces are embedded by default only if the type S promotes to T.
embedding(A::Type{GeometricSpace{S}}, B::Type{GeometricSpace{T}}) where {S,T} = _embedding(A, B, promotes_to(S,T))

# Check whether S promotes to T. If so, return True, otherwise return False.
promotes_to(::Type{S}, ::Type{T}) where {S,T} = _promotes_to(S, T, promote_type(S,T))
_promotes_to(::Type{S}, ::Type{T}, ::Type{T}) where {S,T} = True
_promotes_to(::Type{S}, ::Type{T}, ::Type{U}) where {S,T,U} = False

_embedding(A::Type{GeometricSpace{S}}, B::Type{GeometricSpace{T}}, ::Type{True}) where {S,T} = True
_embedding(A::Type{GeometricSpace{S}}, B::Type{GeometricSpace{T}}, ::Type{False}) where {S,T} =
    embedding(superspace(A), B)

# We add the rules below to make sure that the recursion above ends
embedding(::Type{AnySpace}, ::Type{AnySpace}) = True
embedding(::Type{AnySpace}, ::Type{GeometricSpace{T}}) where {T} = False


"""
Returns true if `A` is embedded in `B`.
"""
embedded(a::GeometricSpace, b::GeometricSpace) = embedded(typeof(a), typeof(b))
embedded(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} = result(embedding(A, B))

"The symbol ↪ (\hookrightarrow) is a synonym for `embedded`."
↪ = embedded



###############
# Isomorphisms
###############

isomorphism(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} =
    embedding(A,B) & embedding(B,A)

"""
Returns true if `A` is isomorphic to `B`.
"""
isomorphic(a::GeometricSpace, b::GeometricSpace) = isomorphic(typeof(a), typeof(b))
isomorphic(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} = result(isomorphism(A, B))

"The symbol ≅ (\cong) is a synonym for `isomorphic`."
≅ = isomorphic


#############
# Promotions
#############

# We mimick Julia's promotion system, with the key difference that the types
# a space can be promoted to are not actual supertypes, but superspaces.

# Identical spaces need no promotion
promote_space_rule(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{T}}) where {T} = GeometricSpace{T}

# Spaces are different: check for possible embedding
promote_space_rule(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} = _promote_space_rule(A, B, embedding(A,B), embedding(B,A))

# - There are no embeddings
_promote_space_rule(A, B, ::Type{False}, ::Type{False}) = BottomSpace
# - There is one embedding: choose that one
_promote_space_rule(A, B, ::Type{True}, ::Type{False}) = B
_promote_space_rule(A, B, ::Type{False}, ::Type{True}) = A
# - There is a bidirectional embedding: invoke the isomorphism_promotion_rule
_promote_space_rule(A, B, ::Type{True}, ::Type{True}) = isomorphism_promotion_rule(A, B)

isomorphism_promotion_rule(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S} =
    error("Geometric spaces ", A, " and ", B, " are isomophic: implement isomorphism_promotion_rule(A,B) to select a preferred promotion space.")

promote_space() = ()
promote_space(x) = (x,)
function promote_space(x::T, y::S) where {T,S}
    (convert_space(promote_space_type(spaceof(x),spaceof(y)), x), convert_space(promote_space_type(spaceof(x),spaceof(y)), y))
end

promote_space_type() = ()
promote_space_type(A::Type{GeometricSpace{T}}) where {T} = A
promote_space_type(A, B, C, D...) = promote_space_type(A, promote_space_type(B, C, D...))

# We can stop as soon as there is a BottomSpace
promote_space_type(::Type{BottomSpace}, ::Type{BottomSpace}) = AnySpace
promote_space_type(::Type{GeometricSpace{T}}, ::Type{BottomSpace}) where {T} = GeometricSpace{T}
promote_space_type(::Type{BottomSpace}, ::Type{GeometricSpace{T}}) where {T} = GeometricSpace{T}
promote_space_type(::Type{GeometricSpace{T}}, ::Type{GeometricSpace{T}}) where {T} = GeometricSpace{T}

# We can also stop as soon as there is an AnySpace
promote_space_type(::Type{AnySpace}, ::Type{AnySpace}) = AnySpace
promote_space_type(::Type{AnySpace}, ::Type{BottomSpace}) = AnySpace
promote_space_type(::Type{AnySpace}, ::Type{GeometricSpace{T}}) where {T} = AnySpace
promote_space_type(::Type{BottomSpace}, ::Type{AnySpace}) = AnySpace
promote_space_type(::Type{GeometricSpace{T}}, ::Type{AnySpace}) where {T} = AnySpace


# Logic of promote_space_type(A,B) for different A and B:
# - we first try the corresponding promote_rule of A and B
# - if that succeeds, we end by one of the above lines
# - if that fails
function promote_space_type(A::Type{GeometricSpace{T}}, B::Type{GeometricSpace{S}}) where {T,S}
    # Here we invoke promote_space_rule in both directions
    promote_space_result(A, B, promote_space_rule(A,B), promote_space_rule(B,A))
end

# If the promote_rule succeeds in either direction, then this line will end it
promote_space_result(A, B, C, D) = promote_type(C, D)

# If the promote rule failed, we have two BottomSpace's and we try again using the supertype
function promote_space_result(A, B, ::Type{BottomSpace}, ::Type{BottomSpace})
    promote_superspace_result(A, B, promote_space_type(superspace(A), B), promote_space_type(A, superspace(B)))
end

promote_superspace_result(A, B, ::Type{AnySpace}, ::Type{AnySpace}) = AnySpace
promote_superspace_result(A, B, ::Type{GeometricSpace{T}}, ::Type{AnySpace}) where {T} = GeometricSpace{T}
promote_superspace_result(A, B, ::Type{AnySpace}, ::Type{GeometricSpace{T}}) where {T} = GeometricSpace{T}
promote_superspace_result(A, B, ::Type{GeometricSpace{T}}, ::Type{GeometricSpace{S}}) where {T,S} =
    error("Promotion ambiguity for ", A, " and ", B)
