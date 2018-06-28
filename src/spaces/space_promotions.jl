# space_promotions.jl


# Find the first space in a list of arguments that is not AnySpace
lcd(::Type{GSpace{T}}) where {T} = GSpace{T}
lcd(::Type{AnySpace}, ::Type{AnySpace}) = AnySpace
lcd(::Type{AnySpace}, ::Type{GSpace{T}}) where {T} = GSpace{T}
lcd(::Type{GSpace{T}}, ::Type{AnySpace}) where {T} = GSpace{T}
lcd(::Type{GSpace{T}}, ::Type{GSpace{S}}) where {T,S} = GSpace{T}
lcd(A::Type{GSpace{T}}, B::Type{GSpace{S}}, C::Type{GSpace{U}}, D...) where {T,S,U} = lcd(A, lcd(B, C, D...))



###############
# Isomorphisms
###############

"""
Return `True` if `A` and `B` are isomorphic. This implies that any point in `A`
can be converted to a point in `B`, and vice-versa, with an exact inverse.

The functionality is implemented recursively in terms of the
`isomorphism_reduction` function.
"""
isomorphism(::Type{GSpace{T}}, ::Type{GSpace{T}}) where {T} = True

isomorphism(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} =
    isomorphism_reduction_result(A, B, isomorphism_reduction(A,B)..., isomorphism_reduction(B,A)...)

"""
Indicate a reduction rule for isomorphism. The function either returns `()` or
a 2-tuple. If `isomorphism_reduction(A,B)` returns `(C,D)`, this means that `A`
and `B` are isomorphic if `C` and `D` are. Here, `C` and `D` are the subeltypes
of `A` and `B` respectively.

For example, `VectorSpace{2,T}` and `ComplexSpace{S}` are isomorphic if
`GeometricSpace{T}` and `GeometricSpace{S}` are isomorphic.

The isomorphism only has to be declared in one direction, and is automatically
valid in both directions. For space promotions, in case of ambiguity the first
space `A` here is chosen as the preferred promotion space.
"""
isomorphism_reduction(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} = ()

isomorphism_reduction_result(A, B) = False
isomorphism_reduction_result(A, B, C, D) = isomorphism(C, D)
isomorphism_reduction_result(A, B, C, D, E, F) = _duplicate_reduction_result(A, B, C, D, E, F)
# We have a situation here: isomorphism_reduction leads to a result in both
# directions.This is fine, as long as the outcome is the same. This actually
# happens for composite types, because they may differ only in their subtypes.
# Return an error if the results are not consistent, because that may never be
# the case.
_duplicate_reduction_result(A, B, C::Type{T}, D::Type{T}, E::Type{T}, F::Type{T}) where {T} = isomorphism(C, D)
_duplicate_reduction_result(A, B, C::Type{T}, D::Type{S}, E::Type{T}, F::Type{S}) where {T,S} = isomorphism(C, D)
_duplicate_reduction_result(A, B, C::Type{T}, D::Type{S}, E::Type{S}, F::Type{T}) where {T,S} = isomorphism(C, D)
_duplicate_reduction_result(A, B, C, D, E, F) = error("Duplicate isomorphism_reduction defined for ", A, " and ", B, ", but with inconsistent results.")

"""
Returns true if `A` is isomorphic to `B`, and false otherwise.
"""
isomorphic(a::GSpace, b::GSpace) = isomorphic(typeof(a), typeof(b))
isomorphic(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} = result(isomorphism(A, B))

"The symbol ≅ (\\cong) is a synonym for `isomorphic`."
≅ = isomorphic




#############
# Embeddings
#############


"""
The function `embedding` describes in the type domain whether a space with
type `T` is embedded into a space with type `B`. If so, it returns True,
otherwise it returns False.

Embeddings are the result of the following rules:
1) A{T} is embedded in A{S} if T promotes to S in the Julia type system
2) If a rule embedding_reduction(A,B) = (C,D) has been defined, then A is
   embedded in B if C is embedded in D.
3) Say A and B are isomorphic if C and D are. This automatically results in
   embedding reduction rules from A to B and from B to A.
4) A is embedded in B if the superspace of A is embedded in B.
"""
embedding(::Type{GeometricSpace{T}}, ::Type{GeometricSpace{T}}) where {T} = True

# Rules for AnySpace, necessary for termination of the recursive procedure using
# superspace in what follows:
# All spaces are embedded in AnySpace, but AnySpace is only embedded in itself
embedding(::Type{GeometricSpace{T}}, ::Type{AnySpace}) where {T} = True
embedding(::Type{AnySpace}, ::Type{AnySpace}) = True
embedding(::Type{AnySpace}, ::Type{GeometricSpace{T}}) where {T} = False

# We try all different cases. For isomorphisms, we have to check isomorphism_reduction
# with both possible orders of the arguments (A,B) and (B,A).
embedding(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} =
    one_of(
        embedding_via_promotion(A, B),
        embedding_via_reduction(A, B),
        embedding_via_superspace(A, B))


embedding_via_promotion(A, B) =
    _embedding_via_promotion(A, B, promote_type(eltype(A),eltype(B)))
_embedding_via_promotion(::Type{GSpace{T}}, ::Type{GSpace{S}}, ::Type{S}) where {T,S} = True
_embedding_via_promotion(::Type{GSpace{T}}, ::Type{GSpace{S}}, ::Type{U}) where {T,S,U} = False

# By default, two spaces have no embedding reduction rules.
# However, if there is an isomorphic_reduction rule, then there is automatically
# an embedding_reduction possible in both directions.
embedding_reduction(::Type{GSpace{T}}, ::Type{GSpace{T}}) where {T} = ()
embedding_reduction(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} =
   _embedding_reduction(A, B, isomorphism_reduction(A,B)..., 1, isomorphism_reduction(B,A)...)
_embedding_reduction(A, B, ::Int) = ()
_embedding_reduction(A, B, C, D, ::Int) = (C,D)
_embedding_reduction(A, B, ::Int, E, F) = (F,E)
_embedding_reduction(A, B, C, D, ::Int, E, F) = (C,D)

embedding_via_reduction(A, B) = _embedding_via_reduction(A, B, embedding_reduction(A, B)...)
# There was no specific rule for A and B
_embedding_via_reduction(A, B) = False
# A is embedded in B if C is embedded in D
_embedding_via_reduction(A, B, C, D) = embedding(C, D)

embedding_via_superspace(A, B) = embedding(superspace(A), B)



"""
Returns true if `A` is embedded in `B`.
"""
embedded(a::GSpace, b::GSpace) = embedded(typeof(a), typeof(b))
embedded(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} = result(embedding(A, B))

"The symbol ↪ (\\hookrightarrow) is a synonym for `embedded`."
↪ = embedded


#############
# Conversion
#############

"""
Convert the variable `x` to an element of the space `B`. This is possible if
the space of `x` is embedded in `B`.
"""# We don't need to do anything if the space of x is B
convert_space(B::Type{GSpace{T}}, x::T) where {T} = x

# If it isn't, then we dispatch on the type of embedding
convert_space(B::Type{GSpace{T}}, x) where {T} = convert_spaces(x, spaceof(x), B)

convert_spaces(x, A, B) = convert_spaces(x, A, B,
    embedding_via_promotion(A, B),
    embedding_via_reduction(A, B),
    embedding_via_superspace(A, B))

convert_spaces(x, A, B, ::Type{True}, d2, d3) =
    convert_space_via_promotion(x, A, B)
convert_spaces(x, A, B, ::Type{False}, ::Type{True}, d3) =
    convert_space_via_reduction(x, A, B)
convert_spaces(x, A, B, ::Type{False}, ::Type{False}, ::Type{True}) =
    convert_space_via_superspace(x, A, B)
convert_spaces(x, A, B, ::Type{False}, ::Type{False}, ::Type{False}) =
    throw(InexactError(:convert_spaces, B, x))


# Embedding via promotion: promote the type of x using convert
convert_space_via_promotion(x, A, B) = convert(eltype(B), x)

# Embedding via reduction: if A is embedded in B if C is embedded in D,
# then we have to upgrade the space A{C} to A{D}
convert_space_via_reduction(x, A, B) = _convert_space_via_reduction(x, A, B, embedding_reduction(A, B)...)
_convert_space_via_reduction(x, A, B, C, D) = convert_space(B, convert_space(similar_space(A, eltype(D)), x))

# Embedding via superspace: we convert x first to the superspace of A and then continue
convert_space_via_superspace(x, A, B) = convert_space(B, convert_space(superspace(A), x))


## And now for the converse: restrict_space

"""
Restrict the variable `x` to an element of the space `B`. This is possible if
`B` is embedded in the space of `x`.

Mathematically, `restrict_space(A, y)` for `y` in space `B` is a left inverse
of `convert_space(B, x)` for `x` in space `A`. This means that
`restrict_space(A, convert_space(B, x)) == x` for any `x` in space `A`. However,
`restrict_space(A, y)` for a `y` not in the range of `convert_space(B, x)` could
take any value.
"""# We don't need to do anything if the space of x is B
restrict_space(B::Type{GSpace{T}}, x::T) where {T} = x

# If it isn't, then we dispatch on the type of embedding
restrict_space(B::Type{GSpace{T}}, x) where {T} = restrict_spaces1(x, spaceof(x), B)

# If A and B are isomorphic, then we can use convert_space rather than restrict_space.
restrict_spaces1(x, A, B) = _restrict_spaces1(x, A, B, isomorphism(A, B))
_restrict_spaces1(x, A, B, ::Type{True}) = convert_space(B, x)
_restrict_spaces1(x, A, B, ::Type{False}) = restrict_spaces2(x, A, B)

# We have to follow the inverse logic compared to the case of conversions.
# Here, B is embedded in A and we figure out by which path.
restrict_spaces2(x, A, B) = restrict_spaces2(x, A, B,
    embedding_via_promotion(B, A),
    embedding_via_reduction(B, A),
    embedding_via_superspace(B, A))

restrict_spaces2(x, A, B, ::Type{True}, d2, d3) =
    restrict_space_via_promotion(x, A, B)
restrict_spaces2(x, A, B, ::Type{False}, ::Type{True}, d3) =
    restrict_space_via_reduction(x, A, B)
restrict_spaces2(x, A, B, ::Type{False}, ::Type{False}, ::Type{True}) =
    restrict_space_via_superspace(x, A, B)
restrict_spaces2(x, A, B, ::Type{False}, ::Type{False}, ::Type{False}) =
    throw(InexactError(:restrict_spaces, B, x))


# Embedding via promotion
restrict_space_via_promotion(x, A, B) = demote(eltype(B), x)

"""
The function `demote(S, y::T)` is a left inverse of `convert(T, x::S)`, for the
case where type `S` promotes to type `T` (i.e. `promote_type(S,T) == T`).

This means that `demote(S, convert(T, x::S)) == x`, while
`demote(S, y)` may be anything for `y::T` not in the range of
`convert(T, x::S)`.

Note that Julia's `convert` function is its own left inverse, in the sense that
`convert(S, convert(T,x)) == x` usually holds if `x::S`. However, the `convert`
function generally throws an `InexactError()` for elements of type `T` that are
not in (or close to) the range of `convert(T, x::S)`. In those cases, `demote`
does not throw an error, but there may be an arbitrarily large difference between
`convert(T, demote(S, y))` and `y` itself. In other words, `demote` may differ
wildly from a right inverse of `convert`.
"""
demote(::Type{T}, x) where {T} = convert(T, x)

# We need this rule, because Julia only converts a complex number to a real
# number if the imaginary part is less than some threshold.
demote(::Type{T}, x::Complex{T}) where {T} = real(x)

# Embedding via reduction: if B is embedded in A if C is embedded in D,
# then we have to downgrade the space A{D} to A{C}
restrict_space_via_reduction(x, A, B) = _restrict_space_via_reduction(x, A, B, embedding_reduction(B, A)...)
_restrict_space_via_reduction(x, A, B, C, D) = restrict_space(B, restrict_space(similar_space(A, eltype(C)), x))

# Embedding via superspace: we convert x first to the superspace of A and then continue
restrict_space_via_superspace(x, A, B) = restrict_space(B, restrict_space(superspace(B), x))


#############
# Promotions
#############

# We mimick Julia's promotion system, with the key difference that the types
# a space can be promoted to are not actual supertypes, but superspaces.

# Implementation:
# - the user calls promote_space(a,b)
# - the space to promote to is decided by promote_space_type(spaceof(a),spaceof(b))
# - promote_space_type is implemented in terms of promote_space_rule. This function
#   checks whether one of the spaces is embedded into the other.
# - If that is not the case, then both spaces (one at a time) are upgraded to
#   their superspaces and promote_space_type is called again.
# - It ends as soon as a joined promotable type is found. If that is not the case,
#   all branches end in AnySpace. In that case, a and b are interpreted as
#   elements of AnySpace and they remain unaltered.

"""
Promote the geometric spaces of the arguments to a joined space, to which all
arguments can be converted using embeddings. If no such concrete space exists,
the joined space is AnySpace and the arguments remain unaltered.
"""
promote_space() = ()
promote_space(x) = (x,)
promote_space(x::T, y::T) where {T} = (x,y)
function promote_space(x::T, y::S) where {T,S}
    space_type = promote_space_type(spaceof(x),spaceof(y))
    if space_type == AnySpace
      error("No conversion to AnySpace possible")
    end
    (convert_space(space_type, x), convert_space(space_type, y))
end

promote_space_type(::Type{AnySpace}, ::Type{AnySpace}) = AnySpace
promote_space_type(::Type{AnySpace}, ::Type{GSpace{T}}) where {T} = AnySpace
promote_space_type(::Type{GSpace{T}}, ::Type{AnySpace}) where {T} = AnySpace
promote_space_type(::Type{GSpace{T}}, ::Type{GSpace{T}}) where {T} = GSpace{T}

promote_space_type(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} =
    lcd(
        promote_via_promotion(A,B),
        promote_via_embedding_reduction(A,B),
        promote_via_isomorphism_reduction(A,B),
        promote_via_superspace(A,B))

promote_via_promotion(A::Type{GSpace{T}}, B::Type{GSpace{S}}) where {T,S} =
    _promote_via_promotion(A, B, T, S, promote_type(T,S))

# S promotes to T: return A
_promote_via_promotion(A, B, ::Type{T}, ::Type{S}, ::Type{T}) where {T,S} = A
# T promotes to S: return B
_promote_via_promotion(A, B, ::Type{T}, ::Type{S}, ::Type{S}) where {T,S} = B
# T and S are equal: this should have been caught before!
_promote_via_promotion(A, B, ::Type{T}, ::Type{T}, ::Type{T}) where {T} = A
# neither is the case: return AnySpace
_promote_via_promotion(A, B, ::Type{T}, ::Type{S}, ::Type{U}) where {T,S,U} = AnySpace


promote_via_superspace(A, B) = lcd(promote_space_type(superspace(A),B), promote_space_type(A,superspace(B)))

promote_via_embedding_reduction(A, B) = _promote_via_embedding_reduction(A, B, embedding_reduction(A,B)..., 1, embedding_reduction(B,A)...)
# We passed a dummy Int to distinguish between the cases below
_promote_via_embedding_reduction(A, B, ::Int) = AnySpace
# This line: A is embedded in B via reduction if C is in D: we try to promote C and D
_promote_via_embedding_reduction(A, B, C, D, ::Int) = _promote_via_embedding_reduction1(A, B, C, D, promote_space_type(C,D))
# - Not succesful: we bail out
_promote_via_embedding_reduction1(A, B, C, D, E::Type{AnySpace}) = AnySpace
# - Successful: we promote B{D} to B{E}
_promote_via_embedding_reduction1(A, B, C, D, E::Type{GSpace{T}}) where {T} = similar_space(B, T)
# This case is similar, but with B embedded in A
_promote_via_embedding_reduction(A, B, ::Int, C, D) = _promote_via_embedding_reduction2(A, B, C, D, promote_space_type(C,D))
_promote_via_embedding_reduction2(A, B, C, D, E::Type{AnySpace}) = AnySpace
# - in case promotion of C and D is succesful, we promote A{C} to A{E}
_promote_via_embedding_reduction2(A, B, C, D, E::Type{GSpace{T}}) where {T} = similar_space(A, T)
# Third case: A is embedded in B and B is in A: we rely on isomorphism_reduction
# in order to get the right preferential promotion type.
# The only risk we take here is that A is embedded in B and B is embedded in A, but
# for some reason the isomorphism between A and B is not defined.
_promote_via_embedding_reduction(A, B, C, D, ::Int, E, F) = AnySpace


promote_via_isomorphism_reduction(A, B) = _promote_via_isomorphism_reduction(A, B, isomorphism_reduction(A,B)..., 1, isomorphism_reduction(B,A)...)
_promote_via_isomorphism_reduction(A, B, ::Int) = AnySpace
# Here, we can choose either to promote A or B. We choose A because it was listed first in the reduction rule
_promote_via_isomorphism_reduction(A, B, C, D, ::Int, args...) = similar_space(A, eltype(promote_space_type(C,D)))
# In this second case B was listed first, hence we continue promotion with B.
_promote_via_isomorphism_reduction(A, B, ::Int, C, D) = similar_space(B, eltype(promote_space_type(C,D)))
