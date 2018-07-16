# common.jl


##############################
# Booleans in the type domain
##############################

# We introduce these types to compute with booleans in a type context.
# They are not exported, but they are necessary when an external package wishes
# to extend the embeddings and promotions of spaces in this package.
const True = Val{true}
const False = Val{false}

# Simple boolean operations on the new types
(|)(::Type{True}, ::Type{True}) = True
(|)(::Type{True}, ::Type{False}) = True
(|)(::Type{False}, ::Type{True}) = True
(|)(::Type{False}, ::Type{False}) = False

# Return True if one of the arguments is True
one_of(::Type{True}) = True
one_of(::Type{False}) = False
one_of(a::Type{Val{A}}, b::Type{Val{B}}) where {A,B} = |(one_of(a),one_of(b))
one_of(a::Type{Val{A}}, b::Type{Val{B}}, c::Type{Val{C}}, d...) where {A,B,C} = one_of(a, one_of(b, c, d...))

# Convert the boolean type to a boolean value
result(::Type{True}) = true
result(::Type{False}) = false



##############################
# Promotion helper functions
##############################

"Return True if S promotes to T, i.e., if promote_type(S,T) == T."
promotes_to(S, T) = _promotes_to(S, T, promote_type(S,T))
_promotes_to(::Type{S}, ::Type{T}, ::Type{T}) where {S,T} = True
_promotes_to(::Type{S}, ::Type{T}, ::Type{U}) where {S,T,U} = False



#######################
# Composite structures
#######################

"""
Some types have composite structure, e.g. product domains, a union of domains.
These types contain a list of domains.

It is often undesirable to use `getindex` to access the elements of the composite
type. For this reason we introduce the `elements` functions. Composite types
can implement `elements` and provide a generic way to access their components.

`elements(t)`: returns the elements making up the composite type `t`

`element(t, i)`: return the `i`-th element of the composite type `t`

`numelements(t)`: return the number of elements of the composite type `t`
"""
elements() = nothing

"""
Return the i-th element of a composite structure.

See also: `elements`.
"""
element(t, i) = elements(t)[i]
# By default, we index elements(o)

"""
Return the number of elements of a composite structure.

See also: `elements`.
"""# By default, we return length(elements(t))
numelements(t) = length(elements(t))


###############
# Type factory
###############

"""
A `TypeFactory{T}` is a convenience type to simplify construction of a type.

Having `t = TypeFactory{T}` overrides `getindex` such that `t[a]` invokes `T(a)`.

For example:
```
v = TypeFactory{SVector}
v[0.1,0.2]
```
makes an `SVector{2,Float64}`.
"""
struct TypeFactory{T}
end

Base.getindex(v::TypeFactory{T}, args...) where {T} = T(args...)

const v = TypeFactory{SVector}()

###############
# Type conversion
###############
Base.convert(::Type{SVector}, ::Type{NTuple{N,T}}) where {N,T} = SVector{N,T}
