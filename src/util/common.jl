
##############################
# Promotion helper functions
##############################

const True = Val{true}
const False = Val{false}

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
"""
numelements(t) = length(elements(t))

"Expand all arguments of type C into their components."
expand(::Type{C}) where {C} = ()
expand(::Type{C}, domain, domains...) where {C} =
    (domain, expand(C, domains...)...)
expand(::Type{C}, domain::C, domains...) where {C} =
    (elements(domain)..., expand(C, domains...)...)

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



#################
# Precision type
#################

"The floating point precision type associated with the argument."
prectype(x) = prectype(typeof(x))
prectype(::Type{<:Complex{T}}) where {T} = prectype(T)
prectype(::Type{<:AbstractArray{T}}) where {T} = prectype(T)
prectype(::Type{NTuple{N,T}}) where {N,T} = prectype(T)
prectype(::Type{Tuple{A}}) where {A} = prectype(A)
prectype(::Type{Tuple{A,B}}) where {A,B} = prectype(A,B)
prectype(::Type{Tuple{A,B,C}}) where {A,B,C} = prectype(A,B,C)
@generated function prectype(T::Type{<:Tuple{Vararg}})
    quote $(promote_type(map(prectype, T.parameters[1].parameters)...)) end
end
prectype(::Type{T}) where {T<:AbstractFloat} = T
prectype(::Type{T}) where {T} = prectype(float(T))

prectype(a, b) = promote_type(prectype(a), prectype(b))
prectype(a, b, c...) = prectype(prectype(a, b), c...)

#################
# Numeric type
#################

"The numeric element type of x in a Euclidean space."
numtype(x) = numtype(typeof(x))
numtype(::Type{T}) where {T<:Number} = T
numtype(::Type{T}) where {T} = eltype(T)
numtype(::Type{NTuple{N,T}}) where {N,T} = T
numtype(::Type{Tuple{A,B}}) where {A,B} = promote_type(numtype(A), numtype(B))
numtype(::Type{Tuple{A,B,C}}) where {A,B,C} = promote_type(numtype(A), numtype(B), numtype(C))
numtype(::Type{Tuple{A,B,C,D}}) where {A,B,C,D} = promote_type(numtype(A), numtype(B), numtype(C), numtype(D))
@generated function numtype(T::Type{<:Tuple{Vararg}})
    quote $(promote_type(T.parameters[1].parameters...)) end
end

numtype(a...) = promote_type(map(numtype, a)...)
