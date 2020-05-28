

isreal(::Type{<:Real}) = true
isreal(::Type{<:Complex}) = false
isreal(::Type{T}) where {T} = isreal(eltype(T))


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
elements(t) = ()

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

"""
Is `x` composed of different elements?

See also: `elements`.
"""
iscomposite(t) = length(elements(t)) > 0

"Expand all arguments of type C into their components."
expand(::Type{C}) where {C} = ()
expand(::Type{C}, domain, domains...) where {C} =
    (domain, expand(C, domains...)...)
expand(::Type{C}, domain::C, domains...) where {C} =
    (elements(domain)..., expand(C, domains...)...)



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
prectype(::Type{T}) where {T<:Number} = prectype(float(T))

prectype(a...) = promote_type(map(prectype, a)...)

"Convert `x` such that its `prectype` equals `U`."
convert_prectype(x, ::Type{U}) where {U} = convert(convert_prectype(typeof(x),U), x)

convert_prectype(::Type{T}, ::Type{U}) where {T,U} = error("Don't know how to convert the numtype of $(T) to $(U).")
convert_prectype(::Type{T}, ::Type{U}) where {T <: Real,U <: Real} = U
convert_prectype(::Type{Complex{T}}, ::Type{U}) where {T <: Real,U <: Real} = Complex{U}
convert_prectype(::Type{SVector{N,T}}, ::Type{U}) where {N,T,U} = SVector{N,convert_prectype(T,U)}
convert_prectype(::Type{Vector{T}}, ::Type{U}) where {T,U} = Vector{convert_prectype(T,U)}

"Promote the precision types of the arguments to a joined supertype."
promote_prectype(a) = a
promote_prectype(a, b) = _promote_prectype(prectype(a,b), a, b)
promote_prectype(a, b, c...) = _promote_prectype(prectype(a,b,c...), a, b, c...)
_promote_prectype(U, a) = convert_prectype(a, U)
_promote_prectype(U, a, b) = convert_prectype(a, U), convert_prectype(b, U)
_promote_prectype(U, a, b, c...) =
    (convert_prectype(a, U), convert_prectype(b, U), _promote_prectype(U, c...)...)

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

"Convert `x` such that its `numtype` equals `U`."
convert_numtype(x, ::Type{U}) where {U} = convert(convert_numtype(typeof(x),U), x)

convert_numtype(::Type{T}, ::Type{U}) where {T,U} = error("Don't know how to convert the numtype of $(T) to $(U).")
convert_numtype(::Type{T}, ::Type{U}) where {T <: Number,U <: Number} = U
convert_numtype(::Type{SVector{N,T}}, ::Type{U}) where {N,T,U} = SVector{N,convert_numtype(T,U)}
convert_numtype(::Type{Vector{T}}, ::Type{U}) where {T,U} = Vector{convert_numtype(T,U)}

"Promote the numeric types of the arguments to a joined supertype."
promote_numtype(a) = a
promote_numtype(a, b) = _promote_numtype(numtype(a,b), a, b)
promote_numtype(a, b, c...) = _promote_numtype(numtype(a,b,c...), a, b, c...)
_promote_numtype(U, a) = convert_numtype(a, U)
_promote_numtype(U, a, b) = convert_numtype(a, U), convert_numtype(b, U)
_promote_numtype(U, a, b, c...) =
    (convert_numtype(a, U), convert_numtype(b, U), _promote_numtype(U, c...)...)
