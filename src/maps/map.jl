
"An `AbstractMap` represents a function `y=f(x)` of a single variable."
abstract type AbstractMap end

"A `Map{T}` is a map of a single variable of type `T`."
abstract type Map{T} <: AbstractMap end

"A `TypedMap{T,U}` maps a variable of type `T` to a variable of type `U`."
abstract type TypedMap{T,U} <: Map{T} end

const EuclideanMap{N,T} = Map{<:StaticVector{N,T}}
const VectorMap{T} = Map{Vector{T}}

CompositeTypes.Display.displaysymbol(m::Map) = 'F'

domaintype(m::AbstractMap) = domaintype(typeof(m))
domaintype(::Type{<:AbstractMap}) = Any
domaintype(::Type{<:Map{T}}) where {T} = T

codomaintype(m::AbstractMap) = codomaintype(typeof(m))
codomaintype(::Type{<:AbstractMap}) = Any
codomaintype(M::Type{<:Map{T}}) where {T} = Base.promote_op(applymap, M, T)
codomaintype(::Type{<:TypedMap{T,U}}) where {T,U} = U

# What is the output type given an argument of type S?
codomaintype(m::AbstractMap, ::Type{S}) where {S} = codomaintype(typeof(m), S)
codomaintype(M::Type{<:AbstractMap}, ::Type{S}) where {S} = Base.promote_op(applymap, M, S)
codomaintype(::Type{<:TypedMap{T,U}}, ::Type{T}) where {T,U} = U

isreal(m::Map) = isreal(domaintype(m)) && isreal(codomaintype(m))
isreal(::UniformScaling{T}) where {T} = isreal(T)
isreal(::Type{UniformScaling{T}}) where {T} = isreal(T)

numtype(::Type{<:Map{T}}) where {T} = numtype(T)
prectype(::Type{<:Map{T}}) where {T} = prectype(T)

convert(::Type{AbstractMap}, m::AbstractMap) = m
convert(::Type{Map{T}}, m::Map{T}) where {T} = m
convert(::Type{Map{T}}, m::Map{S}) where {S,T} = similarmap(m, T)
convert(::Type{TypedMap{T,U}}, m::TypedMap{T,U}) where {T,U} = m
convert(::Type{TypedMap{T,U}}, m::TypedMap) where {T,U} = similarmap(m, T, U)

convert_numtype(map::Map{T}, ::Type{U}) where {T,U} = convert(Map{to_numtype(T,U)}, map)
convert_prectype(map::Map{T}, ::Type{U}) where {T,U} = convert(Map{to_prectype(T,U)}, map)

compatible_domaintype(d1, d2) = isconcretetype(promote_type(domaintype(d1),domaintype(d2)))

# Users may call a map, concrete subtypes specialize the `applymap` function
(m::AbstractMap)(x) = applymap(m, x)

# For Map{T}, we allow invocation with multiple arguments by conversion to T
(m::Map{T})(x) where {T} = apply(m, x)
(m::Map{T})(x...) where {T} = apply(m, convert(T, x))

"Promote map and point to compatible types."
promote_map_point_pair(m, x) = (m, x)
promote_map_point_pair(m::AbstractMap, x) = (m, x)
promote_map_point_pair(m::Map, x) = _promote_map_point_pair(m, x, promote_type(domaintype(m), typeof(x)))
# This is the line where we promote both the map and the point:
_promote_map_point_pair(m, x, ::Type{T}) where {T} =
    convert(Map{T}, m), convert(T, x)
_promote_map_point_pair(m, x, ::Type{Any}) = m, x

# Some exceptions:
# - types match, do nothing
promote_map_point_pair(m::Map{T}, x::T) where {T} = m, x
# - an Any map, do nothing
promote_map_point_pair(m::Map{Any}, x) = m, x
# - tuples: these are typically composite maps, promotion may happen later
promote_map_point_pair(m::Map{<:Tuple}, x::Tuple) = m, x
# - abstract vectors: promotion may be expensive
promote_map_point_pair(m::Map{<:AbstractVector}, x::AbstractVector) = m, x
# - SVector: promotion is likely cheap
promote_map_point_pair(m::EuclideanMap{N,T}, x::AbstractVector{S}) where {N,S,T} =
    _promote_map_point_pair(m, x, SVector{N,promote_type(S,T)})

# For maps of type Map{T}, we call promote_map_point_pair and then applymap
apply(m::Map, x) = applymap(promote_map_point_pair(m,x)...)
applymap!(y, m::AbstractMap, x) = y .= m(x)


isvectorvalued_type(::Type{T}) where {T<:Number} = true
isvectorvalued_type(::Type{T}) where {T<:AbstractVector} = true
isvectorvalued_type(::Type{T}) where {T} = false

"Is the map a vector-valued function, i.e., a function from Rn to Rm?"
isvectorvalued(m::Map{T}) where {T} =
    isvectorvalued_type(T) && isvectorvalued_type(codomaintype(m))

# size should be defined for vector valued maps
# The size of a map equals the size of its jacobian
# The jacobian can be a number, a vector, an adjoint vector, or a matrix
size(m::AbstractMap, i) = _map_size(m, i, size(m))
_map_size(m::AbstractMap, i, size::Tuple{Int,Int}) = i <= 2 ? size[i] : 1
_map_size(m::AbstractMap, i, size::Tuple{Int}) = i <= 1 ? size[i] : 1
_map_size(m::AbstractMap, i, size::Tuple{}) = 1

"Is the given map a square map?"
issquare(m::AbstractMap) = isvectorvalued(m) && (size(m,1) == size(m,2))

isoverdetermined(m::AbstractMap) = size(m,1) >= size(m,2)
isunderdetermined(m::AbstractMap) = size(m,1) <= size(m,2)

is_scalar_to_vector(m::Map) = size(m) isa Tuple{Int}
is_vector_to_scalar(m::Map) = size(m) isa Tuple{Int,Int} && codomaintype(m)<:Number
is_scalar_to_scalar(m::Map) = size(m) == ()
is_vector_to_vector(m::Map) = size(m) isa Tuple{Int,Int} && !is_vector_to_scalar(m)


# Display routines
map_stencil(m::AbstractMap, x) = [SymbolObject(m), '(', x, ')']
map_stencil_broadcast(m::AbstractMap, x) = [SymbolObject(m), ".(", x, ')']
