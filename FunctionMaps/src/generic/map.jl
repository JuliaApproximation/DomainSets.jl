
"An `AbstractMap` represents a function `y=f(x)` of a single variable."
abstract type AbstractMap end

"A `Map{T}` is a map of a single variable of type `T`."
abstract type Map{T} <: AbstractMap end

Map(m) = convert(Map, m)
Map{T}(m) where {T} = convert(Map{T}, m)

"A `TypedMap{T,U}` maps a variable of type `T` to a variable of type `U`."
abstract type TypedMap{T,U} <: Map{T} end

const EuclideanMap{N,T} = Map{<:StaticVector{N,T}}
const VectorMap{T} = Map{Vector{T}}

CompositeTypes.Display.displaysymbol(m::Map) = 'F'

"What is the expected type of a point in the domain of the function map `m`?"
domaintype(m) = domaintype(typeof(m))
domaintype(::Type{M}) where {M} = Any
domaintype(::Type{<:Map{T}}) where {T} = T

"""
    codomaintype(m[, T])

What is the codomain type of the function map `m`, given that `T` is its domain type?
"""
codomaintype(m) = codomaintype(m, domaintype(m))
codomaintype(m, ::Type{T}) where {T} = codomaintype(typeof(m), T)

codomaintype(::Type{M}, ::Type{T}) where {M,T} = Any
codomaintype(::Type{M}, ::Type{Any}) where {M} = Any
codomaintype(M::Type{<:AbstractMap}, ::Type{T}) where {T} = Base.promote_op(applymap, M, T)
codomaintype(M::Type{<:AbstractMap}, ::Type{Any}) = Any
codomaintype(M::Type{<:TypedMap{T,U}}, ::Type{T}) where {T,U} = U

prectype(::Type{<:Map{T}}) where T = prectype(T)
numtype(::Type{<:Map{T}}) where T = numtype(T)

isreal(m::AbstractMap) = isreal(domaintype(m)) && isreal(codomaintype(m))
isreal(::UniformScaling{T}) where {T} = isreal(T)
isreal(::Type{UniformScaling{T}}) where {T} = isreal(T)

convert(::Type{AbstractMap}, m::AbstractMap) = m
convert(::Type{Map{T}}, m::Map{T}) where {T} = m
convert(::Type{Map{T}}, m::Map{S}) where {S,T} = similarmap(m, T)
convert(::Type{TypedMap{T,U}}, m::TypedMap{T,U}) where {T,U} = m
convert(::Type{TypedMap{T,U}}, m::TypedMap) where {T,U} = similarmap(m, T, U)

convert_domaintype(::Type{T}, map::Map{T}) where {T} = map
convert_domaintype(::Type{U}, map::Map{T}) where {T,U} = convert(Map{U}, map)
convert_domaintype(::Type{Any}, map) = map
convert_domaintype(::Type{Any}, map::Map{T}) where T = map

convert_numtype(::Type{U}, map::Map{T}) where {T,U} = convert(Map{to_numtype(U,T)}, map)
convert_prectype(::Type{U}, map::Map{T}) where {T,U} = convert(Map{to_prectype(U,T)}, map)

convert_codomaintype(::Type{U}, map) where U =
    _convert_codomaintype(U, map, domaintype(map), codomaintype(map))
# types match: we can return map
_convert_codomaintype(::Type{U}, map, ::Type{T}, ::Type{U}) where {T,U} = map
# types don't match: we first convert the numtype
_convert_codomaintype(::Type{U}, map, ::Type{T}, ::Type{V}) where {T,U,V} =
    _convert_codomaintype2(U, convert_numtype(numtype(U), map))
_convert_codomaintype2(::Type{U}, map::Map) where U =
    _convert_codomaintype2(U, map, domaintype(map), codomaintype(map))
# types match this time around: we can return map
_convert_codomaintype2(::Type{U}, map, ::Type{T}, ::Type{U}) where {T,U} = map
# types don't match: we can't do this automatically
_convert_codomaintype2(::Type{U}, map, ::Type{T}, ::Type{V}) where {T,U,V} =
    throw(ArgumentError("Don't know how to convert the codomain type of $(map) to $(U)."))

convert_codomaintype(::Type{Any}, map) = map


# Users may call a map, concrete subtypes specialize the `applymap` function
(m::AbstractMap)(x) = applymap(m, x)

# For Map{T}, we allow invocation with multiple arguments by conversion to T
(m::Map{T})(x) where {T} = promote_and_apply(m, x)
(m::Map{T})(x...) where {T} = promote_and_apply(m, convert(T, x))

"Promote map and point to compatible types."
promote_map_point_pair(m, x) = (m, x)
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
promote_and_apply(m::Map, x) = applymap(promote_map_point_pair(m,x)...)

applymap!(y, m, x) = y .= m(x)

# Fallback for functions that are not of type AbstractMap
applymap(m, x) = m(x)
applymap(m::AbstractMap, x) = error("Please implement applymap for map $(m)")


"Can the maps be promoted to a common domain type without throwing an error?"
promotable_maps(maps...) = promotable_eltypes(map(domaintype, maps)...)

"Promote the given maps to have a common domain type."
promote_maps() = ()
promote_maps(m1) = m1
function promote_maps(m1, m2)
    T = promote_type(domaintype(m1), domaintype(m2))
    convert_domaintype(T, m1), convert_domaintype(T, m2)
end
function promote_maps(m1, m2, m3, maps...)
    T = mapreduce(domaintype, promote_type, (m1,m2,m3,maps...))
    convert_domaintype.(T, (m1,m2,m3,maps...))
end

isvectorvalued_type(::Type{T}) where {T<:Number} = true
isvectorvalued_type(::Type{T}) where {T<:AbstractVector} = true
isvectorvalued_type(::Type{T}) where {T} = false


"Is the map a vector-valued function, i.e., a function from Rn to Rm?"
isvectorvalued(m) =
    isvectorvalued_type(domaintype(m)) && isvectorvalued_type(codomaintype(m))

# mapsize should be defined for vector valued maps
# The size of a map equals the size of its jacobian
# The jacobian can be a number, a vector, an adjoint vector, or a matrix
mapsize(m, i) = _mapsize(m, i, mapsize(m))
_mapsize(m, i, size::Tuple{Int,Int}) = i <= 2 ? size[i] : 1
_mapsize(m, i, size::Tuple{Int}) = i <= 1 ? size[i] : 1
_mapsize(m, i, size::Tuple{}) = 1

"Is the given map a square map?"
issquaremap(m) = isvectorvalued(m) && (mapsize(m,1) == mapsize(m,2))

isoverdetermined(m) = mapsize(m,1) >= mapsize(m,2)
isunderdetermined(m) = mapsize(m,1) <= mapsize(m,2)

is_scalar_to_vector(m) = mapsize(m) isa Tuple{Int}
is_vector_to_scalar(m) = mapsize(m) isa Tuple{Int,Int} && codomaintype(m)<:Number
is_scalar_to_scalar(m) = mapsize(m) == ()
is_vector_to_vector(m) = mapsize(m) isa Tuple{Int,Int} && !is_vector_to_scalar(m)

# Display routines
map_stencil(m, x) = [Display.SymbolObject(m), '(', x, ')']
map_stencil_broadcast(m, x) = [Display.SymbolObject(m), ".(", x, ')']
