
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

numtype(::Type{<:Map{T}}) where {T} = numtype(T)
prectype(::Type{<:Map{T}}) where {T} = prectype(T)

convert(::Type{AbstractMap}, m::AbstractMap) = m
convert(::Type{Map{T}}, m::Map{T}) where {T} = m
convert(::Type{Map{T}}, m::Map{S}) where {S,T} = similarmap(m, T)
convert(::Type{TypedMap{T,U}}, m::TypedMap{T,U}) where {T,U} = m
convert(::Type{TypedMap{T,U}}, m::TypedMap) where {T,U} = similarmap(m, T, U)

convert_numtype(map::Map{T}, ::Type{U}) where {T,U} = convert(Map{to_numtype(T,U)}, map)
convert_prectype(map::Map{T}, ::Type{U}) where {T,U} = convert(Map{to_prectype(T,U)}, map)

dimension(m::Map{T}) where {T} = euclideandimension(T)

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

# There is a difference between a map being invertible, and the map
# knowing its inverse explicitly and implementing `inv`.
inv(m::AbstractMap) = error("Map ", m, " does not have a known inverse map.")

# For inverse functions, we supply a two-argument version because finding the
# inverse for a given point is often cheaper than finding the inverse and then
# applying it to a point.

"""
    inverse(m::AbstractMap[, x])

Return the inverse of `m`. The two-argument function evaluates the inverse
at the point `x`.
"""
inverse(m::AbstractMap) = inv(m)
inverse(m::AbstractMap, x) = inverse(m)(x)

(\)(m::AbstractMap, y) = inverse(m, y)

"""
    leftinverse(m::AbstractMap[, x])

Return a left inverse of the given map. This left inverse `mli` is not unique,
but in any case it is such that `(mli ∘ m) * x = x` for each `x` in the domain
of `m`.

The two-argument function applies the left inverse to the point `x`.
"""
leftinverse(m::AbstractMap) = inverse(m)
leftinverse(m::AbstractMap, x) = leftinverse(m)(x)

"""
    rightinverse(m::AbstractMap[, x])

Return a right inverse of the given map. This right inverse `mri` is not unique,
but in any case it is such that `(m ∘ mri) * y = y` for each `y` in the range
of `m`.

The two-argument function applies the right inverse to the point `x`.
"""
rightinverse(m::AbstractMap) = inverse(m)
rightinverse(m::AbstractMap, x) = rightinverse(m)(x)


"""
    jacobian(m::AbstractMap[, x])

Return the jacobian map. The two-argument version evaluates the jacobian
at a point `x`.
"""
jacobian(m::AbstractMap) = error("Map ", m, " does not have a known Jacobian map.")
jacobian(m::AbstractMap, x) = jacobian(m)(x)

jacobian!(y, m::AbstractMap, x) = y .= jacobian(m, x)

"""
    jacdet(m::Map, x)

The Jacobian determinant of the map at a point `x`.
"""
jacdet(m::AbstractMap, x) = det(jacobian(m, x))

# Display routines
map_stencil(m::AbstractMap, x) = [SymbolObject(m), '(', x, ')']
map_stencil_broadcast(m::AbstractMap, x) = [SymbolObject(m), ".(", x, ')']

"Return the expected type of the jacobian matrix of the map"
matrixtype(m::AbstractMap) = matrixtype(domaintype(m), codomaintype(m))

# We check all combinations of scalars, static vector and abstract vectors
matrixtype(::Type{T}, ::Type{U}) where {T,U} = promote_type(T,U)
matrixtype(::Type{T}, ::Type{U}) where {T<:Number,U<:Number} = promote_type(T,U)
matrixtype(::Type{SVector{N,T}}, ::Type{SVector{M,U}}) where {N,M,T,U} = SMatrix{M,N,promote_type(T,U)}
matrixtype(::Type{SVector{N,T}}, ::Type{U}) where {N,T,U<:Number} = SMatrix{1,N,promote_type(T,U)}
matrixtype(::Type{T}, ::Type{SVector{M,U}}) where {T<:Number,M,U} = SVector{M,promote_type(T,U)}
matrixtype(::Type{<:AbstractVector{T}}, ::Type{<:AbstractVector{U}}) where {T,U} = Array{promote_type{T,U},2}
matrixtype(::Type{<:AbstractVector{T}}, ::Type{U}) where {T,U<:Number} = Array{promote_type{T,U},2}
matrixtype(::Type{T}, ::Type{<:AbstractVector{U}}) where {T<:Number,U} = Array{promote_type(T,U),1}

# size may not be defined for all maps, because not all maps have a fixed size
size(m::AbstractMap, i) = i <= 2 ? size(m)[i] : 1
