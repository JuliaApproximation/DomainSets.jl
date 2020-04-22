
"An `AbstractMap` represents a function `y=f(x)` of a single variable."
abstract type AbstractMap end

"A `Map{T}` is a map of a single variable of type `T`."
abstract type Map{T} <: AbstractMap end

"A `TypedMap{T,U}` maps a variable of type `T` to a variable of type `U`."
abstract type TypedMap{T,U} <: Map{T} end

domaintype(m::AbstractMap) = domaintype(typeof(m))
domaintype(::Type{<:AbstractMap}) = Any
domaintype(::Type{<:Map{T}}) where {T} = T

codomaintype(m::AbstractMap) = codomaintype(typeof(m))
codomaintype(::Type{<:AbstractMap}) = Any
codomaintype(M::Type{<:Map{T}}) where {T} = Base.promote_op(applymap, M, T)
codomaintype(::Type{<:TypedMap{T,U}}) where {T,U} = U

convert(::Type{AbstractMap}, m::AbstractMap) = m
convert(::Type{Map{T}}, m::Map{T}) where {T} = m
convert(::Type{TypedMap{T,U}}, m::TypedMap{T,U}) where {T,U} = m

(m::AbstractMap)(x) = applymap(m, x)

# For maps of type Map{T}, we convert the point x or the map or both
(m::Map)(x) = _applymap(m, x)
_applymap(m::Map{T}, x::T) where {T} = applymap(m, x)
_applymap(m::Map{S}, x::T) where {S,T} = _applymap(m, x, promote_type(S,T))
_applymap(m::Map{S}, x::T, ::Type{V}) where {S,T,V} =
    applymap(convert(Map{V}, m), convert(V, x))

applymap!(y, m::AbstractMap, x) = y .= m(x)

# Note that there is a difference between a map being invertible, and the map
# knowing its inverse explicitly and implementing `inv`.
inv(m::AbstractMap) = error("Map ", m, " does not have a known inverse.")

(\)(map::AbstractMap, y) = inv(map)(y)

"""
Return a left inverse of the given map. This left inverse `mli` is not unique,
but in any case it is such that `(mli ∘ m) * x = x` for each `x` in the domain
of `m`.
"""
leftinv(m::AbstractMap) = inv(m)

"""
Return a right inverse of the given map. This right inverse `mri` is not unique,
but in any case it is such that `(m ∘ mri) * y = y` for each `y` in the range
of `m`.
"""
rightinv(m::AbstractMap) = inv(m)

"""
    jacobian(m::AbstractMap[, x])

Return the jacobian map. The two-argument version evaluates the jacobian
at a point `x`.
"""
jacobian(m::AbstractMap) = error("Map ", m, " does not have a known Jacobian.")
jacobian(m::AbstractMap, x) = jacobian(m)(x)

jacobian!(y, m::AbstractMap, x) = y .= jacobian(m, x)

"""
    jacdet(m::Map, x)

The Jacobian determinant of the map at a point `x`.
"""
jacdet(m::AbstractMap, x) = det(jacobian(m, x))



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

# size is not defined for all maps, because not all maps have a fixed size
# it is defined for maps that act on vectors, where the size information is not
# available from the type
size(m::AbstractMap, i) = i == i <= 2 ? size(m)[i] : 1