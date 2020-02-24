
"A `Map{T}` represents a function `y=f(x)` of a single variable of type `T`."
abstract type Map{T} end

"Supertype of maps `y=f(x)` where `x` has type `T` and `y` has type `U`."
abstract type TypedMap{T,U} <: Map{T} end

domaintype(m::Map) = domaintype(typeof(m))
domaintype(::Type{<:Map{T}}) where {T} = T

codomaintype(map::TypedMap) = codomaintype(typeof(map))
codomaintype(::Type{<:TypedMap{T,U}}) where {T,U} = U

convert(::Type{Map{T}}, m::Map{T}) where {T} = m
convert(::Type{TypedMap{T,U}}, m::TypedMap{T,U}) where {T,U} = m

(m::Map)(x) = _applymap(m, x)
_applymap(m::Map{T}, x::T) where {T} = applymap(m, x)
_applymap(m::Map{S}, x::T) where {S,T} = _applymap(m, x, promote_type(S,T))
_applymap(m::Map{S}, x::T, ::Type{U}) where {S,T,U} =
    applymap(convert(Map{U}, m), convert(U, x))

call!(y, m::Map, x) = _call!(y, m, x)
_call!(y, m::Map{T}, x::T) where {T} = applymap!(y, m, x)
_call!(y, m::Map{S}, x::T) where {S,T} = _call!(y, m, x, promote_type(S,T))
_call!(y, m::Map{S}, x::T, ::Type{U}) where {S,T,U} =
    applymap!(y, convert(Map{U}, m), convert(U, x))

# Note that there is a difference between a map being invertible, and the map
# knowing its inverse explicitly and implementing `inv`.
inv(m::Map) = error("Map ", m, " does not have a known inverse.")

(\)(map::Map, y) = inv(map)(y)

"""
Return a left inverse of the given map. This left inverse `mli` is not unique,
but in any case it is such that `(mli ∘ m) * x = x` for each `x` in the domain
of `m`.
"""# By default we invoke inv
leftinv(m::Map) = inv(m)

"""
Return a right inverse of the given map. This right inverse `mri` is not unique,
but in any case it is such that `(m ∘ mri) * y = y` for each `y` in the range
of `m`.
"""# By default we invoke inv
rightinv(m::Map) = inv(m)
