
abstract type LazyMap{T} <: Map{T} end

abstract type CompositeLazyMap{T} <: LazyMap{T} end
abstract type SingleLazyMap{T} <: LazyMap{T} end

supermap(m::SingleLazyMap) = m.map
elements(m::CompositeLazyMap) = m.maps

isreal(m::SingleLazyMap) = isreal(supermap(m))
isreal(m::CompositeLazyMap) = all(map(isreal, elements(m)))

abstract type DerivedMap{T} <: SingleLazyMap{T} end

applymap(m::DerivedMap, x) = supermap(m)(x)
appymap!(y, m::DerivedMap, x) = applymap!(y, supermap(m), x)

jacobian(m::DerivedMap) = jacobian(supermap(m))


"A `WrappedMap{T}` takes any object and turns it into a `Map{T}`."
struct WrappedMap{T,M} <: DerivedMap{T}
    map ::  M
end
WrappedMap{T}(map) where {T} = WrappedMap{T,typeof(map)}(map)
WrappedMap(map) = WrappedMap{Float64}(map)

convert(::Type{Map{T}}, m::WrappedMap{T}) where {T} = m
convert(::Type{Map{T}}, m::WrappedMap) where {T} = WrappedMap{T}(m.map)

convert(::Type{Map}, m::Map) = m
convert(::Type{Map}, m) = WrappedMap(m)
convert(::Type{Map{T}}, m) where {T} = WrappedMap{T}(m)


"""
A lazy Jacobian `l` stores a map and implements `jacobian(l)` by invoking `jacobian(m, x)`
on the stored map.
"""
struct LazyJacobian{T,M} <: Map{T}
	map	::	M
end

LazyJacobian(m::Map{T}) where {T} = LazyJacobian{T,typeof(m)}(m)

applymap(m::LazyJacobian, x) = jacobian(m.map, x)
