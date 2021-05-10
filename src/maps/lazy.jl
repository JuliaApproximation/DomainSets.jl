
"""
A lazy map has an action that is defined in terms of other maps. Those maps are
stored internally, and the action of the lazy map is computed on-the-fly and only
when invoked.
"""
abstract type LazyMap{T} <: Map{T} end

"A composite lazy map is defined in terms of several other maps."
abstract type CompositeLazyMap{T} <: LazyMap{T} end
"A simple lazy map derives from a single other map."
abstract type SimpleLazyMap{T} <: LazyMap{T} end

supermap(m::SimpleLazyMap) = m.map
components(m::CompositeLazyMap) = m.maps

isreal(m::SimpleLazyMap) = isreal(supermap(m))
isreal(m::CompositeLazyMap) = all(map(isreal, components(m)))

"A `DerivedMap` inherits all of its properties from another map, but has its own type."
abstract type DerivedMap{T} <: SimpleLazyMap{T} end

applymap(m::DerivedMap, x) = supermap(m)(x)
appymap!(y, m::DerivedMap, x) = applymap!(y, supermap(m), x)

jacobian(m::DerivedMap) = jacobian(supermap(m))


"A `WrappedMap{T}` takes any object and turns it into a `Map{T}`."
struct WrappedMap{T,M} <: DerivedMap{T}
    map ::  M
end
WrappedMap{T}(map) where {T} = WrappedMap{T,typeof(map)}(map)
WrappedMap(map) = WrappedMap{Float64}(map)

similarmap(m::WrappedMap, ::Type{T}) where {T} = WrappedMap{T}(m)

convert(::Type{Map}, m::Map) = m
convert(::Type{Map}, m) = WrappedMap(m)
convert(::Type{Map{T}}, m) where {T} = WrappedMap{T}(m)

==(m1::WrappedMap, m2::Function) = m1.map == m2
==(m1::Function, m2::WrappedMap) = m1 == m2.map

Display.displaystencil(m::WrappedMap{T}) where {T} =
	["WrappedMap{$T}(", supermap(m), ")"]
show(io::IO, mime::MIME"text/plain", m::WrappedMap) = composite_show(io, mime, m)
