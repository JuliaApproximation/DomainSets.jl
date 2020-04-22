
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
jacobian(m::DerivedMap, x) = jacobian(supermap(m), x)
jacobian!(y, m::DerivedMap, x) = jacobian!(y, supermap(m), x)
jacdet(m::DerivedMap, x) = jacdet(supermap(m), x)

struct WrappedMap{T,M} <: DerivedMap{T}
    map ::  M
end
WrappedMap{T}(map) where {T} = WrappedMap{T,typeof(map)}(map)
WrappedMap(map) = WrappedMap{Float64}(map)

convert(::Type{Map{T}}, m::WrappedMap{T}) where {T} = m
convert(::Type{Map{T}}, m::WrappedMap) where {T} = WrappedMap{T}(m.map)

convert(::Type{Map}, m::Map) = m
convert(::Type{Map}, m) = WrappedMap(m)