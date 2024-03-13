
"""
    MapStyle(m)

Trait to indicate whether or not `m` implements the map interface.
"""
abstract type MapStyle end

"""
    IsMap()

indicates an object implements the map interface.
"""
struct IsMap <: MapStyle end

"""
    NotMap()

indicates an object does not implement the Map interface.
"""
struct NotMap <: MapStyle end


MapStyle(m) = MapStyle(typeof(m))
# the default is no map
MapStyle(::Type) = NotMap()
# subtypes of Map are maps
MapStyle(::Type{<:Map}) = IsMap()

"""
    functionmap(m)

Return a map associated with the object `m`.
"""
functionmap(m::Map) = m

"""
    MapRef(m)

A reference to a map.

In a function call, `MapRef(x)` can be used to indicate that `x` should be
treated as a map, e.g., `foo(x, MapRef(m))`.
"""
struct MapRef{M}
    map  ::  M
end

functionmap(m::MapRef) = m.map
domaintype(m::MapRef) = domaintype(functionmap(m))


"""
`AnyMap` is the union of `Map` and `MapRef`.

In both cases `map(m::AnyMap)` returns the map itself.
"""
const AnyMap = Union{Map,MapRef}

"""
   checkmap(m)

Checks that `m` is a map or refers to a map and if so returns that map,
throws an error otherwise.
"""
checkmap(m::Map) = m
# we trust the explicit intention of a user providing a map reference
checkmap(m::MapRef) = functionmap(m)
# for other objects we check MapStyle
checkmap(m) = _checkmap(m, MapStyle(m))
_checkmap(m, ::IsMap) = m
_checkmap(m, ::NotMap) =
    throw(ArgumentError("Map does not implement map interface as indicated by MapStyle."))
