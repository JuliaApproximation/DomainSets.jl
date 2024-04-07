"""
    canonicalmap([ctype::CanonicalType, ]map)

Return an associated canonical map, if any, of the given map.

Optionally, a canonical type argument may specify an alternative canonical map.
Canonical maps help with converting between equal maps of different types.
"""
canonicalmap(m) = m

"Does the map have a canonical map?"
hascanonicalmap(m) = !(m === canonicalmap(m))


"Supertype of different kinds of canonical objects."
abstract type CanonicalType end

canonicalmap(ctype::CanonicalType, m) = m
hascanonicalmap(ctype::CanonicalType, m) = !(m === canonicalmap(ctype, m))


"""
    Equal <: CanonicalType

A canonical object that is equal but simpler.
"""
struct Equal <: CanonicalType end

canonicalmap(::Equal, m) = m
equalmap(m) = canonicalmap(Equal(), m)
hasequalmap(m) = hascanonicalmap(Equal(), m)

"""
    simplify(m)

Simplify the given map to an equal map.
"""
simplify(m) = equalmap(m)
"""
    simplifies(m)

Does the map simplify?
"""
simplifies(m) = hasequalmap(m)

"Convert the given map to a map defined in FunctionMaps.jl."
tofunctionmap(m::Map) = m
tofunctionmap(m::MapRef) = tofunctionmap(functionmap(m))


"""
    Equivalent <: CanonicalType

A canonical object that is equivalent but may have different type.
"""
struct Equivalent <: CanonicalType end

canonicalmap(::Equivalent, m) = canonicalmap(Equal(), m)

canonicalmap(::Equivalent, m::Map{<:StaticVector{1,T}}) where {T} = convert(Map{T}, m)
canonicalmap(::Equivalent, m::Map{NTuple{N,T}}) where {N,T} = convert(Map{SVector{N,T}}, m)

equivalentmap(m) = canonicalmap(Equivalent(), m)
hasequivalentmap(m) = hascanonicalmap(Equivalent(), m)

"""
    CanonicalExtensionType <: CanonicalType

Canonical types used to translate between packages.
"""
abstract type CanonicalExtensionType <: CanonicalType
end

"Return the extension type associated with the given object."
canonicalextensiontype(m) = canonicalextensiontype(typeof(m))


==(m1::Map, m2::Map) = isequalmap(m1, m2)   # Method from Base

"""
    isequalmap(map1, map2)

Are the two given maps equal?
"""
isequalmap(m1, m2) = isequalmap1(m1, m2)
isequalmap1(m1, m2) = simplifies(m1) ? isequalmap(simplify(m1), m2) : isequalmap2(m1, m2)
isequalmap2(m1, m2) = simplifies(m2) ? isequalmap(m1, simplify(m2)) : default_isequalmap(m1, m2)
default_isequalmap(m1, m2) = m1===m2

# Associated with == we have to define the hashes of maps
Base.hash(m::Map, h::UInt) = map_hash(m, h)

map_hash(m) = map_hash(m, zero(UInt))
map_hash(m, h::UInt) = invoke(hash, Tuple{Any,UInt}, m, h)
