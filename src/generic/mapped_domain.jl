
"""
A `MappedDomain` represents the mapping of a domain.

The map of a domain `d` under the mapping `y=f(x)` consists of all points `f(x)`
with `x ∈ d`. The characteristic function of a mapped domain is defined in
terms of the inverse map `g = inv(f)`:
```
x ∈ m ⟺ g(x) ∈ d
```
"""
abstract type AbstractMappedDomain{T} <: SingleLazyDomain{T} end

superdomain(d::AbstractMappedDomain) = d.domain

show(io::IO, d::AbstractMappedDomain) =  print(io, "A mapped domain based on ", superdomain(d))

tointernalpoint(d::AbstractMappedDomain, x) = inverse_map(d)(x)
toexternalpoint(d::AbstractMappedDomain, y) = forward_map(d)(y)

const MappedVectorDomain{T} = AbstractMappedDomain{Vector{T}}

# TODO: check whether the map alters the dimension
dimension(d::MappedVectorDomain) = dimension(superdomain(d))

# TODO: check whether the map affects these properties
isempty(d::AbstractMappedDomain) = isempty(superdomain(d))
isopenset(d::AbstractMappedDomain) = isopenset(superdomain(d))
isclosedset(d::AbstractMappedDomain) = isclosedset(superdomain(d))


"A `MappedDomain` stores the inverse map of a mapped domain."
struct MappedDomain{T,D,F} <: AbstractMappedDomain{T}
    domain  ::  D
    invmap  ::  F
end

MappedDomain(domain::Domain{T}, invmap) where {T} = MappedDomain{T}(domain, invmap)
MappedDomain(domain::Domain{S}, invmap::Map{T}) where {S,T} = MappedDomain{T}(domain, invmap)

MappedDomain{T}(domain::Domain, invmap) where {T} =
    MappedDomain{T,typeof(domain),typeof(invmap)}(domain, invmap)

# If the map is a Map{T}, make sure it matches the T of the MappedDomain
MappedDomain{T}(domain::Domain, invmap::Map{T}) where {T} =
    MappedDomain{T,typeof(domain),typeof(invmap)}(domain, invmap)
MappedDomain{T}(domain::Domain, invmap::Map{S}) where {S,T} =
    MappedDomain{T}(domain, convert(Map{T}, invmap))

convert(::Type{Domain{T}}, d::MappedDomain) where {T} = MappedDomain{T}(d.domain, d.invmap)

forward_map(d::MappedDomain) = inv(d.invmap)
forward_map(d::MappedDomain, x) = inv(d.invmap)(x)

inverse_map(d::MappedDomain) = d.invmap
inverse_map(d::MappedDomain, y) = d.invmap(y)

"Map a domain with the inverse of the given map"
map_domain(map::AbstractMap, domain::Domain) = mapped_domain(inv(map), domain)

"Make a mapped domain with the given inverse map"
mapped_domain(invmap, domain::Domain) = MappedDomain(domain, invmap)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
mapped_domain(invmap, d::MappedDomain) = mapped_domain(inverse_map(d) ∘ invmap, superdomain(d))
