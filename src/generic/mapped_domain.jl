
"""
A `MappedDomain` represents the mapping of a domain.

The characteristic function of the mapped domain is defined in terms of the
inverse of the map `f`, i.e.:
```
in(x, mappeddomain) = in(inv(f)(x), domain)
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

MappedDomain(domain::Domain, invmap::Map{T}) where {T} = MappedDomain{T}(domain, invmap)

MappedDomain{T}(domain::Domain, invmap) where {T} =
    MappedDomain{T,typeof(domain),typeof(invmap)}(domain, invmap)

forward_map(d::MappedDomain) = inv(d.invmap)
forward_map(d::MappedDomain, x) = inv(d.invmap)(x)

inverse_map(d::MappedDomain) = d.invmap
inverse_map(d::MappedDomain, y) = d.invmap(y)

inversemap_domain(invmap, domain::Domain) = MappedDomain(domain, invmap)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
inversemap_domain(invmap, d::MappedDomain) = inversemap_domain(inverse_map(d) ∘ invmap, superdomain(d))
