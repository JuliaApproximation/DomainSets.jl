
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

tointernalpoint(d::AbstractMappedDomain, x) = inverse_map(d, x)
toexternalpoint(d::AbstractMappedDomain, y) = forward_map(d, y)


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

# In the constructor, we have to decide which T to use for the MappedDomain.
# - we don't know anything about invmap: deduce T from the given domain
MappedDomain(domain::Domain{T}, invmap) where {T} = MappedDomain{T}(domain, invmap)
# - if the map is a Map{T}, use that T for the MappedDomain
MappedDomain(domain::Domain{S}, invmap::Map{T}) where {S,T} = MappedDomain{T}(domain, invmap)
# - if T is given in the constructor, by all means we use that
MappedDomain{T}(domain::Domain, invmap) where {T} =
    MappedDomain{T,typeof(domain),typeof(invmap)}(domain, invmap)
# - in that case, if the map is a Map{S}, make sure that S matches T
MappedDomain{T}(domain::Domain, invmap::Map{T}) where {T} =
    MappedDomain{T,typeof(domain),typeof(invmap)}(domain, invmap)
MappedDomain{T}(domain::Domain, invmap::Map{S}) where {S,T} =
    MappedDomain{T}(domain, convert(Map{T}, invmap))

similardomain(d::MappedDomain, ::Type{T}) where {T} =
    MappedDomain{T}(d.domain, d.invmap)

forward_map(d::MappedDomain) = inv(d.invmap)
forward_map(d::MappedDomain, x) = inverse(d.invmap, x)

inverse_map(d::MappedDomain) = d.invmap
inverse_map(d::MappedDomain, y) = d.invmap(y)

"Map a domain with the inverse of the given map"
map_domain(map, domain::Domain) = _map_domain(map, domain)

# Fallback: we don't know anything about map, just try to invert
_map_domain(map, domain) = mapped_domain(inv(map), domain)
# If map is a Map{T}, then verify and if necessary update T
function _map_domain(map::Map, domain)
    U = codomaintype(map, eltype(domain))
    mapped_domain(inv(convert(Map{U}, map)), domain)
end

==(a::MappedDomain, b::MappedDomain) = (a.invmap == b.invmap) && (superdomain(a) == superdomain(b))


"Make a mapped domain with the given inverse map"
mapped_domain(invmap, domain::Domain) = _mapped_domain(invmap, domain)

# We face the same task as in the constructor: attempt to identify T
# Here, we are more flexible, and attempt to do more conversions. We assume
# that users invoking the MappedDomain constructor know what they are doing,
# but users invoking mapped_domain just expect it to work.

# - we don't know anything about invmap, just pass it on
_mapped_domain(invmap, domain) = MappedDomain(domain, invmap)
# - invmap is a Map{T}: its codomaintype should match the eltype of the domain
# -- first, update the numtype
_mapped_domain(invmap::Map, domain) =
    _mapped_domain(invmap, domain, promote_type(numtype(invmap),numtype(domain)))
_mapped_domain(invmap::Map{T}, domain::Domain{S}, ::Type{U}) where {S,T,U} =
    _mapped_domain2(convert_numtype(invmap,U), convert_numtype(domain,U))
# -- then, ensure the codomaintype of the map equals the element type of the domain
_mapped_domain2(invmap, domain) = _mapped_domain2(invmap, domain, codomaintype(invmap), eltype(domain))
# --- it's okay
_mapped_domain2(invmap, domain, ::Type{T}, ::Type{T}) where {T} =
    MappedDomain(domain, invmap)
# --- it's not okay: attempt to convert the map (triggers e.g. when combining a scalar
#        LinearMap with a vector domain)
_mapped_domain2(invmap, domain, ::Type{S}, ::Type{T}) where {S,T} =
    MappedDomain(domain, convert(Map{T}, invmap))


(∘)(domain::Domain, invmap::Union{Function,AbstractMap}) = mapped_domain(invmap, domain)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
mapped_domain(invmap, d::MappedDomain) = mapped_domain(inverse_map(d) ∘ invmap, superdomain(d))

canonicaldomain(d::MappedDomain) = superdomain(d)
fromcanonical(d::MappedDomain) = forward_map(d)
tocanonical(d::MappedDomain) = inverse_map(d)

# bijection1(d1::MappedDomain, d2) =
#     bijection(superdomain(d1), d2) ∘ inverse_map(d1)
# bijection2(d1, d2::MappedDomain) =
#     forward_map(d2) ∘ bijection(d1, superdomain(d2))
