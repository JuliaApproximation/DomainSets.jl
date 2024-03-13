
"""
A `MappedDomain` represents the mapping of a domain.

The map of a domain `d` under the mapping `y=f(x)` consists of all points `f(x)`
with `x ∈ d`. The characteristic function of a mapped domain is defined in
terms of the inverse map `g = inverse(f)`:
```
x ∈ m ⟺ g(x) ∈ d
```
"""
abstract type AbstractMappedDomain{T} <: SimpleLazyDomain{T} end

superdomain(d::AbstractMappedDomain) = d.domain

const MappedVectorDomain{T} = AbstractMappedDomain{Vector{T}}

tointernalpoint(d::AbstractMappedDomain, x) = inverse_map(d, x)
toexternalpoint(d::AbstractMappedDomain, y) = forward_map(d, y)

canonicaldomain(d::AbstractMappedDomain) = canonicaldomain(superdomain(d))
mapfrom_canonical(d::AbstractMappedDomain) = forward_map(d) ∘ mapfrom_canonical(superdomain(d))
mapto_canonical(d::AbstractMappedDomain) = mapto_canonical(superdomain(d)) ∘ inverse_map(d)


# TODO: check whether the map alters the dimension
dimension(d::MappedVectorDomain) = dimension(superdomain(d))

# TODO: check whether the map affects these properties
isempty(d::AbstractMappedDomain) = isempty(superdomain(d))
isopenset(d::AbstractMappedDomain) = isopenset(superdomain(d))
isclosedset(d::AbstractMappedDomain) = isclosedset(superdomain(d))

corners(d::AbstractMappedDomain) = [forward_map(d, x) for x in corners(superdomain(d))]

## I/O functionality

show(io::IO, mime::MIME"text/plain", d::AbstractMappedDomain) = composite_show(io, mime, d)
Display.displaystencil(d::AbstractMappedDomain) =
    map_stencil_broadcast(forward_map(d), superdomain(d))
Display.object_parentheses(d::AbstractMappedDomain) =
    Display.object_parentheses(forward_map(d))
Display.stencil_parentheses(d::AbstractMappedDomain) =
    Display.stencil_parentheses(forward_map(d))



"A `MappedDomain` stores the inverse map of a mapped domain."
struct MappedDomain{T,F,D} <: AbstractMappedDomain{T}
    invmap  ::  F
    domain  ::  D
end

# In the constructor, we have to decide which T to use for the MappedDomain.
# - if we don't know anything about invmap: deduce T from the given domain
MappedDomain(invmap, domain) =
    MappedDomain{domaineltype(domain)}(invmap, domain)
# - if the map is a Map{T}, use that T for the MappedDomain
MappedDomain(invmap::Map{T}, domain) where {T} =
    MappedDomain{T}(invmap, domain)
# If T is given in the constructor, by all means we use that:
MappedDomain{T}(invmap, domain) where {T} =
    _MappedDomain(T, invmap, checkdomain(domain))
# - in that case, if the map is a Map{S}, make sure that S matches T
MappedDomain{T}(invmap::Map{T}, domain) where {T} =
    _MappedDomain(T, invmap, checkdomain(domain))
MappedDomain{T}(invmap::Map{S}, domain) where {S,T} =
    MappedDomain{T}(convert(Map{T}, invmap), domain)
# invoke the constructor
_MappedDomain(::Type{T}, invmap, domain) where {T} =
    MappedDomain{T,typeof(invmap),typeof(domain)}(invmap, domain)

similardomain(d::MappedDomain, ::Type{T}) where {T} =
    MappedDomain{T}(d.invmap, d.domain)

forward_map(d::MappedDomain) = rightinverse(d.invmap)
forward_map(d::MappedDomain, x) = rightinverse(d.invmap, x)

inverse_map(d::MappedDomain) = d.invmap
inverse_map(d::MappedDomain, y) = d.invmap(y)

"Map a domain with the inverse of the given map."
map_domain(map, domain) = map_domain1(map, domain)
map_domain1(map, domain) = map_domain2(map, domain)
map_domain2(map, domain) = default_map_domain(map, domain)
default_map_domain(map, domain) = mapped_domain(leftinverse(map), domain)

isequaldomain(a::MappedDomain, b::MappedDomain) =
    isequalmap(a.invmap, b.invmap) && isequaldomain(superdomain(a), superdomain(b))
domainhash(d::MappedDomain, h::UInt) = hashrec("MappedDomain", d.invmap, hash(superdomain(d)))

"Make a mapped domain with the given inverse map."
mapped_domain(invmap, domain) = mapped_domain1(invmap, domain)
mapped_domain1(invmap, domain) = mapped_domain2(invmap, domain)
mapped_domain2(invmap, domain) = default_mapped_domain(invmap, domain)

function default_mapped_domain(invmap, domain)
    T = promote_type(codomaintype(invmap), domaineltype(domain))
    MappedDomain{T}(convert_codomaintype(T, invmap), domain)
end

# We make a special case for linear mappings like 2 .* d, when d contains vectors
mapped_domain1(invmap::LinearMap{<:Number}, domain::Domain{T}) where {T<:AbstractVector} =
    mapped_domain(Map{T}(invmap), domain)
# allow some mixing of SVector and Vector
mapped_domain1(invmap::Map{SVector{N,T}}, domain) where {N,T} =
    _vector_mapped_domain1(invmap, domain)
_vector_mapped_domain1(invmap::Map{SVector{N,T}}, domain::Domain{Vector{U}}) where {N,T,U} =
    mapped_domain(invmap, convert(Domain{SVector{N,promote_type(T,U)}}, domain))
_vector_mapped_domain1(invmap::Map{SVector{N,T}}, domain) where {N,T} =
    mapped_domain2(invmap, domain)

mapped_domain2(invmap, domain::Domain{SVector{N,T}}) where {N,T} =
    _vector_mapped_domain2(invmap, domain)
_vector_mapped_domain2(invmap::Map{Vector{U}}, domain::Domain{SVector{N,T}}) where {N,T,U} =
    mapped_domain(convert(Map{SVector{N,promote_type(T,U)}}, invmap), domain)
_vector_mapped_domain2(invmap, domain::Domain{SVector{N,T}}) where {N,T} =
    default_mapped_domain(invmap, domain)

# Avoid nested mapping domains, construct a composite map instead
mapped_domain2(invmap, d::MappedDomain) = mapped_domain(inverse_map(d) ∘ invmap, superdomain(d))

boundary(d::MappedDomain) = _boundary(d, boundary(superdomain(d)), inverse_map(d))
_boundary(d::MappedDomain, superbnd, invmap) = MappedDomain(invmap, superbnd)
_boundary(d::MappedDomain, superbnd::UnionDomain, invmap) =
    UnionDomain(map(t->mapped_domain(invmap, t), components(superbnd)))

interior(d::MappedDomain) = _interior(d, superdomain(d), inverse_map(d))
_interior(d::MappedDomain, superdomain, invmap) = MappedDomain(invmap, interior(superdomain))
closure(d::MappedDomain) = _closure(d, superdomain(d), inverse_map(d))
_closure(d::MappedDomain, superdomain, invmap) = MappedDomain(invmap, closure(superdomain))

convert(::Type{MappedDomain}, d::Domain{T}) where {T} =
    MappedDomain{T}(mapto_canonical(d), canonicaldomain(d))
convert(::Type{MappedDomain{T}}, d::Domain) where {T} =
    MappedDomain{T}(mapto_canonical(d), canonicaldomain(d))


"A `ParametricDomain` stores the forward map of a mapped domain."
struct ParametricDomain{T,F,D} <: AbstractMappedDomain{T}
    fmap    ::  F
    domain  ::  D
end

ParametricDomain(fmap, domain) =
    ParametricDomain{codomaintype(fmap)}(fmap, domain)
ParametricDomain{T}(fmap, domain::Domain) where {T} =
    ParametricDomain{T,typeof(fmap),typeof(domain)}(fmap, domain)
ParametricDomain{T}(fmap, domain) where {T} =
    _ParametricDomain(T, fmap, checkdomain(domain))
_ParametricDomain(::Type{T}, fmap, domain) where {T} =
    ParametricDomain{T,typeof(fmap),typeof(domain)}(fmap, domain)

similardomain(d::ParametricDomain, ::Type{T}) where {T} =
    ParametricDomain{T}(d.fmap, d.domain)

forward_map(d::ParametricDomain) = d.fmap
forward_map(d::ParametricDomain, x) = d.fmap(x)

function indomain(x, d::ParametricDomain)
    # To check for membership, we can't use the inverse map because it may not exist
    # We assume a left inverse exists, but the left inverse may be many-to-one.
    # So we also have to check whether the left-inverse-point maps back to x
    y = leftinverse(d.fmap, x)
    x2 = forward_map(d, y)
    isapprox(x, x2)
end

isequaldomain(d1::ParametricDomain, d2::ParametricDomain) =
    isequalmap(forward_map(d1), forward_map(d2)) && isequaldomain(superdomain(d1), superdomain(d2))
domainhash(d::ParametricDomain, h::UInt) = hashrec("ParametricDomain", forward_map(d), superdomain(d), h)

"Return the domain that results from mapping the given domain."
parametric_domain(fmap, domain) = ParametricDomain(fmap, domain)
parametric_domain(fmap, domain::ParametricDomain) =
    parametric_domain(fmap ∘ forward_map(domain), superdomain(domain))

boundary(d::ParametricDomain) = _boundary(d, boundary(superdomain(d)), forward_map(d))
_boundary(d::ParametricDomain, superbnd, fmap) = ParametricDomain(fmap, superbnd)
_boundary(d::ParametricDomain, superbnd::UnionDomain, fmap) =
    UnionDomain(map(t -> parametric_domain(fmap, t), components(superbnd)))

interior(d::ParametricDomain) = _interior(d, superdomain(d), forward_map(d))
_interior(d::ParametricDomain, superdomain, fmap) = ParametricDomain(fmap, interior(superdomain))
closure(d::ParametricDomain) = _closure(d, superdomain(d), forward_map(d))
_closure(d::ParametricDomain, superdomain, fmap) = ParametricDomain(fmap, closure(superdomain))
