# mapped_domain.jl

##################
# Mapped domains
##################

"""
A MappedDomain consists of a domain and a bidirectional map. The forward map
maps the domain onto the mapped domain, and the inverse map maps it back.
A point lies in the mapped domain, if the inverse map of that point lies in the
original domain.
"""
struct MappedDomain{DOMAIN <: Domain,MAP,T} <: Domain{T}
    domain  ::  DOMAIN
    # The forward map, from the underlying domain to the mapped domain
    fmap    ::  MAP

    # With this inner constructor we enforce that N is the dimension of the domain
    MappedDomain{DOMAIN,MAP,T}(domain::Domain{T}, fmap) where {DOMAIN <: Domain,MAP,T} = new{DOMAIN,MAP,T}(domain, fmap)
end

MappedDomain(domain::Domain{T}, fmap) where {T} =
    MappedDomain{typeof(domain),typeof(fmap),T}(domain, fmap)

domain(d::MappedDomain) = d.domain

mapping(d::MappedDomain) = d.fmap

indomain(x, d::MappedDomain) = indomain(inverse_map(mapping(d), x), domain(d))

apply_map(domain::Domain, map::AbstractMap) = MappedDomain(domain, map)

apply_map(d::MappedDomain, map::AbstractMap) = MappedDomain(domain(d), map*mapping(d))

(*)(map::AbstractMap, domain::Domain) = apply_map(domain, map)

(*)(domain::Domain, a::Number) = scaling_map(a*diagm(ones(eltype(domain)))) * domain

# TODO: revise
(+)(d::Domain, x::SVector{N,T}) where {N,T} = AffineMap(eye(SMatrix{N,N,T}),x) * d
# (+){N}(d::Domain{N}, x::AbstractVector) = d + SVector{N}(x)

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", domain(d))
