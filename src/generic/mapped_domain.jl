# mapped_domain.jl

################
# Preliminaries
################


#######################
# Main type definition
#######################

abstract type AbstractMappedDomain{T} <: Domain{T}
end

superdomain(d::AbstractMappedDomain) = d.superdomain


struct ForwardMappedDomain{D,F,T} <: AbstractMappedDomain{T}
    superdomain ::  D
    fwmap       ::  F
end

ForwardMappedDomain(d::Domain{S}, f::AbstractMap{T,S}) where {T,S} =
    ForwardMappedDomain{typeof(d), typeof(f), T}(d, f)

function indomain(x, d::ForwardMappedDomain)
    t = left_inverse(d.fwmap) * x
    if t ∈ superdomain(d)
        x2 = d.fwmap * t
        x == x2
    else
        false
    end
end

function approx_indomain(x, d::ForwardMappedDomain, tolerance)
    t = left_inverse(d.fwmap) * x
    if approx_indomain(t, superdomain(d), tolerance)
        x2 = d.fwmap * t
        norm(x-x2) <= tolerance
    else
        false
    end
end

forwardmap_domain(fwmap, domain::Domain) = ForwardMappedDomain(domain, fwmap)


"""
A `MappedDomain` consists of a domain `d` and a map `f` such that
```
in(x, d::MappedDomain) = in(f(x), d)
```

If the function `f` maps a variable of type `S` to a variable of type `T`, then
the underlying domain should have type `S` and the mapped domain has type `T`.
"""
struct MappedDomain{D,F,T} <: AbstractMappedDomain{T}
    superdomain ::  D
    f           ::  F
end

MappedDomain(domain::Domain{T}, f::AbstractMap{T,S}) where {S,T} =
    MappedDomain{typeof(domain),typeof(f),S}(domain, f)

function MappedDomain(domain::Domain{T}, f) where {T}
    S = domaintype(f)
    return_type(f, S) == T || error("Return type ", return_type(f,S), " of map ", f, " does not match element type ", T, " of domain ", domain, ".")
    MappedDomain{typeof(domain),typeof(f),S}(domain, f)
end

mapping(d::MappedDomain) = d.f

indomain(x, d::MappedDomain) = indomain(mapping(d) * x, superdomain(d))

approx_indomain(x, d::MappedDomain, tolerance) = approx_indomain(mapping(d) * x, superdomain(d), tolerance)

map_domain(f, domain::Domain) = MappedDomain(domain, f)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
map_domain(f, domain::MappedDomain) = map_domain(mapping(domain) ∘ f, superdomain(domain))

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", superdomain(d))
