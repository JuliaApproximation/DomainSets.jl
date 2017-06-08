# mapped_domain.jl

################
# Preliminaries
################


#######################
# Main type definition
#######################


"""
A `MappedDomain` consists of a domain `d` and a map `f` such that
```
in(x, d::MappedDomain) = in(f(x), d)
```

If the function `f` maps a variable of type `S` to a variable of type `T`, then
the underlying domain should have type `S` and the mapped domain has type `T`.
"""
struct MappedDomain{D,F,T} <: DerivedDomain{T}
    superdomain ::  D
    f           ::  F
end

MappedDomain(domain::Domain{S}, f) where {S} =
    MappedDomain{typeof(domain),typeof(f),return_type(f,S)}(domain, f)

mapping(d::MappedDomain) = d.f

indomain(x, d::MappedDomain) = indomain(mapping(f) * x, superdomain(d))

applymap(f, domain::Domain) = MappedDomain(domain, f)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
applymap(f, domain::MappedDomain) = applymap(superdomain(domain), f ∘ mapping(domain))

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", superdomain(d))
