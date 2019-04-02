# mapped_domain.jl

################
# Preliminaries
################


#######################
# Main type definition
#######################

"""
There are three objects involved in the mapping of a domain:
- the original domain (denoted source)
- the map (denoted f)
- the resulting domain (target)

If `f` maps a variable of type `S` to a variable of type `T`, then the `source`
domain has eltype `S` and the target domain has eltype `T`.

The characteristic function of the resulting domain is defined in terms of the
inverse of the map `f`, i.e.:
```
in(x, target) = in(inv(f)*x, source)
```

Concrete mapped domains can be implemented in various ways, e.g. by storing `source`
and `f`, or by storing `source` and `inv(f)`, ...
"""
abstract type MappedDomain{D,T} <: Domain{T}
end

source(d::MappedDomain) = d.source

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", source(d))

point_in_domain(d::MappedDomain) = forward_map(d) * point_in_domain(source(d))

elements(d::MappedDomain{D,T}) where {D<:ProductDomain,T} = map(x->forward_map(d)*x, elements(source(d)))

"""
A `ForwardMappedDomain` stores the `source` domain and the forward map `f`, which
maps `source` to `target`.
"""
struct ForwardMappedDomain{F,D,T} <: MappedDomain{D,T}
    source  ::  D
    fwmap   ::  F
end

ForwardMappedDomain(source::Domain, fwmap::AbstractMap{S,T}) where {S,T} =
    ForwardMappedDomain{typeof(fwmap),typeof(source),T}(source, fwmap)


forward_map(d::ForwardMappedDomain) = d.fwmap
inverse_map(d::ForwardMappedDomain) = inv(d.fwmap)

# Rationale:
# The point x that is given can be any point in T. Yet, the map from S to T does
# not have to be fully invertible on T. It is sufficient that it is invertible
# on its range, and that it has a left inverse that is defined on T.
# Then we can proceed as follows:
# - we compute the left inverse of the point x ∈ T and check that the resulting
#   point t ∈ S lies in the source domain
# - but that is not enough: then we map the point t back to T using the forward
#   map and check that it agrees with the original point x.
function indomain(x, d::ForwardMappedDomain)
    t = left_inverse(d.fwmap) * x
    if t ∈ source(d)
        x2 = d.fwmap * t
        # We can be strict in the first condition above, but have no choice but
        # to be approximate in the second one
        norm(x-x2) <= default_tolerance(d)
    else
        false
    end
end

# Reasoning is the same as above. However, for both tests we now use a tolerance.
# Note that the meaning of the tolerance may be different in space S and space T.
function approx_indomain(x, d::ForwardMappedDomain, tolerance)
    t = left_inverse(d.fwmap) * x
    if approx_indomain(t, source(d), tolerance)
        x2 = d.fwmap * t
        norm(x-x2) <= tolerance
    else
        false
    end
end

# we must map any point somewhere
isempty(d::ForwardMappedDomain) = isempty(d.source)

forwardmap_domain(fwmap, domain::Domain) = ForwardMappedDomain(domain, fwmap)

# Avoid chain of multiple mappings, use composite map instead
forwardmap_domain(fwmap, domain::ForwardMappedDomain) = forwardmap_domain(fwmap ∘ forward_map(domain), source(domain))


"""
An `InverseMappedDomain` stores the `source` and the inverse of the map `f`.
"""
struct InverseMappedDomain{F,D,T} <: MappedDomain{D,T}
    source  ::  D
    invmap  ::  F
end

InverseMappedDomain(source::Domain{T}, invmap::AbstractMap{S,T}) where {S,T} =
    InverseMappedDomain{typeof(invmap),typeof(source),S}(source, invmap)

function InverseMappedDomain(source::Domain{T}, invmap) where {T}
    S = domaintype(invmap)
    return_type(invmap, S) == T ||
        error("Return type ", return_type(invmap,S), " of map ", invmap, " does not match element type ", T, " of domain ", source, ".")
    InverseMappedDomain{typeof(invmap),typeof(source),S}(source, invmap)
end

forward_map(d::InverseMappedDomain) = inv(d.invmap)
inverse_map(d::InverseMappedDomain) = d.invmap

indomain(x, d::InverseMappedDomain) = in(d.invmap * x, source(d))

approx_indomain(x, d::InverseMappedDomain, tolerance) = approx_indomain(d.invmap * x, source(d), tolerance)

inversemap_domain(invmap, source::Domain) = InverseMappedDomain(source, invmap)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
inversemap_domain(invmap, d::MappedDomain) = inversemap_domain(inverse_map(d) ∘ invmap, source(d))



"""
A `BidirectionalMappedDomain` stores the `source` domain and both the map `f`
and a left inverse.
"""
struct BidirectionalMappedDomain{F1,F2,D,T} <: MappedDomain{D,T}
    source  ::  D
    fwmap   ::  F1
    invmap  ::  F2
end

BidirectionalMappedDomain(source::Domain, fwmap::AbstractMap) = BidirectionalMappedDomain(source, fwmap, inv(fwmap))

BidirectionalMappedDomain(source::Domain{S}, fwmap::AbstractMap{S,T}, invmap::AbstractMap{T,S}) where {S,T} =
    BidirectionalMappedDomain{typeof(fwmap),typeof(invmap),typeof(source),T}(source, fwmap, invmap)

forward_map(d::BidirectionalMappedDomain) = d.fwmap
inverse_map(d::BidirectionalMappedDomain) = d.invmap

# We proceed as we do for ForwardMappedDomain, except that we can use the stored
# left inverse. This should be refactored.
function indomain(x, d::BidirectionalMappedDomain)
    t = d.invmap * x
    if t ∈ source(d)
        x2 = d.fwmap * t
        norm(x-x2) <= default_tolerance(d)
    else
        false
    end
end

function approx_indomain(x, d::BidirectionalMappedDomain, tolerance)
    t = d.invmap * x
    if approx_indomain(t, source(d), tolerance)
        x2 = d.fwmap * t
        norm(x-x2) <= tolerance
    else
        false
    end
end
