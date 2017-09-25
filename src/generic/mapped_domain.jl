# mapped_domain.jl

################
# Preliminaries
################


#######################
# Main type definition
#######################

"""
There are three objects involved in the mapping of a domain:
- the original domain (denoted src)
- the map (denoted f)
- the resulting domain (target)

If `f` maps a variable of type `S` to a variable of type `T`, then the `src`
domain has eltype `S` and the target domain has eltype `T`.

The characteristic function of the resulting domain is defined in terms of the
inverse of the map `f`, i.e.:
```
in(x, target) = in(inv(f)*x, src)
```

Concrete mapped domains can be implemented in various ways, e.g. by storing `src`
and `f`, or by storing `src` and `inv(f)`, ...
"""
abstract type MappedDomain{D,T} <: Domain{T}
end

src(d::MappedDomain) = d.src

show(io::IO, d::MappedDomain) =  print(io, "A mapped domain based on ", src(d))

point_in_domain(d::MappedDomain) = forward_map(d) * point_in_domain(src(d))

"""
A `ForwardMappedDomain` stores the `src` domain and the forward map `f`, which
maps `src` to `target`.
"""
struct ForwardMappedDomain{F,D,T} <: MappedDomain{D,T}
    src     ::  D
    fwmap   ::  F
end

ForwardMappedDomain(src::Domain{S}, fwmap::AbstractMap{T,S}) where {T,S} =
    ForwardMappedDomain{typeof(fwmap),typeof(src),T}(src, fwmap)

forward_map(d::ForwardMappedDomain) = d.fwmap
inverse_map(d::ForwardMappedDomain) = inv(d.fwmap)

# Rationale:
# The point x that is given can be any point in T. Yet, the map from S to T does
# not have to be fully invertible on T. It is sufficient that it is invertible
# on its range, and that it has a left inverse that is defined on T.
# that is defined on the space of T.
# Then we can proceed as follows:
# - we compute the left inverse of the point x ∈ T and check that the resulting
#   point t ∈ S lies in the src domain
# - but that is not enough: then we map the point t back to T using the forward
#   map and check that it agrees with the original point x.
function indomain(x, d::ForwardMappedDomain)
    t = left_inverse(d.fwmap) * x
    if t ∈ src(d)
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
    if approx_indomain(t, src(d), tolerance)
        x2 = d.fwmap * t
        norm(x-x2) <= tolerance
    else
        false
    end
end

forwardmap_domain(fwmap, domain::Domain) = ForwardMappedDomain(domain, fwmap)

# Avoid chain of multiple mappings, use composite map instead
forwardmap_domain(fwmap, domain::ForwardMappedDomain) = forwardmap_domain(fwmap ∘ forward_map(domain), src(domain))


"""
An `InverseMappedDomain` stores the `src` and the inverse of the map `f`.
"""
struct InverseMappedDomain{F,D,T} <: MappedDomain{D,T}
    src     ::  D
    invmap  ::  F
end

InverseMappedDomain(src::Domain{T}, invmap::AbstractMap{T,S}) where {S,T} =
    InverseMappedDomain{typeof(invmap),typeof(src),S}(src, invmap)

function InverseMappedDomain(src::Domain{T}, invmap) where {T}
    S = domaintype(invmap)
    return_type(invmap, S) == T ||
        error("Return type ", return_type(invmap,S), " of map ", invmap, " does not match element type ", T, " of domain ", src, ".")
    InverseMappedDomain{typeof(invmap),typeof(src),S}(src, invmap)
end

forward_map(d::InverseMappedDomain) = inv(d.invwmap)
inverse_map(d::InverseMappedDomain) = d.invmap

indomain(x, d::InverseMappedDomain) = indomain(d.invmap * x, src(d))

approx_indomain(x, d::InverseMappedDomain, tolerance) = approx_indomain(d.invmap * x, src(d), tolerance)

inversemap_domain(invmap, src::Domain) = InverseMappedDomain(src, invmap)

# Avoid nested mapping domains, construct a composite map instead
# This assumes that the map types can be combined using \circ
inversemap_domain(invmap, d::MappedDomain) = inversemap_domain(inverse_map(d) ∘ invmap, src(d))



"""
A `BidirectionalMappedDomain` stores the `src` domain and both the map `f` and its inverse.
"""
struct BidirectionalMappedDomain{F1,F2,D,T} <: MappedDomain{D,T}
    src     ::  D
    fwmap   ::  F1
    invmap  ::  F2
end

BidirectionalMappedDomain(src::Domain, fwmap::AbstractMap) = BidirectionalMappedDomain(src, fwmap, inv(fwmap))

BidirectionalMappedDomain(src::Domain{S}, fwmap::AbstractMap{T,S}, invmap::AbstractMap{S,T}) where {S,T} =
    BidirectionalMappedDomain{typeof(fwmap),typeof(invmap),typeof(src),T}(src, fwmap, invmap)

forward_map(d::BidirectionalMappedDomain) = d.fwmap
inverse_map(d::BidirectionalMappedDomain) = d.invmap
