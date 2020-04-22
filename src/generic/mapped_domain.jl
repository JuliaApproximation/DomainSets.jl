
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
<<<<<<< HEAD
dimension(d::MappedVectorDomain) = dimension(superdomain(d))
=======
dimension(d::MappedVectorDomain) = dimension(source(d))

# isopenset(d::MappedDomain) = isopenset(source(d))
# isclosedset(d::MappedDomain) = isclosedset(source(d))
#
# TODO: we can't really define open and close in general this way.
# We need more properties of the map to be able to conclude whether the mapped
# domain is open or closed.

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

# TODO: rethink maps and mapped domains and avoid using the default_tolerance
# function below.

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
>>>>>>> master

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
