# derived_domain.jl
# Generic code for a domain that is implemented in terms of an underlying domain

"""
DerivedDomain is an abstract supertype for domains that are implemented in terms
of another domain using composition. DerivedDomain transfers the interface of
Domain to this superdomain.

A concrete subtype that inherits from DerivedDomain and stores a `superdomain`
is functionally equivalent to that superdomain. Any properties of the superdomain
can be modified by overriding a suitable function. For example, `in` of the
concrete domain may be implemented in terms of the `in` of the superdomain.
"""
abstract type DerivedDomain{T} <: Domain{T}
end

# We assume that the underlying domain is stored in a field called superdomain
superdomain(d::DerivedDomain) = d.superdomain

"Return the eltype of the superdomain."
supereltype(d::DerivedDomain) = eltype(superdomain(d))

indomain(x, d::DerivedDomain) = in(x, superdomain(d))
