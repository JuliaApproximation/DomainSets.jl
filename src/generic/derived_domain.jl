# derived_domain.jl
# Generic code for a domain that is implemented in terms of an underlying domain

"""
DerivedDomain is an abstract supertype for domains that are implemented in terms
of another domain using composition. DerivedDomain transfers the interface of
Domain to this superdomain. Concrete subtypes may modify any function, but without
modifications the domain acts exactly like the superdomain.
"""
abstract type DerivedDomain{T} <: Domain{T}
end

# We assume that the underlying domain is stored in a field called superdomain
superdomain(d::DerivedDomain) = d.superdomain

indomain(x, d::DerivedDomain) = indomain(x, superdomain(d))
