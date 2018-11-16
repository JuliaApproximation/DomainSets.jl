
"""
A `WrappedDomain` is a wrapper around an object that implements the domain
interface, and that is itself a domain.
"""
struct WrappedDomain{D,T} <: Domain{T}
    domain  ::  D
end

WrappedDomain(domain) = WrappedDomain{typeof(domain),eltype(domain)}(domain)

indomain(x, d::WrappedDomain) = in(x, d.domain)

# Anything can be converted to a domain by wrapping it. An error will be thrown
# if the object does not support `eltype`.
convert(::Type{Domain}, v::Domain) = v
convert(::Type{Domain}, v) = WrappedDomain(v)
