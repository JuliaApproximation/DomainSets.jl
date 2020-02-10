
"""
A lazy domain evaluates its membership function on the fly in terms of that of
other domains.
"""
abstract type LazyDomain{T} <: Domain{T} end

elements(d::LazyDomain) = d.domains

abstract type Composition end
struct NoComposition <: Composition end
struct Combination <: Composition end
struct Product <: Composition end

preprocess(d::LazyDomain, x) = x
composition(d::LazyDomain) = NoComposition()

indomain(x, d::LazyDomain) = _indomain(preprocess(d, x), d, composition(d), elements(d))
_indomain(x, d, ::NoComposition, domains) = in(x, domains[1])
_indomain(x, d, ::Combination, domains) = combine(d, map(d->in(x, d), domains))
_indomain(x, d, ::Product, domains) = mapreduce(in, &, x, domains)

approx_indomain(x, d::LazyDomain, tolerance) =
	_approx_indomain(preprocess(d, x), d, tolerance, composition(d), elements(d))

_approx_indomain(x, d, tolerance, ::NoComposition, domains) =
    approx_in(x, domains[1], tolerance)
_approx_indomain(x, d, tolerance, ::Combination, domains) =
    combine(d, map(d -> approx_in(x, d, tolerance), domains))
_approx_indomain(x, d, tolerance, ::Product, domains) =
    mapreduce((u,v)->approx_in(u, v, tolerance), &, x, domains)

==(a::D, b::D) where {D<:LazyDomain} = elements(a) == elements(b)


"""
A `WrappedDomain` is a wrapper around an object that implements the domain
interface, and that is itself a domain.
"""
struct WrappedDomain{T,D} <: LazyDomain{T}
    domain  ::  D
end

elements(d::WrappedDomain) = (d.domain,)

WrappedDomain(domain::D) where {T,D<:Domain{T}} = WrappedDomain{T,D}(domain)
WrappedDomain(domain::D) where {D} = WrappedDomain{eltype(D),D}(domain)

# Anything can be converted to a domain by wrapping it. An error will be thrown
# if the object does not support `eltype`.
convert(::Type{Domain}, v::Domain) = v
convert(::Type{Domain}, v) = WrappedDomain(v)
