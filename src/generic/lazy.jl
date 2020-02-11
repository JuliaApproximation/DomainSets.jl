
"""
A lazy domain evaluates its membership function on the fly in terms of that of
other domains.

The `in(x, domain::LazyDomain)` applies three types of transformations:
1. Preprocess: `y = preprocess(domain, x)`
2. Distribution of `y` over member domains given by `elements(domain)`
3. Combination of the outputs into a single boolean result.

The distribution step is determined by the result of `composition(domain)`,
see `composition`. The combination is performed by `combine`.
"""
abstract type LazyDomain{T} <: Domain{T} end

elements(d::LazyDomain) = d.domains

"""
Supertype of all compositions of a lazy domain. The composition determines how
the point `x` is distributed to the member domains of a lazy domain.

Three compositions implemented in the package are:
- `NoComposition`: the lazy domain encapsulates a single domain and `x` is passed
	through unaltered
- `Combination`: the lazy domain has several members and `x` is passed to the `in`
	method of all members
- `Product`: the lazy domain has several members and the components of `x` are
	passed to the components of the lazy domain
"""
abstract type Composition end

struct NoComposition <: Composition end
struct Combination <: Composition end
struct Product <: Composition end

preprocess(d::LazyDomain, x) = x
composition(d::LazyDomain) = NoComposition()

indomain(x, d::LazyDomain) = _indomain(preprocess(d, x), d, composition(d), elements(d))
_indomain(x, d, ::NoComposition, domains) = in(x, domains[1])
_indomain(x, d, ::Combination, domains) = combine(d, map(d->in(x, d), domains))
if VERSION < v"1.2"
	_indomain(x, d, ::Product, domains) = reduce(&, map(in, x, domains))
else
	_indomain(x, d, ::Product, domains) = mapreduce(in, &, x, domains)
end

approx_indomain(x, d::LazyDomain, tolerance) =
	_approx_indomain(preprocess(d, x), d, tolerance, composition(d), elements(d))

_approx_indomain(x, d, tolerance, ::NoComposition, domains) =
    approx_in(x, domains[1], tolerance)
_approx_indomain(x, d, tolerance, ::Combination, domains) =
    combine(d, map(d -> approx_in(x, d, tolerance), domains))
if VERSION < v"1.2"
	_approx_indomain(x, d, tolerance, ::Product, domains) =
	    reduce(&, map((u,v)->approx_in(u, v, tolerance), x, domains))
else
	_approx_indomain(x, d, tolerance, ::Product, domains) =
	    mapreduce((u,v)->approx_in(u, v, tolerance), &, x, domains)
end


==(a::D, b::D) where {D<:LazyDomain} = elements(a) == elements(b)

"""
Combine the outputs of `in` of member domains into a single output of the lazy
domain.
"""
combine

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
