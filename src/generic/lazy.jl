
"""
A lazy domain evaluates its membership function on the fly in terms of that of
other domains.

The `in(x, domain::LazyDomain)` applies three types of transformations:
1. Point mapping: `y = tointernalpoint(domain, x)`
2. Distribution of `y` over member domains given by `elements(domain)`
3. Combination of the outputs into a single boolean result.

The distribution step is determined by the result of `composition(domain)`,
see `composition`. The combination is performed by `combine`. Mapping between
points of the lazy domain and points of its member domains is described by
`y = tointernalpoint(domain, x)` and `x = toexternalpoint(domain, y)`.
"""
abstract type LazyDomain{T} <: Domain{T} end

"Translate a point of the lazy domain to a point (or points) of the composing domain."
tointernalpoint(d::LazyDomain, x) = x
"Inverse of `tointernalpoint`."
toexternalpoint(d::LazyDomain, y) = y


"""
A single lazy domain is defined in terms of a single domain.

It has no composition and no combination of its `in` function.
"""
abstract type SingleLazyDomain{T} <: LazyDomain{T} end

superdomain(d::SingleLazyDomain) = d.domain

indomain(x, d::SingleLazyDomain) = in(tointernalpoint(d, x), superdomain(d))
approx_indomain(x, d::SingleLazyDomain, tolerance) = approx_in(tointernalpoint(d, x), superdomain(d), tolerance)


"A composite lazy domain is defined in terms of multiple domains."
abstract type CompositeLazyDomain{T} <: LazyDomain{T} end

elements(d::CompositeLazyDomain) = d.domains

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
abstract type LazyComposition end

struct NoComposition <: LazyComposition end
struct Combination <: LazyComposition end
struct Product <: LazyComposition end

composition(d::CompositeLazyDomain) = NoComposition()

indomain(x, d::CompositeLazyDomain) = _indomain(tointernalpoint(d, x), d, composition(d), elements(d))
_indomain(x, d, ::NoComposition, domains) = in(x, domains[1])
_indomain(x, d, ::Combination, domains) = combine(d, map(d->in(x, d), domains))
_indomain(x, d, ::Product, domains) = mapreduce(in, &, x, domains)

approx_indomain(x, d::CompositeLazyDomain, tolerance) =
	_approx_indomain(tointernalpoint(d, x), d, tolerance, composition(d), elements(d))

_approx_indomain(x, d, tolerance, ::NoComposition, domains) =
    approx_in(x, domains[1], tolerance)
_approx_indomain(x, d, tolerance, ::Combination, domains) =
    combine(d, map(d -> approx_in(x, d, tolerance), domains))
_approx_indomain(x, d, tolerance, ::Product, domains) =
    mapreduce((u,v)->approx_in(u, v, tolerance), &, x, domains)

point_in_domain(d::SingleLazyDomain) = toexternalpoint(d, point_in_domain(superdomain(d)))
point_in_domain(d::CompositeLazyDomain) = toexternalpoint(d, map(point_in_domain, elements(d)))

==(a::D, b::D) where {D<:CompositeLazyDomain} = elements(a) == elements(b)


dimension(d::SingleLazyDomain{Vector{T}}) where {T} = dimension(superdomain(d))
function dimension(d::CompositeLazyDomain{Vector{T}}) where {T}
	dim = dimension(element(d,1))
	@assert all(isequal(dim), map(dimension, elements(d)))
	dim
end

"""
Combine the outputs of `in` of member domains into a single output of the lazy
domain.
"""
combine


"Abstract supertype for domains that wrap another domain."
abstract type DerivedDomain{T} <: SingleLazyDomain{T} end

isempty(d::DerivedDomain) = isempty(superdomain(d))

canonicaldomain(d::DerivedDomain) = superdomain(d)


"""
A `WrappedDomain` is a wrapper around an object that implements the domain
interface, and that is itself a domain.
"""
struct WrappedDomain{T,D} <: DerivedDomain{T}
    domain  ::  D
end

WrappedDomain(domain::Domain{T}) where {T} = WrappedDomain{T}(domain)
WrappedDomain(domain) = WrappedDomain{eltype(domain)}(domain)

WrappedDomain{T}(domain::D) where {T,D<:Domain{T}} = WrappedDomain{T,D}(domain)
WrappedDomain{T}(domain::Domain) where {T} = WrappedDomain{T}(convert(Domain{T}, domain))
WrappedDomain{T}(domain) where {T} = WrappedDomain{T,typeof(domain)}(domain)

similardomain(d::WrappedDomain, ::Type{T}) where {T} = WrappedDomain{T}(d.domain)

# Anything can be converted to a domain by wrapping it. An error will be thrown
# if the object does not support `eltype`.
convert(::Type{Domain}, v::Domain) = v
convert(::Type{Domain}, v) = WrappedDomain(v)

==(d1::WrappedDomain, d2::WrappedDomain) = superdomain(d1)==superdomain(d2)

"Example of a domain that wraps another domain and thus obtains its own type."
struct ExampleNamedDomain{T,D} <: DerivedDomain{T}
	domain	::	D
end
ExampleNamedDomain(domain::D) where {T,D<:Domain{T}} = ExampleNamedDomain{T,D}(domain)
