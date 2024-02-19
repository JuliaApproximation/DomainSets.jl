
"""
    canonicaldomain([ctype::CanonicalType, ]domain)

Return an associated canonical domain, if any, of the given domain.

For example, the canonical domain of an Interval `[a,b]` is the interval `[-1,1]`.

Optionally, a canonical type argument may specify an alternative canonical domain.
Canonical domains help with establishing equality between domains, with finding
maps between domains and with finding parameterizations.

If a domain implements a canonical domain, it should also implement
`mapfrom_canonical` and `mapto_canonical`.
"""
canonicaldomain(d) = d

"Does the domain have a canonical domain?"
hascanonicaldomain(d) = !(d === canonicaldomain(d))

identitymap(d) = IdentityMap{domaineltype(d)}(dimension(d))

"""
    mapfrom_canonical(d[, x])

Return a map to a domain `d` from its canonical domain.

If a second argument `x` is given, the map is evaluated at that point.
The point `x` should be a point in the canonical domain of `d`, and the
result is a point in `d`.
"""
mapfrom_canonical(d) = identitymap(d)
mapfrom_canonical(d, x) = mapfrom_canonical(d)(x)

"Return a map from the domain to its canonical domain."
mapto_canonical(d) = leftinverse(mapfrom_canonical(d))
mapto_canonical(d, x) = mapto_canonical(d)(x)


"Supertype of kinds of canonical domains."
abstract type CanonicalType end

canonicaldomain(ctype::CanonicalType, d) = d
hascanonicaldomain(ctype::CanonicalType, d) = !(d === canonicaldomain(ctype, d))

mapfrom_canonical(ctype::CanonicalType, d) = identitymap(d)
mapto_canonical(ctype::CanonicalType, d) = leftinverse(mapfrom_canonical(ctype, d))

mapfrom_canonical(ctype::CanonicalType, d, x) = mapfrom_canonical(ctype, d)(x)
mapto_canonical(ctype::CanonicalType, d, x) = mapto_canonical(ctype, d)(x)


"A canonical domain that is equal but simpler (e.g. a 1-dimensional ball is an interval)."
struct Equal <: CanonicalType end

canonicaldomain(::Equal, d) = d
mapfrom_canonical(::Equal, d) = identitymap(d)
mapto_canonical(::Equal, d) = leftinverse(mapfrom_canonical(Equal(), d))

simplify(d) = canonicaldomain(Equal(), d)
simplifies(d) = hascanonicaldomain(Equal(), d)

"Convert the given domain to a domain defined in DomainSets.jl."
todomainset(d::Domain) = d
todomainset(d::DomainRef) = todomainset(domain(d))


"A canonical domain that is isomorphic but may have different element type."
struct Isomorphic <: CanonicalType end

canonicaldomain(::Isomorphic, d) = canonicaldomain(Equal(), d)
mapfrom_canonical(::Isomorphic, d) = mapfrom_canonical(Equal(), d)
mapto_canonical(::Isomorphic, d) = leftinverse(mapfrom_canonical(Isomorphic(), d))

canonicaldomain(::Isomorphic, d::Domain{SVector{1,T}}) where {T} =
    convert(Domain{T}, d)
mapfrom_canonical(::Isomorphic, d::Domain{SVector{1,T}}) where {T} =
    NumberToVector{T}()

canonicaldomain(::Isomorphic, d::Domain{NTuple{N,T}}) where {N,T} =
    convert(Domain{SVector{N,T}}, d)
mapfrom_canonical(::Isomorphic, d::Domain{NTuple{N,T}}) where {N,T} =
    VectorToTuple{N,T}()



"A parameter domain that can be mapped to the domain."
struct Parameterization <: CanonicalType end

canonicaldomain(ctype::Parameterization, d) =
    hascanonicaldomain(d) ? canonicaldomain(ctype, canonicaldomain(d)) : d
mapfrom_canonical(ctype::Parameterization, d) =
    hascanonicaldomain(d) ? mapfrom_canonical(d) ∘ mapfrom_canonical(ctype, canonicaldomain(d)) : mapfrom_canonical(d)
mapto_canonical(ctype::Parameterization, d) = leftinverse(mapfrom_canonical(ctype, d))

# We define some convenience functions:
"Return a parameter domain which supports a `parameterization`."
parameterdomain(d) = canonicaldomain(Parameterization(), d)

"Return a parameterization of the given domain."
parameterization(d) = mapfrom_canonical(Parameterization(), d)

"Does the domain have a parameterization?"
hasparameterization(d) = hascanonicaldomain(Parameterization(), d)

mapfrom_parameterdomain(d) = mapfrom_canonical(Parameterization(), d)
mapto_parameterdomain(d) = mapto_canonical(Parameterization(), d)
mapfrom_parameterdomain(d, x) = mapfrom_canonical(Parameterization(), d, x)
mapto_parameterdomain(d, x) = mapto_canonical(Parameterization(), d, x)


"Return a map from domain `d1` to domain `d2`."
mapto(d1, d2) = mapto1(d1, d2)
mapto(d1::D, d2::D) where {D} = d1 == d2 ? identitymap(d1) : mapto1(d1,d2)

# simplify the first argument
mapto1(d1, d2) =
    hasparameterization(d1) ? mapto(parameterdomain(d1), d2) ∘ mapto_parameterdomain(d1) : mapto2(d1, d2)
# simplify the second argument
mapto2(d1, d2) =
    hasparameterization(d2) ? mapfrom_parameterdomain(d2) ∘ mapto(d1, parameterdomain(d2)) : no_known_mapto(d1,d2)

no_known_mapto(d1, d2) = d1 == d2 ? identitymap(d1) : error("No map known between $(d1) and $(d2).")



## Equality for domains
# We generically define the equality of Domains, using the framework
# of canonical domains.

@deprecate isequal1(d1,d2) isequaldomain1(d1,d2)
@deprecate isequal2(d1,d2) isequaldomain2(d1,d2)

"Are the two given domains equal?"
isequaldomain(d1, d2) = isequaldomain1(d1, d2)
# simplify the first argument
isequaldomain1(d1, d2) = simplifies(d1) ? isequaldomain(simplify(d1), d2) : isequaldomain2(d1, d2)
# simplify the second argument
isequaldomain2(d1, d2) = simplifies(d2) ? isequaldomain(d1, simplify(d2)) : d1===d2

==(d1::AnyDomain, d2::AnyDomain) = isequaldomain(domain(d1), domain(d2))
