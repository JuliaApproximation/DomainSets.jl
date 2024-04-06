
"""
    canonicaldomain([ctype::CanonicalType, ]d)

Return an associated canonical domain, if any, of the given domain.

For example, the canonical domain of an Interval `[a,b]` is the interval `[-1,1]`.

Optionally, a canonical type argument may specify an alternative canonical domain.
Canonical domains help with establishing equality between domains, with finding
maps between domains and with finding parameterizations.

If a domain implements a canonical domain, it should also implement
`mapfrom_canonical` and `mapto_canonical`.

See also: [`mapfrom_canonical`](@ref), [`mapto_canonical`](@ref).
"""
canonicaldomain(d) = d

"""
    hascanonicaldomain([ctype::CanonicalType, ]d)

Does the domain have a canonical domain?

See also: [`canonicaldomain`](@ref).
"""
hascanonicaldomain(d) = !(d === canonicaldomain(d))

identitymap(d) = IdentityMap{domaineltype(d)}(dimension(d))

"""
    mapfrom_canonical(d[, x])

Return a map to a domain `d` from its canonical domain.

If a second argument `x` is given, the map is evaluated at that point.
The point `x` should be a point in the canonical domain of `d`, and the
result is a point in `d`.

See also: [`mapto_canonical`](@ref), [`canonicaldomain`](@ref).
"""
mapfrom_canonical(d) = identitymap(d)
mapfrom_canonical(d, x) = mapfrom_canonical(d)(x)

"""
    mapto_canonical(d[, x])

Return a map from a domain `d` to its canonical domain.

If a second argument `x` is given, the map is evaluated at that point.
The point `x` should be a point in the domain `d`, and the
result is a point in the canonical domain.

See also: [`mapfrom_canonical`](@ref), [`canonicaldomain`](@ref).
"""
mapto_canonical(d) = leftinverse(mapfrom_canonical(d))
mapto_canonical(d, x) = mapto_canonical(d)(x)

# For definition of CanonicalType: see FunctionMaps

canonicaldomain(ctype::CanonicalType, d) = d
hascanonicaldomain(ctype::CanonicalType, d) = !(d === canonicaldomain(ctype, d))

mapfrom_canonical(ctype::CanonicalType, d) = identitymap(d)
mapto_canonical(ctype::CanonicalType, d) = leftinverse(mapfrom_canonical(ctype, d))

mapfrom_canonical(ctype::CanonicalType, d, x) = mapfrom_canonical(ctype, d)(x)
mapto_canonical(ctype::CanonicalType, d, x) = mapto_canonical(ctype, d)(x)


canonicaldomain(::Equal, d) = d
mapfrom_canonical(::Equal, d) = identitymap(d)
mapto_canonical(::Equal, d) = leftinverse(mapfrom_canonical(Equal(), d))

"""
    equaldomain(domain)

Return a canonical domain that is equal, but simpler. For example,
a 1-dimensional ball is an interval. This function is equivalent to
`canonicaldomain(Equal(), d)`.

A domain and its `equaldomain` are always equal domains according to
`isequaldomain`.

See also: [`canonicaldomain`](@ref).
"""
equaldomain(d) = canonicaldomain(Equal(), d)

"""
    hasequaldomain(d)

Does the domain have a known equal domain?
"""
hasequaldomain(d) = hascanonicaldomain(Equal(), d)

# these definitions should be removed because equal domains are never mapped
@deprecate mapfrom_equaldomain(d) mapfrom_canonical(Equal(), d)
@deprecate mapto_equaldomain(d) mapto_canonical(Equal(), d)
@deprecate mapfrom_equaldomain(d, x) mapfrom_canonical(Equal(), d, x)
@deprecate mapto_equaldomain(d, x) mapto_canonical(Equal(), d, x)

"""
    simplify(domain)

Simplify the given domain to an equal domain.
"""
function simplify(d)
    if hasequaldomain(d)
        equaldomain(d)
    elseif hasparameterization(d) && hasequaldomain(parameterdomain(d))
        map_domain(mapfrom_parameterdomain(d), simplify(parameterdomain(d)))
    else
        d
    end
end

"""
    simplifies(domain)

Does the domain simplify?
"""
simplifies(d) = hasequaldomain(d) ||
    (hasparameterization(d) && hasequaldomain(parameterdomain(d)))

"Convert the given domain to a domain defined in DomainSets.jl."
todomainset(d::Domain) = d
todomainset(d::DomainRef) = todomainset(domain(d))


"""
    Isomorphic <: CanonicalType

An isomorphic canonical domain is a domain that is the same up to an isomorphism.
"""
struct Isomorphic <: CanonicalType end

canonicaldomain(::Isomorphic, d) = canonicaldomain(Equal(), d)
mapfrom_canonical(::Isomorphic, d) = mapfrom_canonical(Equal(), d)
mapto_canonical(::Isomorphic, d) = leftinverse(mapfrom_canonical(Isomorphic(), d))

canonicaldomain(::Isomorphic, d::Domain{<:StaticVector{1,T}}) where {T} =
    convert(Domain{T}, d)
mapfrom_canonical(::Isomorphic, d::Domain{<:StaticVector{1,T}}) where {T} =
    NumberToVector{T}()

canonicaldomain(::Isomorphic, d::Domain{NTuple{N,T}}) where {N,T} =
    convert(Domain{SVector{N,T}}, d)
mapfrom_canonical(::Isomorphic, d::Domain{NTuple{N,T}}) where {N,T} =
    VectorToTuple{N,T}()



"""
    Parameterization <: CanonicalType

A parametric canonical domain can be parameterized from a simpler domain.
"""
struct Parameterization <: CanonicalType end

canonicaldomain(ctype::Parameterization, d) =
    hascanonicaldomain(d) ? canonicaldomain(ctype, canonicaldomain(d)) : d
mapfrom_canonical(ctype::Parameterization, d) =
    hascanonicaldomain(d) ? mapfrom_canonical(d) ∘ mapfrom_canonical(ctype, canonicaldomain(d)) : mapfrom_canonical(d)
mapto_canonical(ctype::Parameterization, d) = leftinverse(mapfrom_canonical(ctype, d))

# We define some convenience functions:
"""
    parameterdomain(d)

A domain from which `d` can be parameterized.
"""
parameterdomain(d) = canonicaldomain(Parameterization(), d)

"""
    parameterization(d)

A map from the parameter domain of `d` to `d`.
"""
parameterization(d) = mapfrom_canonical(Parameterization(), d)

"""
    hasparameterization(d)

Does the domain have a parameterization?
"""
hasparameterization(d) = hascanonicaldomain(Parameterization(), d)

"""
    mapfrom_parameterdomain(d[, x])

Convenience alias for `mapfrom_canonical(Paramaterization(), d[, x])`.
"""
mapfrom_parameterdomain(d) = mapfrom_canonical(Parameterization(), d)

"""
    mapto_parameterdomain(d[, x])

Convenience alias for `mapto_canonical(Paramaterization(), d[, x])`.

See also: [`mapfrom_parameterdomain`](@ref).
"""
mapto_parameterdomain(d) = mapto_canonical(Parameterization(), d)
mapfrom_parameterdomain(d, x) = mapfrom_canonical(Parameterization(), d, x)
mapto_parameterdomain(d, x) = mapto_canonical(Parameterization(), d, x)


"""
    mapto(d1, d2)

Return a map from domain `d1` to domain `d2`.
"""
mapto(d1, d2) = mapto1(d1, d2)
mapto(d1::D, d2::D) where {D} = d1 == d2 ? identitymap(d1) : mapto1(d1,d2)

# simplify the first argument
mapto1(d1, d2) =
    hasparameterization(d1) ? mapto(parameterdomain(d1), d2) ∘ mapto_parameterdomain(d1) : mapto2(d1, d2)
# simplify the second argument
mapto2(d1, d2) =
    hasparameterization(d2) ? mapfrom_parameterdomain(d2) ∘ mapto(d1, parameterdomain(d2)) : default_mapto(d1,d2)

default_mapto(d1, d2) = d1 == d2 ? identitymap(d1) : throw(ArgumentError("No map known between $(d1) and $(d2)."))



## Equality for domains
# We generically define the equality of Domains, using the framework
# of canonical domains.

"""
    isequaldomain(d1, d2)

Are the two given domains equal?

Domains are considered equal if their membership functions return the same
output for the same input.

It is not always possible to verify this automatically. If the result is `true`,
then the domains are guaranteed to be equal. If the result is `false`, then
either the domains are not equal or they are equal but the implementation fails
to recognize this.
"""
isequaldomain(d1, d2) = isequaldomain1(d1, d2)
isequaldomain1(d1, d2) = simplifies(d1) ? isequaldomain(simplify(d1), d2) : isequaldomain2(d1, d2)
isequaldomain2(d1, d2) = simplifies(d2) ? isequaldomain(d1, simplify(d2)) : default_isequaldomain(d1, d2)
default_isequaldomain(d1, d2) = d1 === d2

# Can't use == for arrays etcetera since order doesn't matter
isequaldomain(d1::BaseDomainType, d2::BaseDomainType) = d1 ⊆ d2 && d2 ⊆ d1

==(d1::AnyDomain, d2::AnyDomain) = isequaldomain(domain(d1), domain(d2))

# Associated with == we have to define the hashes of domains
# Since equality is based on `simplify`, so should hash be. To that end we
# are forced to intercept hash to call simplify. In order to avoid a stack
# overflow we call a different function `domainhash`, for which we provide
# a fallback by invoking the implementation in Base.
hash(d::Domain, h::UInt) = domainhash(simplify(d), h)

domainhash(d) = domainhash(d, zero(UInt))
domainhash(d, h::UInt) = invoke(hash, Tuple{Any,UInt}, d, h)
