
"""
    canonicaldomain(d::Domain[, args...])

Return an associated canonical domain, if any, of the given domain.

Optionally, additional arguments may specify an alternative canonical domain.
Canonical domains help with establishing equality between domains, with finding
maps between domains and with finding parameterizations.

For example, the canonical domain of any interval `[a,b]` is the unit interval
`[0,1]`.
"""
canonicaldomain(d::Domain) = d
canonicaldomain(d, args...) = canonicaldomain(d)

identitymap(d) = IdentityMap{eltype(d)}(dimension(d))

"Return a map to a domain from its canonical domain."
fromcanonical(d) = identitymap(d)
fromcanonical(d, args...) = fromcanonical(d)

"Return a map from the domain to its canonical domain."
tocanonical(d) = leftinverse(fromcanonical(d))
tocanonical(d, args...) = leftinverse(fromcanonical(d, args...))


"Supertype of kinds of canonical domains."
abstract type CanonicalType end

"A canonical domain that is equal but simpler (e.g. a 1-dimensional ball is an interval)."
struct Equal <: CanonicalType end

canonicaldomain(d, ::Equal) = d
fromcanonical(d, ::Equal) = identitymap(d)
tocanonical(d, ::Equal) = leftinverse(tocanonical(d, Equal()))

simplify(d) = canonicaldomain(d, Equal())

"A canonical domain that is isomorphic but may have different element type."
struct Isomorphic <: CanonicalType end

canonicaldomain(d, ::Isomorphic) = canonicaldomain(d, Equal())
fromcanonical(d, ::Isomorphic) = fromcanonical(d, Equal())
tocanonical(d, ::Isomorphic) = leftinverse(tocanonical(d, Isomorphic()))

canonicaldomain(d::Domain{SVector{1,T}}, ::Isomorphic) where {T} =
    convert(Domain{T}, d)
fromcanonical(d::Domain{SVector{1,T}}, ::Isomorphic) where {T} =
    NumberToVector{T}()

canonicaldomain(d::Domain{NTuple{N,T}}, ::Isomorphic) where {N,T} =
    convert(Domain{SVector{N,T}}, d)
fromcanonical(d::Domain{NTuple{N,T}}, ::Isomorphic) where {N,T} =
    VectorToTuple{N,T}()

"A parameter domain that can be mapped to the domain."
struct Parameterization <: CanonicalType end

# We define some convenience functions:
"Return a parameter domain which supports a `parameterization`."
parameterdomain(d::Domain) = canonicaldomain(d, Parameterization())

"Return a parameterization of the given domain."
parameterization(d::Domain) = fromcanonical(d, Parameterization())

from_parameterdomain(d::Domain) = fromcanonical(d, Parameterization())
to_parameterdomain(d::Domain) = tocanonical(d, Parameterization())


"Return a map from domain `d1` to domain `d2`."
mapto(d1, d2) = mapto1(d1, d2)
mapto(d1::D, d2::D) where {D} = d1 == d2 ? identitymap(d1) : mapto1(d1,d2)

# simplify the first argument
mapto1(d1, d2) = _mapto1(d1, d2, parameterdomain(d1))
_mapto1(d1::D, d2, cd::D) where {D} = mapto2(d1, d2)
_mapto1(d1, d2, cd) = mapto(cd, d2) ∘ to_parameterdomain(d1)
# simplify the second argument
mapto2(d1, d2) = _mapto2(d1, d2, parameterdomain(d2))
_mapto2(d1, d2::D, cd::D) where {D} = no_known_mapto(d1, d2)
_mapto2(d1, d2, cd) = from_parameterdomain(d2) ∘ mapto(d1, cd)

no_known_mapto(d1, d2) = d1 == d2 ? identitymap(d1) : error("No map known between $(d1) and $(d2).")



## Equality for domains

==(d1::Domain, d2::Domain) = isequal1(d1, d2)
# simplify the first argument
isequal1(d1, d2) = _isequal1(d1, d2, simplify(d1))
_isequal1(d1::D, d2, sd::D) where {D} = isequal2(d1, d2)
_isequal1(d1, d2, sd) = sd==d2

# simplify the second argument
isequal2(d1, d2) = _isequal2(d1, d2, simplify(d2))
_isequal2(d1, d2::D, sd::D) where {D} = d1===d2
_isequal2(d1, d2, sd) = d1==sd
