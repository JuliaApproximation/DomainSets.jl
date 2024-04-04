# Equality of domains

A domain behaves as much as possible like the mathematical set it represents, irrespective of its type. Thus, for example, two domains are considered
equal if their membership functions agree.

It is not always possible to realize this intended behaviour in practice. Indeed
it may be difficult to discover automatically whether two domains are equal,
especially when their types are different. Deciding whether two domains are equal requires supporting implementation, hence the outcome is not always accurate.
Still, the principle serves as a design goal.

```julia
julia> UnitInterval() == UnitSimplex{Float64}() == 0..1
true

julia> ChebyshevInterval() == UnitBall{Float64}() == -1 .. 1
true
```

Equality of domains is decided by the `isequaldomain` function in general, and
simply by `==` for subtypes of the `Domain` supertype.

```@docs; canonical=false
isequaldomain
```

!!! note
    If two domains are equal according to `==`, then their hashes are also equal.
    This allows the use of domains as keys in a `Dict`, if one is so inclined.
