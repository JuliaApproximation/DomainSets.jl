# Equality of domains

A domain behaves as much as possible like the mathematical set it represents, irrespective of its type. Thus, for example, two domains are considered
equal if their membership functions agree.

It is not always possible to realize this intended behaviour in practice. Indeed
it may be difficult to discover automatically whether two domains are equal,
especially when their types are different. Deciding whether two domains are equal requires supporting implementation, hence the outcome is not always accurate.
Still, the principle serves as a design goal.

## The `isequaldomain` function

Equality of domains is decided by the `isequaldomain` function in general, and
simply by `==` for subtypes of the `Domain` supertype.
```julia
julia> UnitInterval() == UnitSimplex{Float64}() == 0..1
true

julia> ChebyshevInterval() == UnitBall{Float64}() == -1 .. 1
true
```
The methodology behind this example is explained in the section on
[Canonical domains](@ref).

!!! note
    If two domains are equal according to `==`, then their hashes are also equal.
    This allows the use of domains as keys in a `Dict`, if one is so inclined.


```@docs; canonical=false
isequaldomain
```


## Types versus mathematical sets

In some cases the notion of a set is more interesting than the functionality
provided by the membership function. A concrete type offers a representation of
a notion. Having a type for the real numbers (see [`ℝ`](@ref)) allows one to
specify the domain of a function of a real variable.

A type can also be convenient for other reasons. Although a ball with center and radius can be represented in terms of an affine map and a unit ball, not all affine maps
map the unit ball to a ball. Thus, in Julia, such a representation does not allow one to dispatch on the property of being a ball.

Yet, having more types also comes with disadvantages. Dispatching on the type of
a domain such as an `Interval` does not make a function apply to all intervals - because a scalar unit ball has a different type but mathematically is also an interval.

This tension is not resolved in the package. Still, deciding whether two domains are equal is a relevant aspect.
