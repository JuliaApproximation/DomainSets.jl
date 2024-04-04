# The domain interface

Existing types may add the interpretation of being a domain by implementing the
domain interface. They gain the ability to interact with other domains.

## The `in` function

A domain is a set of elements that is possibly continuous. Continuous sets are
defined mathematically, not by an exhaustive list of their elements. In practice membership of the set is defined by the implementation of `in`. The function
call `in(x, domain)` evaluates to true if the domain contains an element `y`
such that `y == x`.

!!! note
    **Incompatible element types.**
    In principle, the function `in(x, domain)` should not throw an exception even
    if the types seem mathematically nonsensical. In that case, the correct return
    value is `false`. This mimicks the behaviour of `in` for finite sets in Julia:
    ```julia
    julia> in(rand(3,3), 1:3)
    false
    ```
    Indeed, a `3x3` matrix is not equal to any of the numbers `1`, `2` or `3`.


## The `domaineltype` function

The defining mathematical condition of a continuous set might be satisfied by
variables of different types. Still, the interface defines the `domaineltype`
of a domain. It is a valid type for elements of the set.

Functions that generate elements of the domain should generate elements
of that type. As a consequence, for finite sets such as an `AbstractArray` or `AbstractSet`, the `domaineltype` agrees with the `eltype` of that set. For
intervals on the real line, the `domaineltype` might be `Float64`. When there is
no clear candidate the `domaineltype` might simply be `Any`.


## Minimal formal interface

The domain interface is formally summarised in the following table:

| Required methods           |                             | Brief description                        |
|:-------------------------- |:--------------------------- |:---------------------------------------- |
| `in(x, d)`                 |                             | Returns `true` when `x` is an element of the domain, `false` otherwise |
| `DomainStyle(d)`           |                             | Returns `IsDomain()` if `d` implements this interface |
| **Optional methods**       | **Default definition**      | **Brief description**                 |
| `domaineltype(d)`          | `eltype(d)`                 | Returns a valid type for elements of the domain |

Several extensions of this minimal interface are defined in the `DomainSets` package.


## The Domain supertype and DomainStyle trait

Domains in this package inherit from the abstract type `Domain{T}`. It is the
supertype of continuous sets with `domaineltype` equal to `T`.

The package also defines the trait `DomainStyle`. Any type can declare to
implement the domain interface by defining
```julia
DomainSets.DomainStyle(d::MyDomain) = IsDomain()
```
Objects of type `Number`, `AbstractArray` and `AbstractSet` are declared to be
domains in this package.
