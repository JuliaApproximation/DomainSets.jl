# Set operations

DomainSets implements a number of standard set operations. In some cases the
operation can be performed analytically and a concrete domain is returned. More
often than not, however, the operation results in a lazy structure.

Functions like `uniondomain` and `intersectdomain` are not restricted to
arguments of type `Domain`. Thus, one can construct the lazy union of any two
objects, regardless of their types.

A distinction is made for set operations between functions in lowercase, such as `uniondomain`, and capitalized functions, in this case `UnionDomain`. The former
attempts to simplify the arguments (see [Canonical domains](@ref)) and returns
a domain that is mathematically equivalent to the union of the two domains.
The latter is the constructor of a type and hence always returns an instance of `UnionDomain`.

!!! note
    For concrete subtypes of `Domain` one can use standard Julia operators and
    functions for common set operations. The relevant symbols and functions are
    `∪` (union), `∩` (intersect) and `∖` (setdiff). For example:
    ```julia
    julia> UnitBall(3) ∩ Point([0,0,0])
    Point([0.0, 0.0, 0.0])
    ```

!!! warning
    Set operations in DomainSets are often **not type-stable**. For example,
    by convention `uniondomain(d1,d2)` returns a domain of the simplest type
    that is mathematically equivalent to the union of `d1` and `d2`. The union
    of two overlapping intervals is a single interval, but the union of
    non-overlapping intervals is a `UnionDomain`:
    ```julia
    julia> uniondomain(1..3, 2..4)
    1 .. 4
    julia> uniondomain(1..3, 5..7)
    (1 .. 3) ∪ (5 .. 7)
    ```
    When type-safety is important, use the corresponding constructor:
    ```julia
    julia> UnionDomain(1..3, 2..4)
    (1 .. 3) ∪ (2 .. 4)
    ```


## Product domains

Product domains are created most easily by invoking the `ProductDomain`
constructor. The constructor of the associated abstract type `ProductDomain`
returns an instance of a suitable concrete subtype. In many cases a user need
not be aware of which type is being returned by `ProductDomain`, as it always
behaves like the requested domain.


```@docs; canonical=false
ProductDomain()
```

A number of concrete product domain types are implemented. They differ in what
the `eltype` of the product domain is. In the most generic case `T` is a tuple,
with each element representing the element type of the corresponding factor. A
`VcatDomain` is a special case for product domains of Euclidean type, i.e.,
one in which each factor has a scalar or statically-sized vector as element type. Finally, a `VectorProductDomain` has a `Vector` element type whose dimension is determined by the number of factors. The `ProductDomain` aims to return the most specialized type of domain, but the individual constructors may be invoked to ensure a
specific one.


```julia
julia> ProductDomain(2..4.5, 3.0..5.0)
(2.0 .. 4.5) × (3.0 .. 5.0)

julia> eltype(ProductDomain(2..4.5, 3.0..5.0))
SVector{2, Float64}

julia> [2.4, 4] ∈ ProductDomain(2..4.5, 3.0..5.0)
true

julia> TupleProductDomain(2..4.5, 3.0..5.0)
(2.0 .. 4.5) × (3.0 .. 5.0)

julia> eltype(TupleProductDomain(2..4.5, 3.0..5.0))
Tuple{Float64, Float64}

julia> (2.4, 4) ∈ TupleProductDomain(2..4.5, 3.0..5.0)
true

julia> ProductDomain([ i..i+1.0 for i in 1:10])
(1.0 .. 2.0) × (2.0 .. 3.0) × (3.0 .. 4.0) × ... × (10.0 .. 11.0)

julia> eltype(ProductDomain([ i..i+1.0 for i in 1:10]))
Vector{Float64}

julia> 1:10 ∈ ProductDomain([ i..i+1.0 for i in 1:10])
true
```

It is noteworthy that a `VcatDomain` can cope with the concatenation of scalars
and vectors. In the example below, a cylinder is represented as the product of
a two-dimensional disk with a one-dimensional interval.
```julia
julia> using DomainSets: ×

julia> cylinder = UnitDisk() × UnitInterval()
UnitDisk() × (0.0 .. 1.0 (Unit))

julia> eltype(cylinder)
SVector{3, Float64}

julia> [0.4,0.2,0.6] ∈ cylinder
true
```

A number of concrete product domains are implemented in `DomainSets`.
```@docs; canonical=false
VcatDomain
VectorProductDomain
TupleProductDomain
Rectangle
```

## Union of sets

The mathematical union of two sets is guaranteed by `uniondomain`, while a lazy
union is returned by `UnionDomain`.

For vectors and sets in Julia, `uniondomain` returns precisely what the standard
`union` operation would do, while `UnionDomain` returns a lazy construct:
```julia
julia> uniondomain(1:3, 10)
4-element Vector{Int64}:
  1
  2
  3
 10
julia> UnionDomain(1:3, 10)
1:3 ∪ 10

julia> 10 ∈ ans
true
```

```@docs; canonical=false
uniondomain
UnionDomain
```

## Set intersection

The mathematical intersection of two sets is guaranteed by `intersectdomain`,
while a lazy intersection is returned by `IntersectDomain`.

```@docs; canonical=false
intersectdomain
IntersectDomain
```


## Set difference

The mathematical difference of two sets is guaranteed by `setdiffdomain`,
while a lazy difference is returned by `SetdiffDomain`.

```@docs; canonical=false
setdiffdomain
SetdiffDomain
```
