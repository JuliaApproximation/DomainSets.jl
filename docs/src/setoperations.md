# Set operations

DomainSets implements a number of standard set operations. In some cases the
operation can be performed analytically and a concrete domain is returned. More
often than not, however, the operation results in a lazy structure.

!!! warning
    Set operations in DomainSets are often not type-stable. By convention,
    `uniondomain(d1,d2)` returns a domain of the simplest type that is
    mathematically equivalent to the union of `d1` and `d2`. The union of two
    overlapping intervals is a single interval, but the union of non-overlapping
    intervals is a `UnionDomain`:
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

Functions like `uniondomain` are not restricted to arguments of type `Domain`.
Thus, one can construct the lazy union of any two objects, regardless of their
types.

!!! note
    For concrete subtypes of `Domain` one can use standard Julia operators and
    functions for common set operations. The relevant symbols and functions are
    `∪` (union), `∩` (intersect) and `∖` (setdiff). For example:
    ```julia
    julia> UnitBall(3) ∩ Point([0,0,0])
    Point([0.0, 0.0, 0.0])
    ```


## Product domains

See:  [`ProductDomain`](@ref)

## Union of sets

See:  [`UnionDomain`](@ref)

## Set intersection

## Set difference
