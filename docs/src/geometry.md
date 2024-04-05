# Euclidean geometry

Domains appearing in geometry are an important special case. A small number of
primitive domains are implemented.

Several of these may be open or closed.

DomainSets favours the definition of canonical types without state, which can
be implemented as singleton types without fields or data. It is made easy to map the
canonical domain to a more specific one, often using
an affine map. See: [Canonical domains](@ref).

## A note on Euclidean dimension

The Euclidean dimension of a domain can be implicit in the element type, if it
is scalar or a static vector. One could have `T=Float64` or
`T = SVector{3,Float64}` for 3D domains. The dimension can also be variable, if the element type is a `Vector{T}`. In the latter case the dimension is typically
explicitly stored in a field of the type. In the former case the dimension is implicit in the type and the type itself can be stateless.

!!! note
    `DomainSets` imposes no preference and enables both representations. Among other
    things, this leads to having more types, sometimes with funny or unexpected names.
    For example, `StaticUnitBall{T}` is a singleton type,
    whose dimension is fixed by the element type `T`. The word "static" in the name of the type refers to this aspect of dimensionality being known at compile-time.
    The alternative is `DynamicUnitBall{T}`. This is a domain with element type `Vector{T}`, whose dimension is explicitly stored in a field.


## Intervals

DomainSets uses [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl)
for the definition and manipulation of intervals. Most notably intervals can be
easily created using the ellipsis notation `a..b`.

`IntervalSets` supports open and closed intervals, finite and infinite intervals,
and intervals of non-numeric types such as dates and times. See the documentation
of the package for the full interface.

A number of intervals with fixed endpoints, known at compile-time, are added in
`DomainSets`, see [Approximation intervals](@ref).

## Rectangles and cubes

Rectangles and cubes are implemented as a `ProductDomain` (see [Product domains](@ref)).
They can be obtained by passing closed intervals to the `ProductDomain` constructor,
or by invoking the `Rectangle` constructor.

Two special cases are the `UnitSquare` and `UnitCube` types. These types have no state.

Some examples:
```julia
julia> ProductDomain(1.5..2.5, 3.5..6.0)
(1.5 .. 2.5) × (3.5 .. 6.0)

julia> [0.5,0.3] ∈ UnitSquare()
true

julia> [0.5,0.3,0.1] ∈ UnitCube()
true

julia> ProductDomain(UnitInterval(), UnitInterval())
UnitSquare()

julia> Rectangle([0.0,1.0,2.0], [1.0,4.0,6.0])
(0.0 .. 1.0) × (1.0 .. 4.0) × (2.0 .. 6.0)
```

##### Relevant functions

```@docs; canonical=false
Rectangle
UnitSquare
UnitCube
```

## Balls and spheres

The unit ball in a space is the volume consisting of all points `x` for which
`norm(x) < 1` or `norm(x) ≤ 1`, depending on whether the ball is open or closed.
The sphere is the boundary of the volume, i.e., the set of points for which
`norm(x) = 1`.

The canonical ball is a `UnitBall`, which can either have a dimension determined
by the element type `T` or, in case of `Vector` elements, a dimension that is
specified explicitly.

The constructor of the abstract type `Ball` returns a suitable concrete type depending on the arguments.

Types and syntax for spheres are analogous to those of balls.
```julia
julia> UnitBall()
UnitBall()

julia> Ball(3.0)
Ball(3.0, [0.0, 0.0, 0.0])

julia> Ball(3.0, [0.4, 2.0, 5.0])
Ball(3.0, [0.4, 2.0, 5.0])

julia> Ball{SVector{3,Float64}}(3, [0, 1, 2])
Ball(3.0, [0.0, 1.0, 2.0])

julia> eltype(ans)
SVector{3, Float64}
```

##### Relevant functions

```@docs; canonical=false
Ball()
Sphere()
UnitBall()
UnitSphere()
```

The unit disk and unit circle are special cases in 2D.
```@docs; canonical=false
UnitDisk
UnitCircle
```


## Simplices

The type hierarchy is similar to that of balls. The abstract supertype
`UnitSimplex` has a constructor which returns a suitable concrete subtype
depending on the arguments.

```@docs; canonical=false
UnitSimplex()
```
