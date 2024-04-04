# Euclidean geometry

Domains appearing in geometry are an important special case. A small number of
primitive domains are implemented.

DomainSets favours the definition of canonical types without state, which can
be implemented as singleton types without fields. It is made easy to map the
canonical domain to a more specific one, often using
an affine map. See: [Canonical domains](@ref), and in particular the functions [`mapfrom_canonical`](@ref) and [`mapto_canonical`](@ref).

For example, each ball, with any given radius and any center,
can be obtained as the translation and scaling of a unit ball with radius one
centered at the origin.


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



## Balls and spheres

The type hierarchy is as follows:
```
abstract Ball
|-> abstract UnitBall: radius is 1
    |-> StaticUnitBall: dimension is part of type
    |-> DynamicUnitBall: dimension is specified by int field
|-> GenericBall: stores center and radius. Here, the dimension can be
    obtained from the center, so there is only one type.
```


```@docs; canonical=false
Ball()
Sphere()
```

## Simplices
