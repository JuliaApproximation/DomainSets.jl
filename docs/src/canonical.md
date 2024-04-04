# Canonical domains

The package favours the implementation of a small number of canonical domains,
in terms of which other domains can be defined. The mechanism to associate domains with
canonical domains, with various meanings of the word "canonical", turns out to be both flexible and productive.

This part of the documentation is oriented towards library developers, rather
than users. As a user, apart from learning applicable syntax, you'll also read
the "why" of it.

## A representative of a class

We will use a ball as a running example.
```julia
julia> using DomainSets, StaticArrays

julia> d = Ball(2.0, SA[0.3, 0.5])
Ball(2.0, [0.3, 0.5])
```
The type of `d` is `GenericBall`, which stores the radius and the center of the
ball. This is a common thing to do: a ball is defined mathematically by its
center and its radius.

The canonical domain is much simpler. In this case it is a unit ball or
`UnitDisk`, the 2D unit ball with radius 1. This domain can be represented without state, as a singleton type without fields. It is easy to reason about and simple to implement functionality for. In turn, much of the functionality for more general balls can be implemented in terms of the map from the unit disk. Better yet, that kind of functionality may itself be generic, i.e., the implementation of the generalized
functionality is not specific to balls. For example, the boundary of a domain is the boundary of the canonical domain, mapped by the same map that takes points of the canonical domain to points of the domain.

The map from the unit ball to a more general one is affine.
```julia
julia> canonicaldomain(d)
UnitDisk()

julia> m = mapfrom_canonical(d)
x -> 2.0 * x + v

v = 2-element SVector{2, Float64} with indices SOneTo(2):
 0.3
 0.5
```

The boundary of a unit ball is a unit sphere.
```julia
julia> boundary(canonicaldomain(d))
UnitCircle()

julia> boundary(d)
Sphere(2.0, [0.3, 0.5])
```
The code does not "know" that the boundary of a generic ball is also a generic
sphere. It simply applies the map to the unit circle:
```julia
julia> map_domain(m, boundary(canonicaldomain(d)))
Sphere(2.0, [0.3, 0.5])
```

Some canonical domains are [`UnitInterval`](@ref), [`UnitBall`](@ref), [`UnitSphere`](@ref) and [`UnitSimplex`](@ref).

##### Relevant functions

```@docs; canonical=false
canonicaldomain
hascanonicaldomain
mapfrom_canonical
mapto_canonical
```

## Equal domains

Different types of canonical domains can be obtained by using a different
`CanonicalType` argument to `canonicaldomain`. The first one we consider is
`Equal()`.

Domains can be the same even if they have different types. This is hard to capture in general without a symbolic engine. The system of canonical domains provides a partial solution by associating domains
with a simpler "equal" domain, whenever that is possible.

The unit ball in one dimension with scalar elements is the same as the interval $[-1,1]$.
```julia
julia> d = UnitBall{Float64}()
UnitBall{Float64}()

julia> canonicaldomain(DomainSets.Equal(), d)
-1.0 .. 1.0 (Chebyshev)

julia> equaldomain(d)
-1.0 .. 1.0 (Chebyshev)

julia> UnitBall{Float64}() == -1..1
true
```
The `equaldomain(d)` function is a convenience shorthand for `canonicaldomain(Equal(), d)`.

Similarly, the one-dimensional unit simplex with scalar element type is really
just the unit interval.
```julia
julia> equaldomain(UnitSimplex{Float64}())
0.0 .. 1.0 (Unit)
```


##### Relevant functions

```@docs; canonical=false
equaldomain
hasequaldomain
```

## Isomorphic domains

Two domains can be isomorphic but not equal. One example results from the
isomorphism between numbers and vectors of length 1, namely $x ↔ [x]$.

The unit ball with length 1 vector elements is isomorphic to the unit ball
with scalar element type:
```julia
julia> eltype(UnitBall(Val(1)))
SVector{1, Float64} (alias for SArray{Tuple{1}, Float64, 1, 1})

julia> canonicaldomain(DomainSets.Isomorphic(), UnitBall(Val(1)))
UnitBall{Float64}()

julia> mapfrom_canonical(DomainSets.Isomorphic(), UnitBall(Val(1)))
x : Float64 -> x : SVector{1, Float64}
```

Another common example is the identification of a region of the complex plane
with a domain in the 2D Euclidean plane, $[a;b] ↔ a+bi$. One such implementation
in the package is the complex unit circle:
```julia
julia> canonicaldomain(DomainSets.Isomorphic(), ComplexUnitCircle())
UnitCircle()

julia> mapto_canonical(DomainSets.Isomorphic(), ComplexUnitCircle())
x : ComplexF64 -> x : StaticArraysCore.SVector{2, Float64}
```

##### Relevant functions

```@docs; canonical=false
DomainSets.Isomorphic
```

The isomorphisms in the two examples above are defined in [FunctionMaps.jl](@ref).

```@docs; canonical=false
FunctionMaps.NumberToVector
FunctionMaps.VectorToNumber
FunctionMaps.ComplexToVector
FunctionMaps.VectorToComplex
```


## Parametric domains

A domain can be parameterised if it is easy to map from a parameter domain to
the domain at hand, but not the other way around. The map might not be invertible,
or might be difficult to invert (numerically or analytically).

We associate a circle with a radius and center with the unit circle. The unit
circle is the canonical domain. However, in turn, the unit circle can be
parameterised from the unit interval. Consequently, by chaining maps, any other
circle can be parameterised from the unit interval as well.

```julia
julia> parameterdomain(UnitCircle())
0.0 .. 1.0 (Unit)

julia> canonicaldomain(DomainSets.Parameterization(), UnitCircle())
0.0 .. 1.0 (Unit)

julia> mapfrom_parameterdomain(UnitCircle())
DomainSets.FunctionMaps.UnitCircleMap{Float64}()

julia> using StaticArrays

julia> canonicaldomain(DomainSets.Parameterization(), Sphere(2.0, SA[1.0, 5.0]))
0.0 .. 1.0 (Unit)

julia> mapfrom_parameterdomain(Sphere(2.0, SA[1.0, 5.0]))
(x -> 2.0 * x + v) ∘ UnitCircleMap{Float64}()

v = 2-element Vector{Float64}:
 1.0
 5.0
```
Note that `parameterdomain` is a convenience function similar to `equaldomain`.


Another common case is a line segment in 2D, which can be described in terms of an
interval and an affine map from ``ℝ`` to ``ℝ^2``.
```julia
julia> using StaticArrays

julia> m = AffineMap(SA[2.0, 3.0], SA[1.0, 2.0]);

julia> d = DomainSets.parametric_domain(m, UnitInterval());

julia> parameterdomain(d)
0.0 .. 1.0 (Unit)

julia> mapfrom_parameterdomain(d)
x -> v₁ * x + v₂

v₁ = 2-element SVector{2, Float64} with indices SOneTo(2):
 2.0
 3.0
v₂ = 2-element SVector{2, Float64} with indices SOneTo(2):
 1.0
 2.0
```

##### Relevant functions

```@docs; canonical=false
parameterdomain
mapfrom_parameterdomain
FunctionMaps.UnitCircleMap
```


## Simplifying a domain

The mechanisms involving various canonical types defined above allow one to
generically simplify a domain. The first possible simplification of a domain
is its [`equaldomain`](@ref).

However, even if the domain does not have an equivalent domain, its canonical
domain might have one. Or its parameter domain, if it exists. The function
[`DomainSets.simplify`](@ref) goes one step further than [`equaldomain`](@ref) and
also checks the latter property.

Since a scalar unit ball is an interval, any affine map of that ball also
corresponds to a similarly mapped interval.
```julia
julia> 2*UnitBall{Float64}()
Ball(2.0, 0.0)

julia> DomainSets.simplify(UnitBall{Float64}())
-1.0 .. 1.0 (Chebyshev)

julia> equaldomain(2*UnitBall{Float64}())
Ball(2.0, 0.0)

julia> DomainSets.simplify(2*UnitBall{Float64}())
-2.0 .. 2.0
```

The simplification of a domain is used in `isequaldomain` to simplify the
arguments before verifying equality. The mechanism above results in the following
positive identification of equality between domains with very different types:
```julia
julia> isequaldomain(2*UnitBall{Float64}() .+ 1, -1..3)
true
```

##### Relevant functions

```@docs; canonical=false
DomainSets.simplify
DomainSets.simplifies
```
