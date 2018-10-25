# DomainSets.jl

[![Build Status](https://travis-ci.org/JuliaApproximation/DomainSets.jl.svg?branch=master)](https://travis-ci.org/JuliaApproximation/DomainSets.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaApproximation/DomainSets.jl/badge.svg)](https://coveralls.io/github/JuliaApproximation/DomainSets.jl)


DomainSets.jl is a package designed to represent simple infinite sets, that
can be used to encode domains of functions. For example, the domain of the
function `log(x::Float64)` is the infinite open interval, which is represented
by the type `HalfLine{Float64}()`.

## Examples

### Intervals

DomainSets.jl uses IntervalSets.jl for Closed and open intervals.

### Rectangles

Rectangles can be constructed as a product of intervals, where the elements of the domain
are `SVector{2}`:

```julia
julia> using DomainSets, StaticArrays; import DomainSets: ×

julia> SVector(1,2) in (-1..1) × (0..3)
true
```

### Circles and Spheres

A `UnitHyperSphere{N,T}`  contains `x::SVector{N,T}` if `norm(x) == one(T)`. `UnitCircle` and `UnitSphere` are two important cases:
```julia
julia> SVector(1,0) in UnitCircle()
true

julia> SVector(1,0,0) in UnitSphere()
true
```

### Disks and Balls

A `UnitHyperBall{N,T}`  contains `x::SVector{N,T}` if `norm(x) ≤ one(T)`. `UnitDisk` and `UnitHyperBall` are two important cases:
```julia
julia> SVector(0.1,0.2) in UnitDisk()
true

julia> SVector(0.1,0.2,0.3) in UnitBall()
true
```


### Union, intersection, and setdiff of domains

Domains can be unioned and intersected together:
```julia
julia> d = UnitCircle() ∪ 2UnitCircle();

julia> in.([SVector(1,0),SVector(0,2), SVector(1.5,1.5)], Ref(d))
3-element BitArray{1}:
  true
  true
 false

julia> d = UnitCircle() ∩ (2UnitCircle() + SVector(1.0,0.0))
the intersection of 2 domains:
	1.	: UnitHyperSphere{2,Float64}()
	2.	: A mapped domain based on UnitHyperSphere{2,Float64}()

julia> SVector(1,0) in d
false

julia> SVector(-1,0) in d
true
```

### The domain interface

A domain is any type that implements the functions `eltype` and `in`. If
`d` is an instance of a type that implements the domain interface, then
the domain consists of all `x` that is an `eltype(d)` such that `x in d`
returns true.

Domains often represent continuous mathematical domains, for example, a domain
`d`  representing the interval `[0,1]` would have `eltype(d) == Int` but still
have `0.2 in d` return true.

### The `Domain` type

DomainSets.jl contains an abstract type `Domain{T}`. All subtypes of `Domain{T}`
must implement the domain interface, and in addition support `convert(Domain{T}, d)`.
