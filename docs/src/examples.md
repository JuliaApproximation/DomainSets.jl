# Examples

```@meta
DocTestSetup = :(using DomainSets, StaticArrays)
```

These are the same examples as those in the README of the package repository.

## Intervals from IntervalSets

DomainSets.jl uses [IntervalSets.jl](https://github.com/JuliaMath/IntervalSets.jl) for closed and open intervals. In addition, it defines a few standard intervals.

```jldoctest
julia> using DomainSets

julia> UnitInterval()
0.0 .. 1.0 (Unit)

julia> ChebyshevInterval()
-1.0 .. 1.0 (Chebyshev)

julia> HalfLine()
0.0 .. Inf (closed-open) (HalfLine)
```

## Rectangles

Rectangles can be constructed as a product of intervals, where the elements of the domain
are `SVector{2}`:

```jldoctest
julia> using DomainSets: ×

julia> (-1..1) × (0..3) × (4.0..5.0)
(-1.0 .. 1.0) × (0.0 .. 3.0) × (4.0 .. 5.0)

julia> [1,2] in (-1..1) × (0..3)
true

julia> UnitInterval()^3
UnitCube()
```

## Circles and Spheres

A `UnitSphere`  contains `x` if `norm(x) == 1`. The unit sphere is N-dimensional,
and its dimension is specified with the constructor. The element types are
`SVector{N,T}` when the dimension is specified as `Val(3)`, and they
are `Vector{T}` when the dimension is specified by an integer value instead:
```jldoctest
julia> using StaticArrays

julia> SA[0,0,1.0] in UnitSphere(Val(3))
true

julia> [0.0,1.0,0.0,0.0] in UnitSphere(4)
true
```
`UnitSphere` itself is an abstract type, hence the examples above return
concrete types `<:UnitSphere`. The intended element type can also be explicitly
specified with the `UnitSphere{T}` constructor:
```jldoctest
julia> typeof(UnitSphere{SVector{3,BigFloat}}())
EuclideanUnitSphere{3, BigFloat} (alias for StaticUnitSphere{SArray{Tuple{3}, BigFloat, 1, 3}})

julia> typeof(UnitSphere{Vector{Float32}}(6))
VectorUnitSphere{Float32} (alias for DynamicUnitSphere{Array{Float32, 1}})
```

Without arguments, `UnitSphere()` defaults to a 3D domain with `SVector{3,Float64}`
elements. Similarly, there is a special case `UnitCircle` in 2D:
```jldoctest
julia> SVector(1,0) in UnitCircle()
true
```



## Disks and Balls

A `UnitBall`  contains `x` if `norm(x) ≤ 1`. As with `UnitSphere`, the dimension
is specified via the constructor by type or by value:
```jldoctest
julia> SVector(0.1,0.2,0.3) in UnitBall(Val(3))
true

julia> [0.1,0.2,0.3,-0.1] in UnitBall(4)
true
```
By default `N=3`, but `UnitDisk` is a special case in 2D, and so are `ComplexUnitDisk` and `ComplexUnitCircle` in the complex plane:
```jldoctest
julia> SVector(0.1,0.2) in UnitDisk()
true

julia> 0.5+0.2im ∈ ComplexUnitDisk()
true
```

`UnitBall` itself is an abstract type, hence the examples above return
concrete types `<:UnitBall`. The types are similar to those associated with
`UnitSphere`. Like intervals, balls can also be open or closed:
```jldoctest
julia> EuclideanUnitBall{3,Float64,:open}()  # the 3-dimensional open unit ball
UnitBall()  (open)
```


## Cartesian products

The cartesian product of domains is constructed with the `ProductDomain` or
`ProductDomain{T}` constructor. This abstract constructor returns concrete types
best adapted to the arguments given.

If `T` is not given, `ProductDomain` makes a suitable choice based on the
arguments. If all arguments are Euclidean, i.e., their element types are numbers
or static vectors, then the product is a Euclidean domain as well:
```jldoctest
julia> ProductDomain(0..2, UnitCircle())
(0.0 .. 2.0) × UnitCircle()

julia> eltype(ans)
SVector{3, Float64} (alias for SArray{Tuple{3}, Float64, 1, 3})
```
The elements of the interval and the unit circle are flattened into a single
vector, much like the `vcat` function. The result is a `VcatDomain`.

If a `Vector` of domains is given, the element type is a `Vector` as well:
```jldoctest
julia> 1:5 in ProductDomain([0..i for i in 1:5])
true
```
In other cases, the points are concatenated into a tuple and membership is
evaluated element-wise:
```jldoctest
julia> ("a", 0.4) ∈ ProductDomain(["a","b"], 0..1)
true
```

Some arguments are recognized and return a more specialized product domain.
Examples are the unit box and more general hyperrectangles:
```jldoctest
julia> ProductDomain(UnitInterval(), UnitInterval())
UnitSquare()

julia> ProductDomain(0..2, 4..5, 6..7.0)
(0.0 .. 2.0) × (4.0 .. 5.0) × (6.0 .. 7.0)

julia> typeof(ans)
Rectangle{SVector{3, Float64}}
```


## Union, intersection, and setdiff of domains

Domains can be unioned and intersected together:
```jldoctest
julia> d = UnitCircle() ∪ 2UnitCircle();

julia> in.([SVector(1,0),SVector(0,2), SVector(1.5,1.5)], d)
3-element BitVector:
 1
 1
 0

julia> d = UnitCircle() ∩ (2UnitCircle() .+ SVector(1.0,0.0))
UnitCircle() ∩ Sphere(2.0, [1.0, 0.0])

julia> SVector(1,0) in d
false

julia> SVector(-1,0) in d
true
```


## Level sets

A domain can be defined by the level sets of a function. The domains of all
points `[x,y]` for which `x*y = 1` or `x*y >= 1` are represented as follows:
```jldoctest
julia> d = LevelSet{SVector{2,Float64}}(prod, 1.0)
level set f(x) = 1.0 with f = prod

julia> [0.5,2] ∈ d
true

julia> SuperlevelSet{SVector{2,Float64}}(prod, 1.0)
superlevel set f(x) >= 1.0 with f = prod
```
There is also `SublevelSet`, and there are the special cases `ZeroSet`,
`SubzeroSet` and `SuperzeroSet`.

## Indicator functions

A domain can be defined by an indicator function or a characteristic function.
This is a function `f(x)` which evaluates to true or false, depending on whether or
not the point `x` belongs to the domain.
```jldoctest
julia> d = IndicatorFunction{Float64}( t ->  cos(t) > 0)
indicator domain defined by function f = #2

julia> 0.5 ∈ d, 3.1 ∈ d
(true, false)
```
This enables generator syntax to define domains:
```jldoctest
julia> d = Domain(x>0 for x in -1..1)
indicator function bounded by: -1 .. 1

julia> 0.5 ∈ d, -0.5 ∈ d
(true, false)

julia> d = Domain( x*y > 0 for (x,y) in UnitDisk())
indicator function bounded by: UnitDisk()

julia> [0.2, 0.3] ∈ d, [0.2, -0.3] ∈ d
(true, false)

julia> d = Domain( x+y+z > 0 for (x,y,z) in ProductDomain(UnitDisk(), 0..1))
indicator function bounded by: UnitDisk() × (0.0 .. 1.0)

julia> [0.3,0.2,0.5] ∈ d
true
```
