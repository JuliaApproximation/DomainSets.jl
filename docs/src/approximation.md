# Approximation domains

DomainSets has its roots in the JuliaApproximation community. It is used for
example in the packages [ApproxFun](https://github.com/JuliaApproximation/ApproxFun.jl)
and [BasisFunctions](https://github.com/JuliaApproximation/BasisFunctions.jl).

## Approximation intervals

There are a number of stateless intervals (i.e. singleton types without data) to
represent standard intervals used in approximation theory. The endpoints of
these intervals are fixed and hence, at least in principle, known to the compiler.

```@docs; canonical=false
UnitInterval
ChebyshevInterval
HalfLine
NegativeHalfLine
```

The types are simple enough that some common set arithmetic operations retain
their special structure:
```julia
julia> HalfLine() ∩ ChebyshevInterval()
0.0 .. 1.0 (Unit)

julia> NegativeHalfLine() \ UnitInterval()
-Inf .. 0.0 (open) (NegativeHalfLine)
```

## The unit cube and other cubes

The product domain associated with any of the fixed intervals remains stateless
if the dimension is fixed as well:
```julia
julia> UnitInterval()^2
UnitSquare()

julia> ChebyshevInterval()^3
(-1.0 .. 1.0 (Chebyshev)) × (-1.0 .. 1.0 (Chebyshev)) × (-1.0 .. 1.0 (Chebyshev))
```

Some special cases are named.
```@docs; canonical=false
UnitSquare
UnitCube
```
