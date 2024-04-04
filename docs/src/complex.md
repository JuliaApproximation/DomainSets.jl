# Complex analysis

## Regions of the complex plane

A number of domains relate to the complex plane:
- [ComplexUnitCircle](@ref)
- [ComplexUnitDisk](@ref)

For a more specific package dealing with complex regions, see [ComplexRegions.jl](https://github.com/complexvariables/ComplexRegions.jl).

## Identification with the plane

Any domain in $ℝ^2$ can be identified with a region of the complex plane. The
relevant maps to use are [`FunctionMaps.VectorToComplex`](@ref) and [`FunctionMaps.ComplexToVector`](@ref).

```@docs; canonical=false
FunctionMaps.VectorToComplex
FunctionMaps.ComplexToVector
```
