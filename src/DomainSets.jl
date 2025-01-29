module DomainSets

using CompositeTypes, CompositeTypes.Display, CompositeTypes.Indexing
using IntervalSets
using LinearAlgebra
using Random
using StaticArrays
using FunctionMaps

################################
## Exhaustive list of imports
################################

import Base:
    ==, isapprox,
    # Set operations
    setdiff, in, isempty, issubset, intersect, union, &, \,
    # Arrays
    eltype, hash,
    # Types, promotions and conversions
    convert, promote,
    # Display
    show

import FunctionMaps:
	convert_eltype,
	convert_prectype,
	convert_numtype,
	prectype,
	numtype,
	factors,
	isrealtype,
	tointernalpoint,
    toexternalpoint,
    compatibleproductdims


import CompositeTypes: component, components

import IntervalSets: (..), endpoints, Domain, AbstractInterval, TypedEndpointsInterval,
                        leftendpoint, rightendpoint, isleftopen, isleftclosed,
                        isrightopen, isrightclosed, isopenset, isclosedset,
                        infimum, supremum
export ..

using FunctionMaps:
	CanonicalType,
	Equal,
	AbstractAffineMap,
	GenericLinearMap,
	GenericAffineMap,
	interval_map,
	UnitCircleMap,
	AngleMap,
	UnitDiskMap,
	StaticTypes,
	hashrec,
	euclideandimension,
	convert_eltype,
	promotable_eltypes,
	promote_prectype,
	promote_numtype,
	convert_fromcartesian,
	convert_tocartesian,
	nfactors,
	factor

################################
## Exhaustive list of exports
################################

## Utils

# from util/common.jl
export prectype, numtype,
    convert_numtype, promote_numtype,
    convert_prectype, promote_prectype,
    iscomposite, component, components, ncomponents

## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, VectorDomain,
    dimension,
    approx_in,
    isopenset, isclosedset, iscompact,
    boundary, ∂,
    interior, closure,
    isrealdomain,
    choice

# from generic/geometry.jl
export boundingbox,
    volume,
    normal, tangents,
    distance_to,
    infimum, supremum

# from generic/canonical.jl
export canonicaldomain, hascanonicaldomain,
    mapto, mapto_canonical, mapfrom_canonical,
    equaldomain, hasequaldomain,
    mapfrom_equaldomain, mapto_equaldomain,
    parameterdomain, parameterization, hasparameterization,
    mapfrom_parameterdomain, mapto_parameterdomain,
    isequaldomain

# from generic/lazy.jl
export superdomain

# from generic/productdomain.jl
export ProductDomain, productdomain, center,
    cartesianproduct,
    VcatDomain, VectorProductDomain, TupleProductDomain,
    factors

# from generic/mapped.jl
export MappedDomain,
    map_domain,
    mapped_domain,
    forward_map,
    inverse_map

# from generic/setoperations.jl
export UnionDomain, uniondomain,
    IntersectDomain, intersectdomain,
    SetdiffDomain, setdiffdomain,
    issubset_domain

# from applications/rotations.jl
export rotate


## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace, TypeDomain,
    emptyspace, fullspace, typedomain
# from domains/numbers.jl
export Integers, RealNumbers, Rationals, ComplexNumbers,
    ℕ, ℤ, ℚ, ℝ, ℂ, ℝ1, ℝ2, ℝ3, ℝ4
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval,
    OpenInterval, ClosedInterval,
    leftendpoint, rightendpoint, isleftopen, isrightopen,
    HalfLine, NegativeHalfLine, RealLine
# from domains/simplex.jl
export UnitSimplex,
    StaticUnitSimplex, DynamicUnitSimplex,
    EuclideanUnitSimplex, VectorUnitSimplex,
    corners
# from domains/point.jl
export Point
# from domains/ball.jl
export Ball, UnitBall,
    center, radius,
    StaticUnitBall, DynamicUnitBall,
    EuclideanUnitBall, VectorUnitBall,
    Sphere, UnitSphere,
    StaticUnitSphere, DynamicUnitSphere,
    VectorUnitSphere, EuclideanUnitSphere,
    Disk, UnitDisk, VectorUnitDisk,
    UnitCircle, VectorUnitCircle,
    ellipse, ellipse_shape, cylinder
# from domains/complex.jl
export ComplexUnitCircle, ComplexUnitDisk
# from domains/cube.jl
export UnitCube,
    StaticUnitCube, DynamicUnitCube,
    EuclideanUnitCube, VectorUnitCube,
    UnitSquare, UnitCube,
    Rectangle
# from domain/levelset.jl
export LevelSet, ZeroSet,
    SublevelSet, SubzeroSet,
    SuperlevelSet, SuperzeroSet,
    pseudolevel
# from domain/indicator.jl
export IndicatorFunction


include("util/common.jl")

include("generic/interface.jl")
include("generic/domain.jl")
include("generic/geometry.jl")
include("generic/canonical.jl")
include("generic/lazy.jl")
include("generic/productdomain.jl")
include("generic/setoperations.jl")
include("generic/mapped.jl")
include("generic/broadcast.jl")

include("domains/trivial.jl")
include("domains/numbers.jl")
include("domains/levelset.jl")
include("domains/point.jl")
include("domains/interval.jl")
include("domains/simplex.jl")
include("domains/ball.jl")
include("domains/complex.jl")
include("domains/cube.jl")
include("domains/indicator.jl")
include("domains/boundingbox.jl")

include("generic/generator.jl")

include("applications/coordinates.jl")
include("applications/random.jl")
include("applications/rotation.jl")

include("deprecated.jl")

end # module
