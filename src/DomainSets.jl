module DomainSets

using StaticArrays
using LinearAlgebra, Statistics
import LinearAlgebra: cross, ×, pinv

using IntervalSets

################################
## Exhaustive list of imports
################################

# Generated functions
import Base: @ncall

# Operator symbols
import Base: *, +, -, /, \, ^,
    |, &,
    ∪, ∩,
    ==, isapprox,
    ∘,
    # Set operations
    intersect, union, setdiff, in, isempty, minimum, maximum,
    issubset,
    # Arrays
    size, length, ndims, getindex, eltype, ndims, hash,
    inv, isreal, zero,
    # Types, promotions and conversions
    convert, widen, promote,
    # Display
    show,
    # Various
    Bool

# IntervalSets
import IntervalSets: (..), endpoints, Domain, AbstractInterval, TypedEndpointsInterval,
                        leftendpoint, rightendpoint, isleftopen, isleftclosed,
                        isrightopen, isrightclosed, isopenset, isclosedset,
                        infimum, supremum
export ..


################################
## Exhaustive list of exports
################################

## Utils

# from util/common.jl
export prectype, numtype,
    convert_numtype, promote_numtype,
    convert_prectype, promote_prectype,
    iscomposite, elements, element, numelements

## Maps

# from maps/map.jl
export AbstractMap, Map, TypedMap,
    applymap,
    domaintype, codomaintype,
    inverse, leftinverse, rightinverse,
    jacobian, jacdet
# from maps/lazy.jl
export WrappedMap
# from maps/composite.jl
export Composition, ∘
# from maps/product.jl
export ProductMap, productmap
# from maps/basic.jl
export IdentityMap,
    StaticIdentityMap, VectorIdentityMap,
    ZeroMap, UnityMap, ConstantMap,
    isconstant, constant
# from maps/affine.jl
export AffineMap, Translation, LinearMap,
    matrix, vector,
    islinear, isaffine


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, VectorDomain,
    dimension,
    approx_in,
    isopenset, isclosedset, iscompact,
    boundary, ∂,
    interior, closure,
    volume,
    point_in_domain,
    canonicaldomain, tocanonical, fromcanonical,
    mapto,
    parameterdomain, parameterization

# from generic/lazy.jl
export DerivedDomain, superdomain, WrappedDomain

# from generic/productdomain.jl
export ProductDomain, productdomain,
    VcatDomain, VectorProductDomain, TupleProductDomain

# from generic/mapped.jl
export MappedDomain,
    map_domain,
    mapped_domain,
    forward_map,
    inverse_map

# from generic/setoperations.jl
export UnionDomain, uniondomain,
    IntersectDomain, intersectdomain,
    SetdiffDomain, setdiffdomain

# from applications/rotations.jl
export rotate

export infimum, supremum


## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace,
    emptyspace, fullspace,
    ℤ, ℚ, ℝ, ℂ, ℝ1, ℝ2, ℝ3, ℝ4
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval,
    OpenInterval, ClosedInterval,
    leftendpoint, rightendpoint, isleftopen, isrightopen,
    cardinality,
    HalfLine, NegativeHalfLine
# from domains/simplex.jl
export EuclideanUnitSimplex, VectorUnitSimplex, UnitSimplex,
    center, corners
# from domains/point.jl
export Point
# from domains/ball.jl
export UnitBall,
    StaticUnitBall, DynamicUnitBall,
    EuclideanUnitBall, VectorUnitBall,
    UnitSphere,
    StaticUnitSphere, DynamicUnitSphere,
    VectorUnitSphere, EuclideanUnitSphere,
    UnitDisk, VectorUnitDisk,
    UnitCircle, VectorUnitCircle,
    ComplexUnitCircle, ComplexUnitDisk,
    ellipse, ellipse_shape, cylinder,
    gradient
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

## Applications
# from applications/rotation.jl
export rotation_map,
    CartToPolarMap, PolarToCartMap

include("util/common.jl")
# include("util/products.jl")

include("maps/map.jl")
include("maps/lazy.jl")
include("maps/composite.jl")
include("maps/product.jl")
include("maps/isomorphism.jl")
include("maps/basic.jl")
include("maps/affine.jl")
include("maps/arithmetics.jl")

include("generic/domain.jl")
include("generic/canonical.jl")
include("generic/lazy.jl")
include("generic/productdomain.jl")
include("generic/setoperations.jl")
include("generic/mapped.jl")
include("generic/broadcast.jl")

include("domains/trivial.jl")
include("domains/levelset.jl")
include("domains/point.jl")
include("domains/interval.jl")
include("domains/simplex.jl")
include("domains/ball.jl")
include("domains/cube.jl")
include("domains/indicator.jl")

include("applications/rotation.jl")

end # module
