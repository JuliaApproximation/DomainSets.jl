module DomainSets

# We use static vectors internally

using StaticArrays
using Base

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
    convert, widen, promote_rule,
    # Display
    show,
    # Various
    Bool

# IntervalSets
import IntervalSets: (..), endpoints, Domain, AbstractInterval, TypedEndpointsInterval,
                        leftendpoint, rightendpoint, isleftopen, isleftclosed,
                        isrightopen, isrightclosed, isopen, isclosed,
                        infimum, supremum
export ..


################################
## Exhaustive list of exports
################################

## Utils

# from util/tensorproducts.jl
export cartesianproduct


## Spaces

# from spaces/space.jl
export spaceof, superspaceof, spacetype, origin, superspace, issubspace, subeltype
export isomorphic, ≅, embedded, ↪


## Maps

# from maps/maps.jl
export AbstractMap, applymap, jacobian, linearize
export domaintype, codomaintype
export left_inverse, right_inverse
export apply_inverse, apply_left_inverse, apply_right_inverse
export image
# from maps/composite_map.jl
export CompositeMap
export ∘
# from maps/productmap.jl
export ProductMap
export tensorproduct
# from maps/basic_maps.jl
export IdentityMap, ConstantMap
# from maps/embedding_map.jl
export embedding_map, restriction_map, isomorphism_map
# from maps/affine_maps.jl
export AffineMap, Translation, LinearMap
export islinear
export linear_map, interval_map, scaling_map, rotation_map, translation_map
# from maps/coordinates.jl
export CartToPolarMap, PolarToCartMap


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain
export dimension
export approx_in
export isclosed, iscompact
export boundary, ∂
export point_in_domain

# from generic/derived_domain.jl
export DerivedDomain
export superdomain

# from generic/productdomain.jl
export ProductDomain

# from generic/setoperations.jl
export UnionDomain, IntersectionDomain, DifferenceDomain
export TranslatedDomain

# from generic/mapped_domain.jl
export source, target


# from generic/arithmetics.jl
export rotate

export infimum, supremum


## Specific domains

# from domains/simple.jl
export EmptySpace, FullSpace
# from domains/trivial.jl
export ℤ, ℚ, ℝ, ℂ, ℝ1, ℝ2, ℝ3, ℝ4
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval,
    OpenInterval, ClosedInterval,
    leftendpoint, rightendpoint, isleftopen, isrightopen,
    cardinality,
    HalfLine, NegativeHalfLine
# from domains/simple.jl
export Cube, Simplex, UnitSimplex, Point
# from domains/circle.jl
export UnitCircle, UnitSphere, UnitHyperSphere,
    UnitDisk, UnitBall, UnitHyperBall,
    EuclideanUnitBall, VectorUnitBall,
    EuclideanUnitSphere, VectorUnitSphere,
    ellipse, ellipse_shape, cylinder,
    parameterization, gradient

include("util/common.jl")
include("util/products.jl")

include("spaces/space.jl")
include("spaces/space_promotions.jl")
include("spaces/basic_spaces.jl")
include("spaces/productspace.jl")

include("maps/maps.jl")
include("maps/composite_map.jl")
include("maps/productmap.jl")
include("maps/basic_maps.jl")
include("maps/embedding_map.jl")
include("maps/affine_map.jl")
include("maps/coordinates.jl")

include("generic/domain.jl")
include("generic/lazy.jl")
include("generic/derived_domain.jl")
include("generic/productdomain.jl")
include("generic/setoperations.jl")
include("generic/mapped_domain.jl")
include("generic/promotion.jl")
include("generic/arithmetics.jl")

include("domains/trivial.jl")
include("domains/interval.jl")
include("domains/simple.jl")
include("domains/circle.jl")

end # module
