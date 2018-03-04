# Domains.jl

module Domains

# We use static vectors internally
using StaticArrays


################################
## Exhaustive list of imports
################################

# Operator symbols
import Base: *, +, -, /, \, ^
import Base: |, &
import Base: ∪, ∩
import Base: ==
import Base: ∘
import Base: cross, ×

# Set operations
import Base: intersect, union, setdiff, in, isempty, minimum, maximum

# Arrays
import Base: size, length, ndims, getindex, eltype, ndims, hash
import Base: inv
import Base: isreal
import Base: zero
import Base: gradient

# Types, promotions and conversions
import Base: convert, widen

# Iteration protocol
import Base: start, next, done

# Display
import Base: show

# Various
import Base: isopen, Bool


################################
## Exhaustive list of exports
################################

## Utils

# from util/common.jl
export elements, element, nb_elements
export TypeFactory

# from util/tensorproducts.jl
export cartesianproduct


## Spaces

# from spaces/space.jl
export GeometricSpace, AnySpace
export spaceof, superspaceof, spacetype, origin, superspace, issubspace, subeltype
# from spaces/space_promotions.jl
export convert_space, restrict_space, demote, promote_space, promote_space_type
export isomorphic, ≅, embedded, ↪
# from spaces/basic_spaces.jl
export IntegerSpace, RationalSpace, RealSpace, ComplexSpace
export VectorSpace, EuclideanSpace, ArraySpace
export ℤ, ℚ, ℝ, ℂ, ℝ1, ℝ2, ℝ3, ℝ4
# from spaces/productspace.jl
export ProductSpace


## Maps

# from maps/maps.jl
export AbstractMap, applymap, apply_inverse, jacobian, linearize
export domaintype, codomaintype
export left_inverse, right_inverse
# from maps/composite_map.jl
export CompositeMap
export ∘
# from maps/productmap.jl
export ProductMap
export tensorproduct
# from maps/basic_maps.jl
export IdentityMap
# from maps/embedding_map.jl
export embedding_map, restriction_map, isomorphism_map
# from maps/affine_maps.jl
export AffineMap, Translation, LinearMap
export islinear, matrix, vector
export linear_map, interval_map, scaling_map, rotation_map, translation_map
# from maps/coordinates.jl
export CartToPolarMap, PolarToCartMap


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, Domain1d, Domain2d, Domain3d, Domain4d
export indomain, dimension
export approx_in, approx_indomain
export isclosed, isopen, iscompact
export boundary, ∂
export point_in_domain

# from generic/derived_domain.jl
export DerivedDomain
export superdomain

# from generic/productdomain.jl
export ProductDomain
export product_eltype

# from generic/setoperations.jl
export UnionDomain, IntersectionDomain, DifferenceDomain
export TranslatedDomain

# from generic/mapped_domain.jl
export MappedDomain, ForwardMappedDomain, InverseMappedDomain, BidirectionalMappedDomain
export source, target, forward_map, inverse_map
export map_domain, forwardmap_domain, inversemap_domain

# from generic/arithmetics.jl
export rotate

export infimum, supremum


## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace, AnyEmptySpace, AnyFullSpace
export euclideanspace, emptyspace, fullspace
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval, HalfOpenLeftInterval, HalfOpenRightInterval, OpenInterval, ClosedInterval
export real_line, halfline, negative_halfline, interval, open_interval,
    closed_interval, unitinterval
export leftendpoint, rightendpoint, open_left, open_right
export similar_interval
# from domains/simple.jl
export UnitBall, Disk, Ball, Cube, Simplex, UnitSimplex, UnitSphere, Point
export disk, ball, cube, simplex, cylinder, rectangle
# from domains/circle.jl
export Circle, Sphere
export circle, sphere
export parameterization, gradient

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
include("generic/derived_domain.jl")
include("generic/productdomain.jl")
include("generic/setoperations.jl")
include("generic/mapped_domain.jl")
include("generic/arithmetics.jl")

include("domains/trivial.jl")
include("domains/interval.jl")
include("domains/simple.jl")
include("domains/circle.jl")

end # module
