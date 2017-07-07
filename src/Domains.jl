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
import Base: intersect, union, setdiff, in

# Arrays
import Base: size, length, ndims, getindex, eltype, ndims
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
# from util/products.jl
export flatten, cartesianproduct


## Spaces

# from spaces/space.jl
export GeometricSpace, AnySpace
export spaceof, spacetype, origin, superspace, issubspace, subeltype
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
export domaintype, rangetype
# from maps/composite_map.jl
export CompositeMap
export ∘
# from maps/productmap.jl
export ProductMap
# from maps/basic_maps.jl
export IdentityMap
# from maps/embedding_map.jl
export embedding_map, restriction_map, isomorphism_map
# from maps/affine_maps.jl
export AffineMap, Translation, LinearMap
export islinear, matrix, vector
export linear_map, interval_map, scaling_map
# from maps/coordinates.jl
export CartToPolarMap, PolarToCartMap


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, Domain1d, Domain2d, Domain3d, Domain4d
export indomain
export isclosed, isopen, iscompact

# from generic/derived_domain.jl
export DerivedDomain
export superdomain

# from generic/productdomain.jl
export ProductDomain

# from generic/setoperations.jl
export UnionDomain, IntersectionDomain, DifferenceDomain
export TranslatedDomain

# from generic/mapped_domain.jl
export MappedDomain
export mapping, map_domain

# from generic/arithmetics.jl
export rotate


## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace, AnyEmptySpace, AnyFullSpace
export euclideanspace, emptyspace, fullspace
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval
export real_line, halfline, negative_halfline
export interval, leftendpoint, rightendpoint
export similar_interval
# from domains/simple.jl
export UnitBall, Disk, Ball, Cube, Simplex, UnitSimplex, UnitSphere
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
