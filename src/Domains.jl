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

# Set operations
import Base: intersect, union, setdiff, in

# Arrays
import Base: size, length, ndims, getindex, eltype, ndims
import Base: inv
import Base: isreal
import Base: zero

# Types, promotions and conversions
import Base: convert, widen

# Iteration protocol
import Base: start, next, done

# Display
import Base: show


################################
## Exhaustive list of exports
################################

## Utils

# from util/common.jl
export elements, element, nb_elements
export TypeFactory
# from util/tensorproducts.jl
export flatten, tensorproduct, ⊗


## Spaces

# from spaces/space.jl
export GeometricSpace, AnySpace
export spaceof, spacetype, origin, superspace, issubspace, subeltype
# from spaces/space_promotions.jl
export convert_space, promote_space, promote_space_type
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
# from maps/affine_maps.jl
export AffineMap, Translation, LinearMap
export linear_map, interval_map, scaling_map
# from maps/productmap.jl
export ProductMap
# from maps/coordinates.jl
export CartToPolarMap, PolarToCartMap
# from maps/basic_maps.jl
export IdentityMap
# from maps/composite_map.jl
export CompositeMap


## Generic domains

# from generic/domain.jl
export Domain, EuclideanDomain, Domain1d, Domain2d, Domain3d, Domain4d
export indomain

# from generic/derived_domain.jl
export DerivedDomain
export superdomain

# from generic/productdomain.jl
export ProductDomain, tensorproduct, ⊗

# from generic/setoperations.jl
export UnionDomain, IntersectionDomain, DifferenceDomain
export TranslatedDomain

# from generic/mapped_domain.jl
export MappedDomain
export mapping, map_domain

# from generic/arithmetics.jl



## Specific domains

# from domains/trivial.jl
export EmptySpace, FullSpace, AnyEmptySpace, AnyFullSpace
export euclideanspace, emptyspace, fullspace
# from domains/interval.jl
export AbstractInterval, Interval, UnitInterval, ChebyshevInterval
export interval, leftendpoint, rightendpoint
# from domains/simple.jl
export UnitBall, Disk, Ball, Cube
export rectangle, cube, cylinder, randomcircles


include("util/common.jl")
include("util/tensorproducts.jl")

include("spaces/space.jl")
include("spaces/space_promotions.jl")
include("spaces/basic_spaces.jl")
include("spaces/productspace.jl")

include("maps/maps.jl")
include("maps/productmap.jl")
include("maps/composite_map.jl")
include("maps/affine_map.jl")
include("maps/basic_maps.jl")
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

end # module
