# Domains.jl

module Domains

# We use static vectors internally
using StaticArrays


################################
## Exhaustive list of imports
################################

# Operator symbols
import Base: *, +, -, /, \, ^, |, &
import Base: ∪, ∩
import Base: ==

# Set operations
import Base: intersect, union, setdiff, in

# Arrays
import Base: size, length, ndims, getindex
import Base: inv

import Base: broadcast

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
# from util/tensorproducts.jl
export flatten, tensorproduct, ⊗
# from util/box.jl
export BBox, BBox1, BBox2, BBox3, BBox4
export ⊂


## Maps

# from maps/maps.jl
export AbstractMap, forward_map, inverse_map, jacobian, linearize
export is_linear, is_compatible
# from maps/affine_maps.jl
export AffineMap, translation, rotation, linear_map, interval_map, scaling_map
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
export Domain, Domain1d, Domain2d, Domain3d, Domain4d
export indomain, boundingbox
export left, right

# from generic/productdomain.jl
export ProductDomain, tensorproduct, ⊗

# from generic/arithmetics.jl
export DomainDifference, DomainUnion, DomainIntersection, RevolvedDomain, RotatedDomain,
    TranslatedDomain
export rotate, revolve

# from generic/collection.jl
export DomainCollection

# from generic/derived_domain.jl
export DerivedDomain
export superdomain

# from generic/mapped_domain.jl
export MappedDomain


## Specific domains

# from domains/trivial.jl
export EmptyDomain, EuclideanSpace
# from domains/interval.jl
export Interval
# from domains/simple.jl
export UnitBall, Disk, Ball, Cube
export rectangle, cube, cylinder, randomcircles
# from domains/fractals.jl
export Mandelbrot, JuliaSet
# from domains/characteristic.jl
export Characteristic
# from domains/atomium.jl
export atomium


include("util/common.jl")
include("util/tensorproducts.jl")
include("util/box.jl")

include("maps/maps.jl")
include("maps/productmap.jl")
include("maps/composite_map.jl")
include("maps/affine_map.jl")
include("maps/basic_maps.jl")
include("maps/coordinates.jl")

include("generic/domain.jl")
include("generic/productdomain.jl")
include("generic/arithmetics.jl")
include("generic/collection.jl")
include("generic/derived_domain.jl")
include("generic/mapped_domain.jl")

include("domains/trivial.jl")
include("domains/interval.jl")
include("domains/simple.jl")
include("domains/fractals.jl")
include("domains/characteristic.jl")
include("domains/atomium.jl")

end # module
