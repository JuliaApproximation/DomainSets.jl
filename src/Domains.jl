# Domains.jl

module Domains

# We use static vectors internally
using StaticArrays


################################
## Exhaustive list of imports
################################

import Base: *, +, -, /, \, |, &, ∪
import Base: ==
import Base: inv

import Base: in

import Base: show

import Base: ndims, getindex

################################
## Exhaustive list of exports
################################

# from util/common.jl
export elements, element, nb_elements, composite_length
# from util/tensorproducts.jl
export flatten, tensorproduct, ⊗
# from util/box.jl
export BBox, BBox1, BBox2, BBox3, BBox4

# from maps/maps.jl
export AbstractMap, forward_map, inverse_map, jacobian, linearize, is_linear
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


# from generic/domain.jl
export Domain, indomain, boundingbox
# from generic/productdomain.jl
export ProductDomain, tensorproduct, ⊗

# from domains/trivial.jl
export EmptyDomain, EuclideanSpace
# from domains/interval.jl
export Interval
# from domains/simple.jl
export Ball, Cube


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
