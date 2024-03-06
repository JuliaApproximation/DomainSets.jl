module FunctionMaps

using CompositeTypes, CompositeTypes.Display, CompositeTypes.Indexing
using LinearAlgebra
using StaticArrays

import Base:
    convert, show,
    ==,             # for the equality of maps
    ∘               # used to compose maps

import CompositeTypes: component, components

# Exhaustive list of exports:

# from util/common.jl
export prectype, numtype

# from generic/interface.jl
export MapRef
# from generic/map.jl
export Map,
    applymap, isequalmap,
    domaintype, codomaintype,
    inverse, leftinverse, rightinverse,
    mapsize, jacobian, jacdet, diffvolume
# from generic/composite.jl
export composedmap
# from generic/product.jl
export productmap

# from types/basic.jl
export IdentityMap,
    StaticIdentityMap, VectorIdentityMap,
    ZeroMap, UnityMap, ConstantMap,
    isconstantmap, mapconstant
# from types/affine.jl
export AffineMap, Translation, LinearMap,
    affine_matrix, affine_vector,
    islinearmap, isaffinemap

include("util/common.jl")

include("generic/map.jl")
include("generic/interface.jl")
include("generic/canonical.jl")
include("generic/lazy.jl")
include("generic/inverse.jl")
include("generic/jacobian.jl")
include("generic/composite.jl")
include("generic/product.jl")
include("generic/isomorphism.jl")
include("types/basic.jl")
include("types/affine.jl")
include("types/arithmetics.jl")

include("changes.jl")

end
