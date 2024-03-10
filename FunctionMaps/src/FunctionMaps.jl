module FunctionMaps

using CompositeTypes, CompositeTypes.Display, CompositeTypes.Indexing
using LinearAlgebra
using StaticArrays

import Base:
    convert, show,
    ==, hash,       # for the equality of maps
    ∘,              # used to compose maps
    isreal          # to check whether maps are real

import CompositeTypes: component, components

# Exhaustive list of exports:

# from util/common.jl
export prectype, numtype

# from generic/map.jl
export Map, MapRef,
    applymap, isequalmap,
    domaintype, codomaintype,
    inverse, leftinverse, rightinverse,
    mapsize, jacobian, jacdet, diffvolume,
    isreal
# from generic/canonical.jl
export canonicalmap
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
    affinematrix, affinevector,
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
include("types/coordinates.jl")
include("types/arithmetics.jl")

include("deprecated.jl")

end
