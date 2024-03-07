module FunctionMapsTests

using Test, LinearAlgebra, StaticArrays
using CompositeTypes, CompositeTypes.Indexing

include("using_fmaps2.jl")

include("test_common.jl")
include("test_interface.jl")

include("test_generic.jl")
include("test_basic.jl")
include("test_affine.jl")
include("test_product.jl")
include("test_maps.jl")

end
