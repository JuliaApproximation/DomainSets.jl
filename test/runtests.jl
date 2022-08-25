using Test, LinearAlgebra, StaticArrays, Random, StableRNGs

using DomainSets
using CompositeTypes.Indexing

include("test_common.jl")
include("test_maps.jl")
include("test_generic_domain.jl")
include("test_specific_domains.jl")
include("test_canonical.jl")
include("test_setoperations.jl")
include("test_applications.jl")
include("test_readme.jl")
