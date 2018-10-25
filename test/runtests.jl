using Test, LinearAlgebra, DomainSets, StaticArrays
import DomainSets: elements, TypeFactory

include("test_spaces.jl")
include("test_maps.jl")
include("test_generic_domain.jl")
include("test_specific_domains.jl")
