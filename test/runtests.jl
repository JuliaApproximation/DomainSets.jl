using Test, LinearAlgebra, StaticArrays, Random, StableRNGs

using DomainSets
using CompositeTypes.Indexing

const io = IOBuffer()
const textmime = MIME"text/plain"()

# First run the test suite of submodule FunctionMaps.jl
println("#############################")
println("# Tests of FunctionMaps.jl")
println("#############################")
include("../FunctionMaps/test/runtests2.jl")

println("#############################")
println("# Tests of DomainSets.jl")
println("#############################")
include("test_common.jl")
include("test_generic_domain.jl")
include("test_specific_domains.jl")
include("test_canonical.jl")
include("test_setoperations.jl")
include("test_interface.jl")
include("test_applications.jl")
include("test_readme.jl")
