using Test, LinearAlgebra, StaticArrays, Random, StableRNGs

using DomainSets, FunctionMaps
using CompositeTypes.Indexing

const io = IOBuffer()
const textmime = MIME"text/plain"()

println("#############################")
println("# Aqua automated tests")
println("#############################")
include("aqua.jl")

include("test_common.jl")
include("test_generic_domain.jl")
include("test_specific_domains.jl")
include("test_canonical.jl")
include("test_setoperations.jl")
include("test_interface.jl")
include("test_applications.jl")
include("test_readme.jl")

if isdefined(Base, :get_extension)
    println("#############################")
    println("# Tests of extensions")
    println("#############################")

    include("test_makieext.jl")
end