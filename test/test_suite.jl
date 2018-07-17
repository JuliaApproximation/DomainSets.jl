module test_suite

if VERSION < v"0.7-"
    using Base.Test
    const ComplexF64 = Complex{Float64}
else
    using Test
    using LinearAlgebra
end

using Domains
using StaticArrays

include("test_utils.jl")
include("test_spaces.jl")
include("test_maps.jl")
include("test_generic_domain.jl")
include("test_specific_domains.jl")

function run_tests()
    delimit("Spaces")
    test_spaces()

    delimit("Maps")
    test_maps()

    delimit("Generic domain tests")
    test_generic_domains()

    delimit("Specific domain tests")
    test_specific_domains()
end

run_tests()

end #module
