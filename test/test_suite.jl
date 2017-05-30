module test_suite


using Base.Test
using Domains
using StaticArrays

include("test_utils.jl")
include("test_generic_domain.jl")
include("test_specific_domains.jl")

function run_tests()
    delimit("Generic tests")
    test_generic_domains()

    delimit("Specific domains")
    test_specific_domains()
end

run_tests()

end #module
