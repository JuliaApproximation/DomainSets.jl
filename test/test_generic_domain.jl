# test_generic_domain.jl


# We test the generic functionality of a domain.
# These tests check whether the given domain correctly implements the
# interface of a domain.
function test_generic_domain(d::Domain)
    @test eltype(eltype(d)) == subeltype(d)
    @test isreal(d) == isreal(subeltype(d))

    if !isempty(d)
        x = point_in_domain(d)
        @test x ∈ d
    else
        try
            x = point_in_domain(d)
            @test false
        catch
        end
    end
end


@testset "generic domains" begin
    domains = [
        0..1,
        UnitInterval()^3
    ]

    @testset "$(rpad("Generic domains",80))" begin
        for domain in domains
            test_generic_domain(domain)
        end
    end
end
