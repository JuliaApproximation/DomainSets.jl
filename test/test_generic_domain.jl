# test_generic_domain.jl

function test_generic_domains()
    domains = [
        interval(),
        cube(),
    ]

    @testset "$(rpad("Generic domains",80))" begin
        for domain in domains
            test_generic_domain(domain)
        end
    end
end

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
