# test_generic_domain.jl

function test_generic_domains()
    domains = [
        interval(),
        cube()
    ]

    @testset "$(rpad("Generic domains",80))" begin
        for domain in domains
            test_generic_domain(domain)
        end
    end
end

function test_generic_domain(d::Domain)

end
