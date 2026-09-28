module DomainSetsRecipesBaseExt

using DomainSets, RecipesBase

# a union is plotted by plotting each of its components, which reuse the color of the
# first and have no legend entry of their own
@recipe function f(d::UnionDomain)
    for (k, c) in enumerate(components(d))
        @series begin
            primary := k == 1
            c
        end
    end
end

end # module
