# collection.jl

"""
A collection of domains.
"""
type DomainCollection{T} <: Domain{T}
    list    ::  Array{Domain{T},1}
end

DomainCollection(d::Domain) = DomainCollection([d])

length(d::DomainCollection) = length(d.list)

domain(d::DomainCollection, i) = d.list[i]

# Iteration over the domain list
start(d::DomainCollection) = start(d.list)

next(d::DomainCollection, state) = next(d.list, state)

done(d::DomainCollection, state) = done(d.list, state)

function indomain(x, dc::DomainCollection)
    z = false
    for d in dc
        z = z || indomain(x, d)
    end
    z
end

push!(dc::DomainCollection, d::Domain) = push!(dc.list, d)

show(io::IO, d::DomainCollection) = print(io, "a collection of ", length(d.list), " domains")
