# collection.jl

"""
A collection of domains.
"""
type DomainCollection{N} <: Domain{N}
    list    ::  Array{Domain{N},1}
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

function indomain_grid!(z, grid, dc::DomainCollection)
    for d in dc
        indomain_grid!(z, grid, d)
    end
    z
end

push!(dc::DomainCollection, d::Domain) = push!(dc.list, d)



function boundingbox(d::DomainCollection)
    ubox = boundingbox(d.list[1])
    for i = 2:length(d.list)
        ubox = union(ubox, boundingbox(d.list[i]))
    end
    ubox
end


show(io::IO, d::DomainCollection) = print(io, "a collection of ", length(d.list), " domains")
