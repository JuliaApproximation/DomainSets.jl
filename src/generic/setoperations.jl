# setoperations.jl
# Code for computing with domains.


############################
# The union of two domains
############################

"""
A `UnionDomain` represents the union of a set of domains.
"""
# DD can be any collection type: a Tuple, or a Vector
struct UnionDomain{DD,T} <: Domain{T}
    domains  ::  DD
end

function UnionDomain(domains::Domain...)
    # TODO: implement promote_space_type for domains and to the promotion properly
    T = eltype(domains[1])
    for d in domains
        @assert eltype(d) == T
    end
    UnionDomain{typeof(domains),T}(domains)
end

UnionDomain(domains::AbstractVector) = UnionDomain{typeof(domains),eltype(eltype(domains))}(domains)
UnionDomain(domains::Tuple) = UnionDomain(domains...)

elements(d::UnionDomain) = d.domains

union(d1::Domain{T}, d2::Domain{T}) where {T} = d1 == d2 ? d1 : UnionDomain(d1, d2)

# Avoid creating nested unions
union(d1::UnionDomain, d2::Domain) = UnionDomain(elements(d1)..., d2)
union(d1::UnionDomain, d2::UnionDomain) = UnionDomain(elements(d1)..., elements(d2)...)
union(d1::Domain, d2::UnionDomain) = UnionDomain(d1, elements(d2)...)


# The union of domains corresponds to a logical OR of their characteristic functions
indomain(x, d::UnionDomain) = mapreduce(d->in(x, d), |, elements(d))

(|)(d1::Domain, d2::Domain) = union(d1, d2)


function ==(a::UnionDomain, b::UnionDomain)
    length(elements(a)) ≠ length(elements(b)) && return false
    for (c, d) in zip(a.domains, b.domains)
        c ≠ d && return false
    end
    return true
end

function show(io::IO, d::UnionDomain)
    print(io, "a union of $(nb_elements(d)) domains:\n")
    for i=1:nb_elements(d)
        print(io, "\t$i.\t: ", element(d,i), "\n")
    end
end


###################################
# The intersection of two domains
###################################

"""
An `IntersectionDomain` represents the intersection of a set of domains.
"""
struct IntersectionDomain{DD,T} <: Domain{T}
    domains ::  DD
end

function IntersectionDomain(domains::Domain...)
    # TODO: implement promote_space_type for domains and to the promotion properly
    T = eltype(domains[1])
    for d in domains
        @assert eltype(d) == T
    end
    IntersectionDomain{typeof(domains),T}(domains)
end

elements(d::IntersectionDomain) = d.domains

intersect(d1::Domain{T}, d2::Domain{T}) where {T} = d1 == d2 ? d1 : IntersectionDomain(d1, d2)

# Avoid creating nested unions
intersect(d1::IntersectionDomain, d2::Domain) = IntersectionDomain(elements(d1)..., d2)
intersect(d1::IntersectionDomain, d2::IntersectionDomain) = IntersectionDomain(elements(d1)..., elements(d2)...)
intersect(d1::Domain, d2::IntersectionDomain) = IntersectionDomain(d1, elements(d2)...)

# The intersection of domains corresponds to a logical AND of their characteristic functions
indomain(x, d::IntersectionDomain) = mapreduce(d->in(x, d), &, elements(d))

(&)(d1::Domain, d2::Domain) = intersect(d1,d2)

function intersect(d1::ProductDomain, d2::ProductDomain)
    if nb_elements(d1) == nb_elements(d2)
        ProductDomain([intersect(element(d1,i), element(d2,i)) for i in 1:nb_elements(d1)]...)
    else
        IntersectionDomain(d1, d2)
    end
end


function show(io::IO, d::IntersectionDomain)
    print(io, "the intersection of $(nb_elements(d)) domains:\n")
    print(io, "\t1.\t: ", element(d,1), "\n")
    print(io, "\t2.\t: ", element(d,2), "\n")
end


################################################################################
### The difference between two domains
################################################################################

struct DifferenceDomain{D1,D2,T} <: Domain{T}
    d1    ::  D1
    d2    ::  D2

    DifferenceDomain{D1,D2,T}(d1::Domain{T}, d2::Domain{T}) where {D1,D2,T} = new{D1,D2,T}(d1, d2)
end

DifferenceDomain(d1::Domain{T}, d2::Domain{T}) where {T} = DifferenceDomain{typeof(d1),typeof(d2),T}(d1,d2)

setdiff(d1::Domain, d2::Domain) = DifferenceDomain(d1, d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
indomain(x, d::DifferenceDomain) = indomain(x, d.d1) && (~indomain(x, d.d2))

(-)(d1::Domain, d2::Domain) = setdiff(d1, d2)
(\)(d1::Domain, d2::Domain) = setdiff(d1, d2)


function show(io::IO, d::DifferenceDomain)
    print(io, "the difference of 2 domains:\n")
    print(io, "\t1.\t: ", d.d1, "\n")
    print(io, "\t2.\t: ", d.d2, "\n")
end
