# setoperations.jl
# Code for computing with domains.


############################
# The union of two domains
############################

# DD can be any collection type: a Tuple, or a Vector
# Eventually it will likely be forced to be an AbstractSet
#
"""
A `UnionDomain` represents the union of a set of domains.
"""
struct UnionDomain{DD,T} <: Domain{T}
    domains  ::  DD
end

UnionDomain(domains::AbstractSet{DD}) where {DD<:Domain{T}} where {T} =
    UnionDomain{DD,T}(domains)
UnionDomain(domains::AbstractSet) = UnionDomain(map(Domain, domains))

function UnionDomain(domains::Domain...)
    # TODO: implement promote_space_type for domains and do the promotion properly
    T = eltype(domains[1])
    for d in domains
        @assert eltype(d) == T
    end
    UnionDomain{typeof(domains),T}(domains)
end

UnionDomain(domains...) = UnionDomain(Domain.(domains)...)

UnionDomain(domains::AbstractVector) =
    UnionDomain{typeof(domains),eltype(eltype(domains))}(domains)
UnionDomain(domains::Tuple) = UnionDomain(domains...)

elements(d::UnionDomain) = d.domains

hash(d::UnionDomain, h::UInt) = hash(Set(d.domains), h)

union(d1::Domain{T}, d2::Domain{T}) where {T} = d1 == d2 ? d1 : UnionDomain(d1, d2)
function union(d1::Domain{S}, d2::Domain{T}) where {S,T}
    TS = promote_type(S, T)
    union(convert(Domain{TS}, d1), convert(Domain{TS}, d2))
end

# Avoid creating nested unions
union(d1::UnionDomain, d2::UnionDomain) = UnionDomain(elements(d1)..., elements(d2)...)
union(d1::UnionDomain, d2::Domain) = UnionDomain(elements(d1)..., d2)
union(d1::Domain, d2::UnionDomain) = UnionDomain(d1, elements(d2)...)


# The union of domains corresponds to a logical OR of their characteristic functions
indomain(x, d::UnionDomain) = mapreduce(d->in(x, d), |, elements(d))
approx_indomain(x, d::UnionDomain, tol = default_tolerance(d)) = mapreduce(d->approx_in(x, d, tol), |, elements(d))

point_in_domain(d::UnionDomain) = point_in_domain(element(d,1))

==(a::UnionDomain, b::UnionDomain) = Set(elements(a)) == Set(elements(b))

function show(io::IO, d::UnionDomain)
    print(io, "a union of $(numelements(d)) domains:\n")
    for (i,e) in enumerate(elements(d))
        print(io, "\t$i.\t: ", e, "\n")
    end
end


## ualgebra

for op in (:+, :-, :*, :/)
    @eval begin
        $op(domain::UnionDomain, x::Number) = UnionDomain(broadcast($op, elements(domain), x))
        $op(x::Number, domain::UnionDomain) = UnionDomain(broadcast($op, x, elements(domain)))
    end
end

\(x::Number, domain::UnionDomain) = UnionDomain(broadcast(\, x, elements(domain)))


for (op, mop) in ((:minimum, :min), (:maximum, :max), (:infimum, :min), (:supremum, :max))
    @eval $op(d::UnionDomain) = mapreduce($op, $mop, elements(d))
end


setdiff(d1::UnionDomain, d2::UnionDomain) = UnionDomain(setdiff.(elements(d1), d2))

function setdiff(d1::UnionDomain, d2::Domain)
    s = Set(elements(d1))
    # check if any element is in d1 and just remove
    s2 = Set(setdiff(s, tuple(d2)))
    s2 ≠ s && return UnionDomain(s2)

    UnionDomain(setdiff.(elements(d1), d2))
end

function setdiff(d1::Domain, d2::UnionDomain)
    ret = d1
    for d in elements(d2)
        ret = setdiff(ret, d)
    end
    ret
end

# use \ as a synomym for setdiff, in the context of domains (though, generically,
# \ means left division in Julia)
\(d1::Domain, d2::Domain) = setdiff(d1, d2)

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
    if numelements(d1) == numelements(d2)
        ProductDomain([intersect(element(d1,i), element(d2,i)) for i in 1:numelements(d1)]...)
    else
        IntersectionDomain(d1, d2)
    end
end


function show(io::IO, d::IntersectionDomain)
    print(io, "the intersection of $(numelements(d)) domains:\n")
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

function setdiff(d1::Domain{T}, d2::Domain) where T
    d1 == d2 && return EmptySpace{T}()
    DifferenceDomain(d1, d2)
end

setdiff(d1::Domain, d2) = setdiff(d1, convert(Domain, d2))

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
indomain(x, d::DifferenceDomain) = indomain(x, d.d1) && (~indomain(x, d.d2))

\(d1::Domain, d2) = setdiff(d1, d2)


function show(io::IO, d::DifferenceDomain)
    print(io, "the difference of 2 domains:\n")
    print(io, "\t1.\t: ", d.d1, "\n")
    print(io, "\t2.\t: ", d.d2, "\n")
end
