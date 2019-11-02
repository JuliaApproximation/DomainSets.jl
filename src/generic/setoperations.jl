# Code for computing with domains.

############################
# The union of two domains
############################

# DD can be any collection type: a Tuple, or a Vector
# Eventually it will likely be forced to be an AbstractSet
#
"""
```using DomainSets```

A `UnionDomain` represents the union of a set of domains.
"""
struct UnionDomain{DD,T} <: Domain{T}
    domains  ::  DD
end

"""
The `UnionDomain` constructor can be invoked in one of three ways:
- with a list of arguments: UnionDomain(d1, d2)
- with a single domain: UnionDomain(d::Domain)
- or with any iterable list of domains: UnionDomain(domains)
"""
UnionDomain(domains...) = UnionDomain(domains)
UnionDomain(d::Domain) = UnionDomain((d,))

UnionDomain(domains) = _UnionDomain(expand_union(promote_domains(domains)...))
_UnionDomain(domains) = UnionDomain{typeof(domains),eltype(first(domains))}(domains)

"Return all arguments in a tuple of domains, with each `UnionDomain` replaced by its elements."
expand_union() = ()
expand_union(domain, domains...) = (domain, expand_union(domains...)...)
expand_union(domain::UnionDomain, domains...) = (elements(domain)..., expand_union(domains...)...)

convert(::Type{Domain}, v::AbstractVector{<:Domain}) = UnionDomain(v)
convert(::Type{Domain}, v::AbstractSet{<:Domain}) = UnionDomain(v)

convert(::Type{Domain{T}}, d::UnionDomain{S}) where {S,T} =
    UnionDomain(convert_domain.(T, elements(d)))

Domain(v::AbstractVector) = convert(Domain, v)
Domain(v::AbstractSet) = convert(Domain, v)

elements(d::UnionDomain) = d.domains

hash(d::UnionDomain, h::UInt) = hash(Set(d.domains), h)

union(domains::Domain...) = UnionDomain(domains...)

# Catch one specific case
union(d1::Domain, d2::Domain)  = d1 == d2 ? d1 : UnionDomain(d1, d2)


# The union of domains corresponds to a logical OR of their characteristic functions
indomain(x, d::UnionDomain) = mapreduce(d->in(x, d), |, elements(d))
approx_indomain(x, d::UnionDomain, tol = default_tolerance(d)) = mapreduce(d->approx_in(x, d, tol), |, elements(d))

point_in_domain(d::UnionDomain) = point_in_domain(element(d,1))

==(a::UnionDomain, b::UnionDomain) = Set(elements(a)) == Set(elements(b))

isempty(d::UnionDomain) = all(isempty, d.domains)

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


setdiff(d1::UnionDomain, d2::UnionDomain) = UnionDomain(setdiff.(elements(d1), Ref(d2)))

function setdiff(d1::UnionDomain, d2::Domain)
    s = Set(elements(d1))
    # check if any element is in d1 and just remove
    s2 = Set(setdiff(s, tuple(d2)))
    s2 ≠ s && return UnionDomain(s2)

    UnionDomain(setdiff.(elements(d1), Ref(d2)))
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


##############################
# The intersection of domains
##############################

"""
An `IntersectionDomain` represents the intersection of a set of domains.
"""
struct IntersectionDomain{DD,T} <: Domain{T}
    domains ::  DD
end

"""
The `IntersectionDomain` constructor can be invoked in one of three ways:
- with a list of arguments: IntersectionDomain(d1, d2)
- with a single domain: IntersectionDomain(d::Domain)
- or with any iterable list of domains: IntersectionDomain(domains)
"""
IntersectionDomain(domains...) = IntersectionDomain(domains)
IntersectionDomain(d::Domain) = IntersectionDomain((d,))

IntersectionDomain(domains) = _IntersectionDomain(expand_intersection(promote_domains(domains)...))
_IntersectionDomain(domains) = IntersectionDomain{typeof(domains),eltype(first(domains))}(domains)

"Return all arguments in a tuple of domains, with each `UnionDomain` replaced by its elements."
expand_intersection() = ()
expand_intersection(domain, domains...) = (domain, expand_intersection(domains...)...)
expand_intersection(domain::IntersectionDomain, domains...) = (elements(domain)..., expand_intersection(domains...)...)

elements(d::IntersectionDomain) = d.domains

intersect(d1::Domain, d2::Domain) = d1 == d2 ? d1 : IntersectionDomain(d1, d2)
function intersect(d1::UnionDomain, d2::UnionDomain)
    d1 == d2 && return d1
    union(intersect.(Ref(d1), elements(d2))...)
end
intersect(d1::UnionDomain, d2::Domain) = union(intersect.(d1.domains, Ref(d2))...)
intersect(d1::Domain, d2::UnionDomain) = union(intersect.(Ref(d1), d2.domains)...)

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

convert(::Type{Domain{T}}, d::IntersectionDomain{S}) where {S,T} =
    IntersectionDomain(convert_domain.(T, elements(d)))


function show(io::IO, d::IntersectionDomain)
    print(io, "the intersection of $(numelements(d)) domains:\n")
    print(io, "\t1.\t: ", element(d,1), "\n")
    print(io, "\t2.\t: ", element(d,2), "\n")
end


#########################################
### The difference between two domains
#########################################

struct DifferenceDomain{D1,D2,T} <: Domain{T}
    d1    ::  D1
    d2    ::  D2
end

DifferenceDomain(d1, d2) = _DifferenceDomain(promote_domain(d1, d2)...)
_DifferenceDomain(d1, d2) = DifferenceDomain{typeof(d1),typeof(d2),eltype(d1)}(d1, d2)

convert(::Type{Domain{T}}, d::DifferenceDomain{S}) where {S,T} =
    DifferenceDomain(convert_domain(T, d.d1), convert_domain(T, d.d2))

function setdiff(d1::Domain{T}, d2::Domain) where T
    d1 == d2 && return EmptySpace{T}()
    DifferenceDomain(d1, d2)
end

setdiff(d1::Domain, d2) = DifferenceDomain(d1, d2)
setdiff(d1, d2::Domain) = DifferenceDomain(d1, d2)

# The difference between two domains corresponds to a logical AND NOT of their characteristic functions
indomain(x, d::DifferenceDomain) = in(x, d.d1) && (~in(x, d.d2))

\(d1::Domain, d2) = setdiff(d1, d2)



function show(io::IO, d::DifferenceDomain)
    print(io, "the difference of 2 domains:\n")
    print(io, "\t1.\t: ", d.d1, "\n")
    print(io, "\t2.\t: ", d.d2, "\n")
end
