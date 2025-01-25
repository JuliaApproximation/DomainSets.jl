
@deprecate convert_prectype(d::Domain, ::Type{T}) where {T} convert_prectype(T, d)
@deprecate convert_numtype(d::Domain, ::Type{T}) where {T} convert_numtype(T, d)

@deprecate point_in_domain(d) choice(d)

@deprecate broadcast_in(A, d::Domain) vectorized_in(A, d)
@deprecate broadcast_approx_in(A, d::Domain) vectorized_approx_in(A, d)
@deprecate broadcast_approx_in(A, d::Domain, tol) vectorized_approx_in(A, d, tol)

@deprecate UnionDomain(domain::Domain) UnionDomain((domain,))
@deprecate UnionDomain{T}(domain::Domain) where {T} UnionDomain{T}((domain,))
@deprecate IntersectDomain(domain::Domain) IntersectDomain((domain,))
@deprecate IntersectDomain{T}(domain::Domain) where T IntersectDomain((domain,))

@deprecate issubset(d1::AnyDomain, d2) issubset(d1, DomainRef(d2))
@deprecate union(d1::AnyDomain, d2, domains...) union(d1, DomainRef(d2), domains...)
@deprecate union(d1, d2::AnyDomain, domains...) union(DomainRef(d1), d2, domains...)
@deprecate intersect(d1, d2::AnyDomain) intersect(DomainRef(d1), d2)
@deprecate intersect(d1::AnyDomain, d2) intersect(d1, DomainRef(d2))
@deprecate \(d1::AnyDomain, d2) d1 \ DomainRef(d2)
@deprecate \(d1, d2::AnyDomain) DomainRef(d1) \ d2
@deprecate setdiff(d1::AnyDomain, d2) setdiff(d1, DomainRef(d2))
@deprecate setdiff(d1, d2::AnyDomain) setdiff(DomainRef(d1), d2)
@deprecate (&)(d1::AnyDomain, d2) d1 & DomainRef(d2)
@deprecate (&)(d1, d2::AnyDomain) DomainRef(d1) & d2

@deprecate isequal1(d1,d2) isequaldomain1(d1,d2)
@deprecate isequal2(d1,d2) isequaldomain2(d1,d2)

# AsDomain was defined in versions 0.7 and 0.7.1, and then replaced by DomainRef
@deprecate AsDomain(d) DomainRef(d)
