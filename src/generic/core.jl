export Domain,
    domain,
    domaineltype,
    DomainStyle,
    IsDomain,
    NotDomain,
    AsDomain,
    AnyDomain,
    checkdomain

# """
# A domain is a set of elements that is possibly continuous.
#
# Examples may be intervals and triangles. These are geometrical shapes, but more
# generally a domain can be any type that supports `in`. Conceptually, a domain
# is the set of all elements `x` for which `in(x, domain)` returns true.
#
# A `Domain{T}` is a domain with elements of type `T`, in analogy with
# `AbstractSet{T}` and `AbstractVector{T}`. Although domains may be defined by a
# mathematical condition such as `a <= x <= b`, irrespective of the type of `x`,
# points generated to belong to the domain have type `T`.
# """
# abstract type Domain{T} end

import IntervalSets: Domain
Base.eltype(::Type{<:Domain{T}}) where {T} = T

"""
    domaineltype(d)

The `domaineltype` of a continuous domain is the element type of any
discretization of that domain.
"""
domaineltype(d) = domaineltype(typeof(d))
domaineltype(::Type{D}) where D = eltype(D)
domaineltype(d::Domain{T}) where T = T      # shortcut definition


abstract type DomainStyle end

"""
    IsDomain()

indicates an object implements the domain interface.
"""
struct IsDomain <: DomainStyle end

"""
    NotDomain()

indicates an object does not implement the domain interface.
"""
struct NotDomain <: DomainStyle end


"""
    DomainStyle(d)

The domain style of `d` is a trait to indicate whether or not `d` implements
the domain interface. If so, the object supports `in(x,d)` and optionally also
implements `domaineltype(d)`.
"""
DomainStyle(d) = DomainStyle(typeof(d))
DomainStyle(::Type) = NotDomain()
DomainStyle(::Type{<:Domain}) = IsDomain()
DomainStyle(::Type{<:Number}) = IsDomain()
DomainStyle(::Type{<:AbstractSet}) = IsDomain()
DomainStyle(::Type{<:AbstractArray}) = IsDomain()

"""
A reference to a domain.

In a function call, `AsDomain(x)` can be used to indicate that `x` should be
treated as a domain in the function, e.g., `foo(x, AsDomain(d))`.
"""
abstract type AsDomain end

"A reference to a specific given domain."
struct DomainRef{D} <: AsDomain
    domain  ::  D
end
domain(d::DomainRef) = d.domain
domain(d::Domain) = d

Base.eltype(::Type{<:DomainRef{D}}) where D = domaineltype(D)

domaineltype(d::AsDomain) = domaineltype(domain(d))
domaineltype(::Type{<:DomainRef{D}}) where D = domaineltype(D)

AsDomain(d) = DomainRef(d)
AsDomain(d::Domain) = d

"""
`AnyDomain` can be a concrete domain or a reference to a domain.

In both cases `domain(d::AnyDomain)` returns the domain itself.
"""
const AnyDomain = Union{Domain,AsDomain}

"""
   checkdomain(d)

Checks that `d` is a domain or refers to a domain and if so returns that domain,
throws an error otherwise.
"""
checkdomain(d::Domain) = d
# we trust the explicit intention of a user providing a domain reference
checkdomain(d::AsDomain) = domain(d)
# for other objects we check DomainStyle
checkdomain(d) = _checkdomain(d, DomainStyle(d))
_checkdomain(d, ::IsDomain) = d
_checkdomain(d, ::NotDomain) =
    error("Domain does not implement domain interface as indicated by DomainStyle.")
