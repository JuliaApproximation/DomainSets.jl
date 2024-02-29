
# The contents of this file are a copy of the preliminary core package
# DomainSetsCore.jl. Once that package is registered and usable,
# this file can be removed.

export Domain,
    domain,
    domaineltype,
    DomainStyle,
    IsDomain,
    NotDomain,
    DomainRef,
    AnyDomain,
    checkdomain

# AsDomain was defined in versions 0.7 and 0.7.1, and then replaced by DomainRef
@deprecate AsDomain(d) DomainRef(d)

"""
    domaineltype(d)

The `domaineltype` of a continuous set is a valid type for elements of that set.
By default it is equal to the `eltype` of `d`, which in turn defaults to `Any`.
"""
domaineltype(d) = eltype(d)


# Definition of Domain commented out as it is defined in IntervalSets.jl
# """
# A `Domain{T}` is a supertype for domains with `domaineltype` equal to `T`.
#
# In addition, the `eltype` of a `Domain{T}` is also equal to `T`.
# """
# abstract type Domain{T} end
#
domaineltype(d::Domain{T}) where T = T
Base.eltype(::Type{<:Domain{T}}) where T = T


"""
    DomainStyle(d)

Trait to indicate whether or not `d` implements the domain interface.
"""
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


DomainStyle(d) = DomainStyle(typeof(d))
# - the default is no domain
DomainStyle(::Type) = NotDomain()
# - subtypes of Domain are domains
DomainStyle(::Type{<:Domain}) = IsDomain()
# - declare Number, AbstractSet and AbstractArray to be valid domain types
DomainStyle(::Type{<:Number}) = IsDomain()
DomainStyle(::Type{<:AbstractSet}) = IsDomain()
DomainStyle(::Type{<:AbstractArray}) = IsDomain()

BaseDomainType = Union{<:Number,<:AbstractSet,<:AbstractArray}

"""
    domain(d)

Return a domain associated with the object `d`.
"""
domain(d::Domain) = d

"""
    DomainRef(d)

A reference to a domain.

In a function call, `DomainRef(x)` can be used to indicate that `x` should be
treated as a domain, e.g., `foo(x, DomainRef(d))`.
"""
struct DomainRef{D}
    domain  ::  D
end

domain(d::DomainRef) = d.domain
domaineltype(d::DomainRef) = domaineltype(domain(d))


"""
`AnyDomain` is the union of `Domain` and `DomainRef`.

In both cases `domain(d::AnyDomain)` returns the domain itself.
"""
const AnyDomain = Union{Domain,DomainRef}

"""
   checkdomain(d)

Checks that `d` is a domain or refers to a domain and if so returns that domain,
throws an error otherwise.
"""
checkdomain(d::Domain) = d
# we trust the explicit intention of a user providing a domain reference
checkdomain(d::DomainRef) = domain(d)
# for other objects we check DomainStyle
checkdomain(d) = _checkdomain(d, DomainStyle(d))
_checkdomain(d, ::IsDomain) = d
_checkdomain(d, ::NotDomain) =
    error("Domain does not implement domain interface as indicated by DomainStyle.")
