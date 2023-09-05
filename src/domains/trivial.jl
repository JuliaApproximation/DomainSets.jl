
# The trivial domains include the empty space and the full space.

"The empty space."
struct EmptySpace{T} <: Domain{T}
end

const AnyEmptySpace = EmptySpace{Any}

EmptySpace() = EmptySpace{Float64}()
EmptySpace(::Type{T}) where {T} = EmptySpace{T}()

similardomain(::EmptySpace, ::Type{T}) where {T} = EmptySpace{T}()

"Return the empty space with the same element type as the given domain."
emptyspace(d) = emptyspace(domaineltype(d))
emptyspace(::Type{T}) where {T} = EmptySpace{T}()

indomain(x::T, d::EmptySpace{T}) where {T} = false
approx_indomain(x, d::EmptySpace, tolerance) = in(x, d)

show(io::IO, d::EmptySpace) = print(io, "{} (empty domain)")

isempty(d::EmptySpace) = true

isopenset(d::EmptySpace) = true
isclosedset(d::EmptySpace) = true

boundary(d::EmptySpace) = d
interior(d::EmptySpace) = d
closure(d::EmptySpace) = d
boundingbox(d::EmptySpace) = d

distance_to(d::EmptySpace, x) = convert(prectype(d), Inf)

# Arithmetic operations

issubset1(d1::EmptySpace, d2) = true
issubset2(d1, d2::EmptySpace) = isempty(d1)

# setdiffdomain2(d1, d2::EmptySpace) = d2

map_domain(map::Map{T}, d::EmptySpace{T}) where {T} = d
mapped_domain(map::Map, d::EmptySpace) = EmptySpace{codomaintype(map)}()

isequaldomain(d1::EmptySpace, d2::EmptySpace) = true
isequaldomain1(d1::EmptySpace, d2) = isempty(d2)
isequaldomain2(d1, d2::EmptySpace) = isempty(d1)
hash(d::EmptySpace, h::UInt) = hash("EmptySpace", h)


"""
A domain that represents the full space.

The element type `T` in `FullSpace{T}` should only be seen as an indication of
the expected types of the elements in the context where the domain is intended
to be used. Due to the default, loose interpretation of `T`, any `FullSpace{T}`
actually contains any `x` regardless of the type of `x`. For a strict domain
of all elements of type `T`, or elements convertible exactly to `T`, use
`TypeDomain{T}`.
"""
struct FullSpace{T} <: Domain{T} end

const AnyFullSpace = FullSpace{Any}

FullSpace() = FullSpace{Float64}()
FullSpace(d) = FullSpace{domaineltype(d)}()

"Return the full space with the same element type as the given domain."
fullspace(d) = fullspace(domaineltype(d))
fullspace(::Type{T}) where {T} = FullSpace{T}()

isfullspace(d::FullSpace) = true
isfullspace(d::Domain) = false
isfullspace(d) = _isfullspace(d, DomainStyle(d))
_isfullspace(d, ::IsDomain) = false
_isfullspace(d, ::NotDomain) = error("isfullspace invoked on non-domain type")

similardomain(::FullSpace, ::Type{T}) where {T} = FullSpace{T}()

euclideanspace(n::Val{N}) where {N} = euclideanspace(n, Float64)
euclideanspace(::Val{N}, ::Type{T}) where {N,T} = FullSpace{SVector{N,T}}()

indomain(x::T, d::FullSpace{T}) where {T} = true
approx_indomain(x, d::FullSpace, tolerance) = in(x, d)

show(io::IO, d::FullSpace) = print(io, "{x} (full space)")

# We choose the origin as a point in the full space
choice(d::FullSpace) = zero(eltype(d))

isempty(::FullSpace) = false

isopenset(d::FullSpace) = true
isclosedset(d::FullSpace) = true

boundary(d::FullSpace{T}) where {T} = EmptySpace{T}()
interior(d::FullSpace) = d
closure(d::FullSpace) = d

distance_to(d::FullSpace, x) = zero(prectype(d))

# Arithmetic operations

issubset1(d1::FullSpace, d2) = isfullspace(d2)
issubset2(d1, d2::FullSpace) = true

map_domain(m::AbstractAffineMap{T}, d::FullSpace{T}) where {T} = d

isequaldomain(d1::FullSpace, d2::FullSpace) = true
isequaldomain1(d1::FullSpace, d2) = isfullspace(d2)
isequaldomain2(d1, d2::FullSpace) = isfullspace(d1)
hash(d::FullSpace, h::UInt) = hash("FullSpace", h)


convert(::Type{Domain}, ::Type{T}) where T = FullSpace{T}()
convert(::Type{Domain{S}}, ::Type{T}) where {T,S} = convert(Domain{S}, convert(Domain, T))

infimum(d::FullSpace{T}) where {T} = typemin(T)
supremum(d::FullSpace{T}) where {T} = typemax(T)



"""
The domain of all objects of type `T` and all objects convertible exactly to
type `T`.
"""
struct TypeDomain{T} <: Domain{T} end

"Return the domain for the element type of the given domain."
typedomain(d) = typedomain(domaineltype(d))
typedomain(::Type{T}) where {T} = TypeDomain{T}()

iscompatiblepair(x::T, d::TypeDomain{T}) where {T} = true
iscompatiblepair(x::T, d::TypeDomain{Any}) where {T} = true
function iscompatiblepair(x::S, d::TypeDomain{T}) where {S,T}
    # There is no generic function to check whether x can be converted to T.
    # A try-catch is not optimal, but it is a failsafe backup.
    # We can at least efficiently check promotion: if x can be converted to
    # T, then S and T must have a joined supertype.
    # For specific T and S, it would be better to specialize this function.
    promote_type(S,T) != Any &&
        try
            convert(T, x) == x
        catch
            false
        end
end

# Suppress the warning about incompatibility because it is not likely to be
# relevant for this domain
compatible_or_false(x, domain::TypeDomain) = iscompatiblepair(x, domain)

indomain(x::T, d::TypeDomain{T}) where {T} = true
approx_indomain(x, d::TypeDomain, tolerance) = in(x, d)

isequaldomain(d1::TypeDomain{S}, d2::TypeDomain{T}) where {S,T} = S==T


# some special cases
iscompatiblepair(x::Irrational, d::TypeDomain{<:Real}) = true
