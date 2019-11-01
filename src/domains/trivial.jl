
"The empty domain with elements of type `T`."
struct EmptyDomain{T} <: Domain{T}
end

# For temporary backwards compatibility, to be removed later
const EmptySpace = EmptyDomain

const AnyEmptyDomain = EmptyDomain{Any}

EmptyDomain() = EmptyDomain{Float64}()
EmptyDomain(::Type{T}) where {T} = EmptyDomain{T}()

EmptyDomain(d::Domain{T}) where {T} = EmptyDomain{T}()

indomain(x::T, d::EmptyDomain{T}) where {T} = false

approx_indomain(x, d::EmptyDomain, tolerance) = in(x, d)

isempty(d::EmptyDomain) = true

# Arithmetic operations

union(d1::EmptyDomain, d2::EmptyDomain) = d1
union(d1::Domain, d2::EmptyDomain) = d1
union(d1::EmptyDomain, d2::Domain) = d2

intersect(d1::EmptyDomain, d2::EmptyDomain) = d1
intersect(d1::Domain, d2::EmptyDomain) = d2
intersect(d1::EmptyDomain, d2::Domain) = d1

setdiff(d1::EmptyDomain, d2::EmptyDomain) = d1
setdiff(d1::EmptyDomain, d2::Domain) = d1
setdiff(d1::Domain, d2::EmptyDomain) = d1

# TODO: verify these - should we restrict x?
(+)(d::EmptyDomain, x::Number) = d
(*)(d::EmptyDomain, x::Number) = d

==(::EmptyDomain, ::EmptyDomain) = true

show(io::IO, d::EmptyDomain) = print(io, "{} (empty domain)")


"The full space of elements of type `T`."
struct FullSpace{T} <: Domain{T} end

const AnyFullSpace = FullSpace{Any}

FullSpace() = FullSpace{Float64}()
FullSpace(d) = FullSpace{eltype(d)}()

euclideanspace(n::Val{N}) where {N} = euclideanspace(n, Float64)
euclideanspace(::Val{N}, ::Type{T}) where {N,T} = FullSpace{SVector{N,T}}()

indomain(x::T, d::FullSpace{T}) where {T} = true
indomain(x::S, d::FullSpace{T}) where {T,S} = promotes_to(S,T) == Val{true}

approx_indomain(x, d::FullSpace, tolerance) = in(x, d)

# We choose the origin as a point in the full space
point_in_domain(d::FullSpace) = zero(eltype(d))

isempty(::FullSpace) = false # constains zero

# Arithmetic operations

union(d1::FullSpace, d2::FullSpace) = d1
union(d1::Domain, d2::FullSpace) = d2
union(d1::FullSpace, d2::Domain) = d1

intersect(d1::FullSpace, d2::FullSpace) = d1
intersect(d1::Domain, d2::FullSpace) = d1
intersect(d1::FullSpace, d2::Domain) = d2


(+)(d::FullSpace{T}, x::T) where {T} = d

(*)(d::FullSpace, x::Number) = d


show(io::IO, d::FullSpace{T}) where {T} = print(io, "{x} (full space)")


convert(::Type{Domain}, ::Type{T}) where T = FullSpace{T}()
convert(::Type{Domain{S}}, ::Type{T}) where {T,S} = convert(Domain{S}, convert(Domain, T))

infimum(d::FullSpace{T}) where {T} = typemin(T)
supremum(d::FullSpace{T}) where {T} = typemax(T)
