# trivial.jl

######################
## The empty domain
######################

struct EmptySpace{T} <: Domain{T}
end

const AnyEmptySpace = EmptySpace{Any}

EmptySpace() = EmptySpace{Float64}()
EmptySpace(::Type{T}) where {T} = EmptySpace{T}()

emptyspace(d::Domain{T}) where {T} = EmptySpace{T}()

indomain(x::T, d::EmptySpace{T}) where {T} = false

approx_indomain(x, d::EmptySpace, tolerance) = indomain(x, d)

isempty(d::EmptySpace) = true

# Arithmetic operations

union(d1::EmptySpace{T}, d2::EmptySpace{T}) where {T} = d1
union(d1::Domain{T}, d2::EmptySpace{T}) where {T} = d1
union(d1::EmptySpace{T}, d2::Domain{T}) where {T} = d2

intersect(d1::EmptySpace{T}, d2::EmptySpace{T}) where {T} = d1
intersect(d1::Domain{T}, d2::EmptySpace{T}) where {T} = d2
intersect(d1::EmptySpace{T}, d2::Domain{T}) where {T} = d1

setdiff(d1::EmptySpace{T}, d2::EmptySpace{T}) where {T} = d1
setdiff(d1::EmptySpace{T}, d2::Domain{T}) where {T} = d1
setdiff(d1::Domain{T}, d2::EmptySpace{T}) where {T} = d1

# TODO: verify these - should we restrict x?
(+)(d::EmptySpace, x::Number) = d
(*)(d::EmptySpace, x::Number) = d

show(io::IO, d::EmptySpace) = print(io, "the empty space with eltype ", eltype(d))


##################################
### The whole space R^N (or C^N)
##################################

struct FullSpace{T} <: Domain{T}
end

const AnyFullSpace = FullSpace{Any}

FullSpace() = FullSpace{Float64}()
FullSpace(::Type{T}) where {T} = FullSpace{T}()

fullspace(d::Domain) = FullSpace{eltype(d)}()

euclideanspace(n::Val{N}) where {N} = euclideanspace(n, Float64)
euclideanspace(::Val{N}, ::Type{T}) where {N,T} = FullSpace(SVector{N,T})

indomain(x::T, d::FullSpace{T}) where {T} = true
indomain(x::S, d::FullSpace{T}) where {T,S} = promotes_to(S,T) == Val{true}

approx_indomain(x, d::FullSpace, tolerance) = indomain(x, d)

# We choose the origin as a point in the full space
point_in_domain(d::FullSpace) = zero(eltype(d))

# Arithmetic operations

union(d1::FullSpace{T}, d2::FullSpace{T}) where {T} = d1
union(d1::Domain{T}, d2::FullSpace{T}) where {T} = d2
union(d1::FullSpace{T}, d2::Domain{T}) where {T} = d1

intersect(d1::FullSpace{T}, d2::FullSpace{T}) where {T} = d1
intersect(d1::Domain{T}, d2::FullSpace{T}) where {T} = d1
intersect(d1::FullSpace{T}, d2::Domain{T}) where {T} = d2


(+)(d::FullSpace{T}, x::T) where {T} = d

(*)(d::FullSpace, x::Number) = d


show(io::IO, d::FullSpace) = print(io, "the full space with eltype ", eltype(d))


convert(::Type{Domain}, ::Type{T}) where T = FullSpace{T}()
convert(::Type{Domain{S}}, ::Type{T}) where {T,S} = convert(Domain{S}, convert(Domain, T))

infimum(d::FullSpace{T}) where {T} = typemin(T)
supremum(d::FullSpace{T}) where {T} = typemax(T)
