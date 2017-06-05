# trivial.jl

######################
## The empty domain
######################

struct EmptySpace{T} <: Domain{T}
end

const AnyEmptySpace = EmptySpace{Any}

EmptySpace() = EmptySpace{Float64}()
EmptySpace(::Type{T}) where {T} = EmptySpace{T}()

emptyspace(d::Domain) = EmptySpace{eltype(d)}()

indomain(x::T, d::EmptySpace{T}) where {T} = false

# Arithmetic operations

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
euclideanspace(::Val{N}, ::Type{T}) where {N,T} = FullSpace(Point{N,T})

indomain(x::T, d::FullSpace{T}) where {T} = true

# Arithmetic operations

(+)(d::FullSpace, x::Number) = d

(*)(d::FullSpace, x::Number) = d


show(io::IO, d::FullSpace) = print(io, "the full space with eltype ", eltype(d))
