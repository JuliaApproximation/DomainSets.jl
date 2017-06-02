# trivial.jl

######################
## The empty domain
######################

struct EmptySpace{N} <: Domain{N}
end

EmptySpace() = EmptySpace{1}()
EmptySpace{N}(::Val{N}) = EmptySpace{N}()

indomain(x, d::EmptySpace) = false

# Arithmetic operations

(+)(d::EmptySpace, x::Number) = d

(*)(d::EmptySpace, x::Number) = d


show(io::IO, d::EmptySpace) = print(io, "the empty space")


##################################
### The whole space R^N (or C^N)
##################################

struct FullSpace{N} <: Domain{N}
end

FullSpace() = FullSpace{1}()
FullSpace{N}(::Val{N}) = FullSpace{N}()

indomain(x, d::FullSpace) = true

# Arithmetic operations

(+)(d::FullSpace, x::Number) = d

(*)(d::FullSpace, x::Number) = d


show(io::IO, e::FullSpace) = print(io, "the ", ndims(e), "-dimensional Euclidean space")
