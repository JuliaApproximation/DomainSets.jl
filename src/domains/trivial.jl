# trivial.jl

######################
## The empty domain
######################

struct EmptyDomain{N} <: Domain{N}
end

EmptyDomain() = EmptyDomain{1}()
EmptyDomain{N}(::Val{N}) = EmptyDomain{N}()

indomain(x, d::EmptyDomain) = false

# Arithmetic operations

(+)(d::EmptyDomain, x::Number) = d

(*)(d::EmptyDomain, x::Number) = d


show(io::IO, d::EmptyDomain) = print(io, "the empty domain")


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
