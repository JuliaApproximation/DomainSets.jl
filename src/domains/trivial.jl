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

struct EuclideanSpace{N} <: Domain{N}
end

EuclideanSpace() = EuclideanSpace{1}()
EuclideanSpace{N}(::Val{N}) = EuclideanSpace{N}()

indomain(x, d::EuclideanSpace) = true

# Arithmetic operations

(+)(d::EuclideanSpace, x::Number) = d

(*)(d::EuclideanSpace, x::Number) = d


show(io::IO, e::EuclideanSpace) = print(io, "the ", ndims(e), "-dimensional Euclidean space")
