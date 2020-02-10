# Routines for a cartesian product of domains

###################
# Helper functions
###################

"""
Create an eltype that is suitable for a product domain. The result is a tuple
type, where each of the elements is the eltype of the corresponding element
of the product domain.
"""
product_eltype(domains::Domain...) = product_eltype(domains)
product_eltype(domains::Tuple) = _product_eltype(Tuple{map(eltype, domains)...})
_product_eltype(t) = t
_product_eltype(::Type{Tuple{T,V}}) where {T<:Real,V<:Real} =
	NTuple{2,promote_type(T,V)}
_product_eltype(::Type{Tuple{T,V}}) where {T<:Complex,V<:Complex} =
	NTuple{2,promote_type(T,V)}


"""
Try to simplify the type of a product domain to a type to which it is isomorphic.
The goal is to automatically embed the product domain in ℝ^N if possible.

Examples of simplifications:
`Tuple{Float64,Float64} -> SVector{2,Float64}`
`Tuple{Tuple{Float64,Float64},Float64} -> SVector{3,Float64}`
"""
simplify_product_eltype(::Type{T}) where {T} = T
simplify_product_eltype(::Type{NTuple{N,T}}) where {N,T} = SVector{N,T}

@generated function simplify_product_eltype(::Type{Tuple{SVector{N,T},T}}) where {N,T}
    M = N+1
    quote SVector{$M,T} end
end
@generated function simplify_product_eltype(::Type{Tuple{T,SVector{N,T}}}) where {N,T}
    M = N+1
    quote SVector{$M,T} end
end
@generated function simplify_product_eltype(::Type{Tuple{SVector{N,T},SVector{K,T}}}) where {N,K,T}
    M = N+K
    quote SVector{$M,T} end
end

@generated function simplify_product_eltype(::Type{Tuple{NTuple{N,T},T}}) where {N,T}
    M = N+1
    quote SVector{$M,T} end
end
@generated function simplify_product_eltype(::Type{Tuple{T,NTuple{N,T}}}) where {N,T}
    M = N+1
    quote SVector{$M,T} end
end
@generated function simplify_product_eltype(::Type{Tuple{NTuple{N,T},NTuple{K,T}}}) where {N,K,T}
    M = N+K
    quote SVector{$M,T} end
end

"""
Convert a vector from a cartesian format to a nested tuple according to the
given dimensions.

For example:
`convert_fromcartesian([1,2,3,4,5], Val{(2,2,1)}()) -> ([1,2],[3,4],5)`
"""
@generated function convert_fromcartesian(x::Tuple, ::Val{DIM}) where {DIM}
	dimsum = [0; cumsum([d for d in DIM])]
	E = Expr(:tuple, [ (dimsum[i+1]-dimsum[i] > 1 ? Expr(:tuple, [:(x[$j]) for j = dimsum[i]+1:dimsum[i+1]]...) : :(x[$(dimsum[i+1])])) for i in 1:length(DIM)]...)
	return quote $(E) end
end

@generated function convert_fromcartesian(x::SVector, ::Val{DIM}) where {DIM}
	dimsum = [0; cumsum([d for d in DIM])]
	E = Expr(:tuple, [ (dimsum[i+1]-dimsum[i] > 1 ? Expr(:call, :SVector, [:(x[$j]) for j = dimsum[i]+1:dimsum[i+1]]...) : :(x[$(dimsum[i+1])])) for i in 1:length(DIM)]...)
	return quote $(E) end
end

"The inverse function of `convert_fromcartesian`."
@generated function convert_tocartesian(x, ::Val{DIM}) where {DIM}
    dimsum = [0; cumsum([d for d in DIM])]
    E = vcat([[:(x[$i][$j]) for j in 1:DIM[i]] for i in 1:length(DIM)]...)
    quote SVector($(E...)) end
end

#######################
# Main type definition
#######################



"""
```using DomainSets```

A `ProductDomain` represents the cartesian product of other domains.
"""
struct ProductDomain{T,DIM,DD} <: LazyDomain{T}
	domains	::	DD
end

const TupleProductDomain{T,DIM,DD<:Tuple} = ProductDomain{T,DIM,DD}
const VectorProductDomain{T,DD} = ProductDomain{T,0,DD}

productdimension(d::TupleProductDomain{T,DIM,DD}) where {T,DIM,DD} = DIM
productdimension(d::VectorProductDomain) = numelements(d)

ProductDomain(domains::Domain...) = ProductDomain(domains)

function ProductDomain(domains::Tuple)
    S = product_eltype(domains...)
	T = simplify_product_eltype(S)
    ProductDomain{T}(domains)
end

function ProductDomain(domains::AbstractVector)
    T = reduce(promote_type, map(eltype, domains))
    ProductDomain{Vector{T}}(domains)
end

ProductDomain{T}(domains::Domain...) where {T} = ProductDomain{T}(domains)

function ProductDomain{T}(domains::Tuple) where {T}
	DD = typeof(domains)
	DIM = map(dimension, domains)
	ProductDomain{T,DIM,DD}(domains)
end

function ProductDomain{T}(domains::Vector) where {T}
	DD = typeof(domains)
	ProductDomain{T,0,DD}(domains)
end

composition(d::ProductDomain) = Product()

preprocess(d::TupleProductDomain{SVector{N,T},DIM,DD}, x) where {N,T,DIM,DD} =
	convert_fromcartesian(x, Val{DIM}())

function preprocess(d::VectorProductDomain, x)
	@assert length(x) == productdimension(d)
	x
end

function convert(::Type{EuclideanDomain{N,T}}, d::TupleProductDomain{S,DIM}) where {N,T,S,DIM}
	@assert sum(DIM) == N
	ProductDomain(convert.(Domain{T}, d.domains)...)
end

convert(::Type{VectorDomain{T}}, d::VectorProductDomain{S}) where {S,T} =
	ProductDomain(convert.(Domain{T}, d.domains))

cross(x::Domain...) = cartesianproduct(x...)

# One can use the cartesianproduct routine to create product domains
cartesianproduct(domains::Domain...) = ProductDomain(expand(ProductDomain, domains...))

^(d::Domain, n::Int) = cartesianproduct(d, n)

point_in_domain(d::ProductDomain) = topoint(d, map(point_in_domain, elements(d)))

topoint(d::ProductDomain, x) = x
topoint(d::TupleProductDomain{<:SVector,DIM}, x) where {DIM} =
	convert_tocartesian(x, Val{DIM}())

infimum(d::ProductDomain) = topoint(d, map(infimum, elements(d)))
supremum(d::ProductDomain) = topoint(d, map(supremum, elements(d)))

isempty(d::ProductDomain) = any(isempty, d.domains)

function show(io::IO, t::ProductDomain)
    L = numelements(t)
    for i in 1:L-1
        show(io, element(t, i))
        print(io, " x ")
    end
    show(io, element(t, L))
end
