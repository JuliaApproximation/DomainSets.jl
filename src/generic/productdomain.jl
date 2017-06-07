# productdomain.jl
# Routines for a tensor product of domains

###################
# Helper functions
###################

"""
Create an eltype that is suitable for a product domain. By default, a tuple type
that contains the element types is returned. However, if the types are homogeneous,
they are collected into an SVector.

for example:
`product_eltype(::Domain{Int}, ::Domain{Float64}) -> Tuple{Int,Float64}`
`product_eltype(::Domain{Float64}, ::Domain{Float64}) -> SVector{2,Float64}`
"""
product_eltype(domains::Domain...) = Tuple{map(eltype, domains)...}

simplify_eltype(::Type{T}) where {T} = T
simplify_eltype(::Type{NTuple{N,T}}) where {N,T} = SVector{N,T}
simplify_eltype(::Type{Tuple{Tuple{T,T},T}}) where {T} = SVector{3,T}


#######################
# Main type definition
#######################

"""
A `ProductDomain` represents the tensor product of other domains.

A product domain has two eltypes, an internal type `S` and an external type `T`.
The internal type `S` is a tuple containing the eltypes of the elements of the
product domain. The external eltype `T` is a type whose associated space is
isomorphic to that of `S`, but which has been simplified.

For example, if `S` is `Tuple{Float64,Float64}`, then `T` is `SVector{2,Float64}`.
"""
struct ProductDomain{D,S,T} <: Domain{T}
	# D is the type of an indexable list of domains, such as a tuple
	domains	::	D

	# Inner constructor to verify that S and T are correct
	function ProductDomain{D,S,T}(domains) where {D,S,T}
		@assert S == product_eltype(domains...)
		@assert isomorphic(spacetype(S),spacetype(T))
		new{D,S,T}(domains)
	end
end

function ProductDomain(domains...)
    D = typeof(domains)
    S = product_eltype(domains...)
	T = simplify_eltype(S)
    ProductDomain{D,S,T}(domains)
end

elements(d::ProductDomain) = d.domains

internal_eltype(::Type{ProductDomain{D,S,T}}) where {D,S,T} = S
internal_eltype(::Type{P}) where {P <: ProductDomain} = S
internal_eltype(d::ProductDomain) = internal_eltype(typeof(d))


# One can use the tensorproduct routine to create product domains
tensorproduct(domains::Domain...) = ProductDomain(domains...)

# We try to avoid creating nested tensor products
tensorproduct(d1::Domain, d2::ProductDomain) = tensorproduct(d1, elements(d2)...)
tensorproduct(d1::ProductDomain, d2::Domain) = tensorproduct(elements(d1)..., d2)
tensorproduct(d1::ProductDomain, d2::ProductDomain) = tensorproduct(elements(d1)..., elements(d2...))


^(d::Domain, n::Int) = tensorproduct(d, n)

indomain(x, d::ProductDomain) = _indomain(convert_space(spacetype(internal_eltype(d)), x), d)

_indomain(x, d::ProductDomain) = reduce(&, map(indomain, x, elements(d)))

# indomain(x::SVector{2}, d1::Domain{1}, d2::Domain{1}) =
# 	indomain(x[1], d1) && indomain(x[2], d2)
#
# indomain(x::SVector{3}, d1::Domain{1}, d2::Domain{1}, d3::Domain{1}) =
# 	indomain(x[1], d1) && indomain(x[2], d2) && indomain(x[3], d3)
#
# indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{1}, d3::Domain{1}, d4::Domain{1}) =
# 	indomain(x[1], d1) && indomain(x[2], d2) && indomain(x[3], d3) && indomain(x[4], d4)
#
# indomain(x::SVector{3}, d1::Domain{1}, d2::Domain{2}) =
# 	indomain(x[1], d1) && indomain(SVector(x[2],x[3]), d2)
#
# indomain(x::SVector{3}, d1::Domain{2}, d2::Domain{1}) =
# 	indomain(SVector(x[1],x[2]), d1) && indomain(x[3], d2)
#
# indomain(x::SVector{4}, d1::Domain{2}, d2::Domain{2}) =
# 	indomain(SVector(x[1],x[2]), d1) && indomain(SVector(x[3],x[4]), d2)
#
# indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{3}) =
# 	indomain(x[1], d1) && indomain(SVector(x[2],x[3],x[4]), d2)
#
# indomain(x::SVector{4}, d1::Domain{3}, d2::Domain{1}) =
# 	indomain(SVector(x[1],x[2],x[3]), d1) && indomain(x[1], d2)
#
# indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{1}, d3::Domain{2}) =
# 	indomain(x[1], d1) && indomain(x[2], d2) && indomain(SVector(x[3],x[4]), d3)
#
# indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{2}, d3::Domain{1}) =
# 	indomain(x[1], d1) && indomain(SVector(x[2],x[3]), d2) && indomain(x[4], d3)
#
# indomain(x::SVector{4}, d1::Domain{2}, d2::Domain{1}, d3::Domain{1}) =
# 	indomain(SVector(x[1],x[2]), d1) && indomain(x[3], d2) && indomain(x[4], d3)


function show(io::IO, t::ProductDomain)
    L = nb_elements(t)
    for i in 1:L-1
        show(io, element(t, i))
        print(io, " x ")
    end
    show(io, element(t, L))
end

(*)(d::ProductDomain, x::Number) = tensorproduct([domain*x for domain in elements(d)]...)
