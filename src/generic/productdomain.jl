# productdomain.jl
# Routines for a tensor product of domains

###################
# Helper functions
###################

"""
Create an eltype that is suitable for a product domain. The result is a tuple
type, where each of the elements is the eltype of the corresponding element
of the product domain.
"""
product_eltype(domains::Domain...) = Tuple{map(eltype, domains)...}

"""
Try to simplify the type of a product domain to a type to which it is isomorphic.
The goal is to automatically embed the product domain in ℝ^N if possible.

Examples of simplifications:
`Tuple{Float64,Float64} -> SVector{2,Float64}`
`Tuple{Tuple{Float64,Float64},Float64} -> SVector{3,Float64}`
"""
simplify_product_eltype(::Type{T}) where {T} = T
simplify_product_eltype(::Type{NTuple{N,T}}) where {N,T} = SVector{N,T}
simplify_product_eltype(::Type{Tuple{Tuple{T,T},T}}) where {T} = SVector{3,T}


#######################
# Main type definition
#######################

"""
A `ProductDomain` represents the tensor product of other domains.

A product domain has two eltypes, an internal type `S` and an external type `T`.
The internal type `S` is a tuple containing the eltypes of the elements of the
product domain. The external eltype `T` is a type whose associated space is
isomorphic to that of `S`, but which has been simplified. (See also
`simplify_product_eltype`).

For example, if `S` is `Tuple{Float64,Float64}`, then `T` is `SVector{2,Float64}`.
"""
struct ProductDomain{DD,S,T} <: Domain{T}
	# D is the type of an indexable list of domains, such as a tuple
	domains	::	DD

	# Inner constructor to verify that S and T are correct
	function ProductDomain{DD,S,T}(domains) where {DD,S,T}
		@assert S == product_eltype(domains...)
		@assert isomorphic(spacetype(S),spacetype(T))
		new{DD,S,T}(domains)
	end
end

function ProductDomain(domains...)
    DD = typeof(domains)
    S = product_eltype(domains...)
	T = simplify_product_eltype(S)
    ProductDomain{DD,S,T}(domains)
end

elements(d::ProductDomain) = d.domains

internal_eltype(::Type{ProductDomain{DD,S,T}}) where {DD,S,T} = S
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

# TODO: check for the efficiency of this operation
_indomain(x, d::ProductDomain) = reduce(&, map(indomain, x, elements(d)))



function show(io::IO, t::ProductDomain)
    L = nb_elements(t)
    for i in 1:L-1
        show(io, element(t, i))
        print(io, " x ")
    end
    show(io, element(t, L))
end
