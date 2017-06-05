# productdomain.jl

###############################
# A tensor product of Domains
###############################

product_eltype(domains::Domain...) = Tuple{map(eltype, domains)...}

"""
A `ProductDomain` represents the tensor product of other domains.

Parameters:
- D is the type of a tuple of domain types
- T is the eltype of the product space
"""
struct ProductDomain{D,T} <: Domain{T}
	domains	::	D

	# Inner constructor to verify that T is correct
	function ProductDomain{D,T}(domains) where {D,T}
		@assert T == product_eltype(domains...)
		new{D,T}(domains)
	end
end

elements(d::ProductDomain) = d.domains

function ProductDomain(domains...)
    D = typeof(domains)
    T = product_eltype(domains...)
    ProductDomain{D,T}(domains)
end

tensorproduct(d::Domain) = d
tensorproduct(d::Domain, n::Int) = tensorproduct([d for i=1:n]...)
tensorproduct(d::Domain...) =
    ProductDomain(flatten(ProductDomain, d...)...)

^(d::Domain, n::Int) = tensorproduct(d, n)

indomain(x, d::ProductDomain) = reduce(&, map(indomain, x, d))

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
#
# # TODO: make this code for indomain more general!
# indomain{N}(x::SVector{N}, d::Vararg{Domain,N}) =
# 	reduce(&, map(indomain, x, d))

# TODO: provide implementation of indomain for tensorproductgrids

boundingbox(d::ProductDomain) = tensorproduct(map(boundingbox, elements(d))...)

function show(io::IO, t::ProductDomain)
    L = nb_elements(t)
    for i in 1:L-1
        show(io, element(t, i))
        print(io, " x ")
    end
    show(io, element(t, L))
end


(*)(d::ProductDomain, x::Number) = tensorproduct([domain*x for domain in elements(d)]...)
