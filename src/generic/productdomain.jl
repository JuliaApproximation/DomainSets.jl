# productdomain.jl

########################################
### A tensor product of Domains
########################################

"""
A ProductDomain represents the tensor product of other domains.

struct ProductDomain{TD,N} <: Domain{N}

Parameters:
- TD is a tuple of (domain) types
- N is the total dimension of this domain
"""
struct ProductDomain{TD,N} <: Domain{N}
	domains	::	TD

	# Inner constructor to verify that N is correct
	function ProductDomain{TD,N}(domains) where {TD,N}
		@assert sum(map(ndims, domains)) == N
		new{TD,N}(domains)
	end
end

elements(d::ProductDomain) = d.domains

function ProductDomain(domains...)
    TD = typeof(domains)
    N = sum(map(ndims, domains))
    ProductDomain{TD,N}(domains)
end

tensorproduct(d::Domain) = d
tensorproduct(d::Domain, n::Int) = tensorproduct([d for i=1:n]...)
tensorproduct(d::Domain...) =
    ProductDomain(flatten(ProductDomain, d...)...)

⊗ = tensorproduct

^(d::Domain, n::Int) = tensorproduct(d, n)

indomain(x, t::ProductDomain) = indomain(x, elements(t)...)

indomain(x::SVector{2}, d1::Domain{1}, d2::Domain{1}) =
	indomain(x[1], d1) && indomain(x[2], d2)

indomain(x::SVector{3}, d1::Domain{1}, d2::Domain{1}, d3::Domain{1}) =
	indomain(x[1], d1) && indomain(x[2], d2) && indomain(x[3], d3)

indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{1}, d3::Domain{1}, d4::Domain{1}) =
	indomain(x[1], d1) && indomain(x[2], d2) && indomain(x[3], d3) && indomain(x[4], d4)

indomain(x::SVector{3}, d1::Domain{1}, d2::Domain{2}) =
	indomain(x[1], d1) && indomain(SVector(x[2],x[3]), d2)

indomain(x::SVector{3}, d1::Domain{2}, d2::Domain{1}) =
	indomain(SVector(x[1],x[2]), d1) && indomain(x[3], d2)

indomain(x::SVector{4}, d1::Domain{2}, d2::Domain{2}) =
	indomain(SVector(x[1],x[2]), d1) && indomain(SVector(x[3],x[4]), d2)

indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{3}) =
	indomain(x[1], d1) && indomain(SVector(x[2],x[3],x[4]), d2)

indomain(x::SVector{4}, d1::Domain{3}, d2::Domain{1}) =
	indomain(SVector(x[1],x[2],x[3]), d1) && indomain(x[1], d2)

indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{1}, d3::Domain{2}) =
	indomain(x[1], d1) && indomain(x[2], d2) && indomain(SVector(x[3],x[4]), d3)

indomain(x::SVector{4}, d1::Domain{1}, d2::Domain{2}, d3::Domain{1}) =
	indomain(x[1], d1) && indomain(SVector(x[2],x[3]), d2) && indomain(x[4], d3)

indomain(x::SVector{4}, d1::Domain{2}, d2::Domain{1}, d3::Domain{1}) =
	indomain(SVector(x[1],x[2]), d1) && indomain(x[3], d2) && indomain(x[4], d3)

# TODO: make this code for indomain more general!
indomain{N}(x::SVector{N}, d::Vararg{Domain,N}) =
	reduce(&, map(indomain, x, d))

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
