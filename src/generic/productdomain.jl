
"""
```using DomainSets```

A `ProductDomain` represents the cartesian product of other domains.
"""
abstract type ProductDomain{T} <: CompositeLazyDomain{T} end

composition(d::ProductDomain) = Product()

elements(d::ProductDomain) = d.domains

==(d1::ProductDomain, d2::ProductDomain) = mapreduce(==, &, elements(d1), elements(d2))

isempty(d::ProductDomain) = any(isempty, elements(d))
isclosedset(d::ProductDomain) = all(isclosedset, elements(d))
isopenset(d::ProductDomain) = all(isopenset, elements(d))

volume(d::ProductDomain) = prod(map(volume, elements(d)))

compatibleproductdims(d1::ProductDomain, d2::ProductDomain) =
	all(map(==, map(dimension, elements(d1)), map(dimension, elements(d2))))

compatibleproduct(d1::ProductDomain, d2::ProductDomain) =
	compatibleproductdims(d1, d2) && compatible_eltype(d1, d2)

function show(io::IO, d::ProductDomain)
    L = numelements(d)
	if L <= 10
	    for i in 1:L-1
	        show(io, element(d, i))
	        print(io, " x ")
	    end
	    show(io, element(d, L))
	else
		for i in 1:5
	        show(io, element(d, i))
	        print(io, " x ")
	    end
	    print(io, "...")
		for i in L-4:L
			print(io, " x ")
	        show(io, element(d, i))
	    end
	end
end

boundary(d::ProductDomain) = _boundary(d, elements(d))
_boundary(d::ProductDomain, domains) =
	UnionDomain(ProductDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))

infimum(d::ProductDomain) = toexternalpoint(d, map(infimum, elements(d)))
supremum(d::ProductDomain) = toexternalpoint(d, map(supremum, elements(d)))

interior(d::ProductDomain) = ProductDomain(map(interior, elements(d)))
closure(d::ProductDomain) = ProductDomain(map(closure, elements(d)))


VcatDomainElement = Union{Domain{<:Number},EuclideanDomain}

ProductDomain(domains...) = _ProductDomain(map(Domain, domains)...)
_ProductDomain(domains...) = TupleProductDomain(domains...)
_ProductDomain(domains::VcatDomainElement...) = VcatDomain(domains...)
ProductDomain(domains::AbstractVector) = VectorProductDomain(domains)
# To create a tuple product domain, invoke ProductDomain{T}. Here, we splat
# and this may end up creating a VcatDomain instead.
ProductDomain(domains::Tuple) = ProductDomain(domains...)

ProductDomain{T}(domains...) where {T} = _TypedProductDomain(T, domains...)
_TypedProductDomain(::Type{SVector{N,T}}, domains...) where {N,T} = VcatDomain{N,T}(domains...)
_TypedProductDomain(::Type{T}, domains...) where {T<:Vector} = VectorProductDomain{T}(domains...)
_TypedProductDomain(::Type{T}, domains...) where {T<:Tuple} = TupleProductDomain{T}(domains...)

cross(x::Domain...) = cartesianproduct(x...)

# One can use the cartesianproduct routine to create product domains
cartesianproduct(domains...) = ProductDomain(expand(ProductDomain, domains...)...)

^(d::Domain, n::Int) = cartesianproduct(d, n)

similardomain(d::ProductDomain, ::Type{T}) where {T} = ProductDomain{T}(elements(d))

canonicaldomain(d::ProductDomain) = ProductDomain(map(canonicaldomain, elements(d)))

tocanonical(d::ProductDomain) = ProductMap(map(tocanonical, elements(d)))
fromcanonical(d::ProductDomain) = ProductMap(map(fromcanonical, elements(d)))


"""
A `VcatDomain` concatenates the element types of its member domains in a single
static vector.
"""
struct VcatDomain{N,T,DIM,DD} <: ProductDomain{SVector{N,T}}
	domains	::	DD
end

VcatDomain(domains::Union{Vector,Tuple}) = VcatDomain(domains...)
function VcatDomain(domains...)
	T = numtype(domains...)
	N = sum(map(dimension, domains))
	VcatDomain{N,T}(domains...)
end

VcatDomain{N,T}(domains::Union{AbstractVector,Tuple}) where {N,T} = VcatDomain{N,T}(domains...)
function VcatDomain{N,T}(domains...) where {N,T}
	DIM = map(dimension,domains)
	VcatDomain{N,T,DIM}(convert_numtype.(domains, T)...)
end

VcatDomain{N,T,DIM}(domains...) where {N,T,DIM} =
	VcatDomain{N,T,DIM,typeof(domains)}(domains)

tointernalpoint(d::VcatDomain{N,T,DIM}, x) where {N,T,DIM} =
	convert_fromcartesian(x, Val{DIM}())
toexternalpoint(d::VcatDomain{N,T,DIM}, y) where {N,T,DIM} =
	convert_tocartesian(y, Val{DIM}())

_boundary(d::VcatDomain, domains) =
	UnionDomain(ProductDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))



"""
A `VectorProductDomain` is a product domain of arbitrary dimension where the
element type is a vector, and all member domains have the same element type.
"""
struct VectorProductDomain{V<:AbstractVector,DD<:AbstractVector} <: ProductDomain{V}
	domains	::	DD

	function VectorProductDomain{V,DD}(domains::DD) where {V,DD}
		@assert eltype(eltype(domains)) == eltype(V)
		new(domains)
	end
end

VectorProductDomain(domains::AbstractVector) =
	VectorProductDomain{Vector{eltype(eltype(domains))}}(domains)

VectorProductDomain{V}(domains::AbstractVector{<:Domain{T}}) where {T,V<:AbstractVector{T}} =
	VectorProductDomain{V,typeof(domains)}(domains)
function VectorProductDomain{V}(domains::AbstractVector) where {T,V<:AbstractVector{T}}
	Tdomains = convert.(Domain{T}, domains)
	VectorProductDomain{V}(Tdomains)
end

# Convenience: allow constructor to be called with multiple arguments, or with
# a container that is not a vector
VectorProductDomain(domains::Domain...) = VectorProductDomain(domains)
VectorProductDomain(domains) = VectorProductDomain(collect(domains))
VectorProductDomain{V}(domains::Domain...) where {V} = VectorProductDomain{V}(domains)
VectorProductDomain{V}(domains) where {V} = VectorProductDomain{V}(collect(domains))

# the dimension equals the number of composite elements
dimension(d::VectorProductDomain) = numelements(d)

tointernalpoint(d::VectorProductDomain, x) =
	(@assert length(x) == dimension(d); x)
toexternalpoint(d::VectorProductDomain, y) =
	(@assert length(y) == dimension(d); y)


_boundary(d::VectorProductDomain, domains) =
	UnionDomain([VectorProductDomain([domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...]) for i in 1:length(domains)])



"""
A `TupleProductDomain` is a product domain that concatenates the elements of
its member domains in a tuple.
"""
struct TupleProductDomain{T,DD} <: ProductDomain{T}
	domains	::	DD
end

TupleProductDomain(domains::Vector) = TupleProductDomain(domains...)
TupleProductDomain(domains::Domain...) = TupleProductDomain(domains)
TupleProductDomain(domains...) = TupleProductDomain(map(Domain, domains)...)
function TupleProductDomain(domains::Tuple)
	T = Tuple{map(eltype, domains)...}
	TupleProductDomain{T}(domains)
end

TupleProductDomain{T}(domains::Vector) where {T} = TupleProductDomain{T}(domains...)
TupleProductDomain{T}(domains...) where {T} = TupleProductDomain{T}(domains)
function TupleProductDomain{T}(domains::Tuple) where {T <: Tuple}
	Tdomains = map((t,d) -> convert(Domain{t},d), tuple(T.parameters...), domains)
	TupleProductDomain{T,typeof(Tdomains)}(Tdomains)
end

_boundary(d::TupleProductDomain, domains) =
	UnionDomain(TupleProductDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))
