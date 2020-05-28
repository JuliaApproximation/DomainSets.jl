
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

function show(io::IO, d::ProductDomain)
    L = numelements(d)
    for i in 1:L-1
        show(io, element(d, i))
        print(io, " x ")
    end
    show(io, element(d, L))
end

boundary(d::ProductDomain) = _boundary(d, elements(d))

VcatDomainElement = Union{Domain{<:Number},EuclideanDomain}

ProductDomain(domains...) = TupleProductDomain(domains...)
ProductDomain(domains::VcatDomainElement...) = VcatDomain(domains...)
ProductDomain(domains::Vector) = VectorProductDomain(domains)

ProductDomain{SVector{N,T}}(domains...) where {N,T} = VcatDomain{N,T}(domains...)
ProductDomain{Vector{T}}(domains...) where {T} = VectorProductDomain{T}(domains...)
ProductDomain{T}(domains...) where {T <: Tuple} = TupleProductDomain{T}(domains...)

cross(x::Domain...) = cartesianproduct(x...)

# One can use the cartesianproduct routine to create product domains
cartesianproduct(domains::Domain...) = ProductDomain(expand(ProductDomain, domains...)...)

^(d::Domain, n::Int) = cartesianproduct(d, n)

# Convert from any product domain to any other product domain. This means that the
# constructors should work for elements(d) of any product domain, i.e., for a tuple
# or a vector of domains.
convert(::Type{Domain{T}}, d::ProductDomain{T}) where {T} = d
convert(::Type{Domain{T}}, d::ProductDomain) where {T} = ProductDomain{T}(elements(d))

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

VcatDomain{N,T}(domains::Union{Vector,Tuple}) where {N,T} = VcatDomain{N,T}(domains...)
function VcatDomain{N,T}(domains...) where {N,T}
	DIM = map(dimension,domains)
	VcatDomain{N,T,DIM}(map(d->convert_numtype(d, T), domains)...)
end

VcatDomain{N,T,DIM}(domains...) where {N,T,DIM} =
	VcatDomain{N,T,DIM,typeof(domains)}(domains)

tointernalpoint(d::VcatDomain{N,T,DIM}, x) where {N,T,DIM} =
	convert_fromcartesian(x, Val{DIM}())
toexternalpoint(d::VcatDomain{N,T,DIM}, y) where {N,T,DIM} =
	convert_tocartesian(y, Val{DIM}())

_boundary(d::VcatDomain, domains) =
	UnionDomain(VcatDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))

"""
A `VectorProductDomain` is a product domain of arbitrary dimension with element
type `Vector{T}`, with all member domains having element type `T`.
"""
struct VectorProductDomain{T,D} <: ProductDomain{Vector{T}}
	domains	::	Vector{D}
end

VectorProductDomain(domains::Domain...) = VectorProductDomain(domains)
VectorProductDomain(domains) = VectorProductDomain(collect(domains))
function VectorProductDomain(domains::Vector)
	T = mapreduce(numtype, promote_type, domains)
	VectorProductDomain{T}(domains)
end

VectorProductDomain{T}(domains::Domain...) where {T} = VectorProductDomain{T}(domains)
VectorProductDomain{T}(domains) where {T} = VectorProductDomain{T}(collect(domains))
function VectorProductDomain{T}(domains::Vector) where {T}
	Tdomains = convert.(Domain{T}, domains)
	VectorProductDomain{T,eltype(Tdomains)}(Tdomains)
end

dimension(d::VectorProductDomain) = numelements(d)

tointernalpoint(d::VectorProductDomain, x) =
	(@assert length(x) == dimension(d); x)
toexternalpoint(d::VectorProductDomain, y) =
	(@assert length(y) == dimension(d); y)

infimum(d::ProductDomain) = toexternalpoint(d, map(infimum, elements(d)))
supremum(d::ProductDomain) = toexternalpoint(d, map(supremum, elements(d)))

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
TupleProductDomain(domain::Domain) = TupleProductDomain((domain,))
function TupleProductDomain(domains::Tuple)
	T = Tuple{map(eltype, domains)...}
	TupleProductDomain{T}(domains)
end

TupleProductDomain{T}(domains::Vector) where {T} = TupleProductDomain{T}(domains...)
TupleProductDomain{T}(domains...) where {T} = TupleProductDomain{T}(domains)
function TupleProductDomain{T}(domains) where {T <: Tuple}
	Tdomains = map((t,d) -> convert(Domain{t},d), tuple(T.parameters...), domains)
	TupleProductDomain{T,typeof(Tdomains)}(Tdomains)
end

_boundary(d::TupleProductDomain, domains) =
	UnionDomain(TupleProductDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))
