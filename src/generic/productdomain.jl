# Routines for a cartesian product of domains

###################
# Helper functions
###################

"""
Convert a vector from a cartesian format to a nested tuple according to the
given dimensions.

For example:
`convert_fromcartesian([1,2,3,4,5], Val{(2,2,1)}()) -> ([1,2],[3,4],5)`
"""
@generated function convert_fromcartesian(x::AbstractVector, ::Val{DIM}) where {DIM}
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
# An alternative is to use "reduce(vcat, x)" (see Julia issue #21672) but the
# generated function is more efficient because the compiler knows the dimensions.


#######################
# Main type definition
#######################

"""
```using DomainSets```

A `ProductDomain` represents the cartesian product of other domains.
"""
abstract type ProductDomain{T} <: LazyDomain{T} end

composition(d::ProductDomain) = Product()

elements(d::ProductDomain) = d.domains

if VERSION >= v"1.2"
	==(d1::ProductDomain, d2::ProductDomain) = mapreduce(==, &, elements(d1), elements(d2))
else
	==(d1::ProductDomain, d2::ProductDomain) = reduce(&, map(==, elements(d1), elements(d2)))
end

isempty(d::ProductDomain) = any(isempty, elements(d))
isclosedset(d::ProductDomain) = all(isclosed, elements(d))
isopenset(d::ProductDomain) = all(isopen, elements(d))

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

ProductDomain(domains) = ProductDomain(domains...)
ProductDomain(domains::VcatDomainElement...) = VcatDomain(domains...)
ProductDomain(domains...) = TupleProductDomain(domains...)
ProductDomain(domains::Vector) = VectorProductDomain(domains)

ProductDomain{T}(args...) where {T <: SVector} = VcatDomain{T}(args...)
ProductDomain{T}(args...) where {T <: Tuple} = TupleProductDomain{T}(args...)
ProductDomain{T}(args...) where {T <: Vector} = VectorProductDomain{T}(args...)

cross(x::Domain...) = cartesianproduct(x...)

# One can use the cartesianproduct routine to create product domains
cartesianproduct(domains::Domain...) = ProductDomain(expand(ProductDomain, domains...)...)

^(d::Domain, n::Int) = cartesianproduct(d, n)



"""
A `VcatDomain` concatenates the element types of its member domains in a single
static vector.
"""
struct VcatDomain{N,T,DIM,DD} <: ProductDomain{SVector{N,T}}
	domains	::	DD
end

VcatDomain(domains...) = VcatDomain(domains)

"Convert the numeric type of a domain."
convert_numtype(::Type{T}, d::Domain{T}) where {T} = d
convert_numtype(::Type{T}, d::Domain{SVector{N,T}}) where {N,T} = d
convert_numtype(::Type{T}, d::Domain{SVector{N,S}}) where {N,T,S} =
	convert(Domain{SVector{N,T}}, d)
convert_numtype(::Type{T}, d::Domain{S}) where {S<:Number,T<:Number} =
	convert(Domain{T}, d)

function VcatDomain(domains)
	T = mapreduce(numtype, promote_type, domains)
	N = sum(map(dimension, domains))
	VcatDomain{N,T}(domains)
end

function VcatDomain{N,T}(domains) where {N,T}
	DIM = map(dimension,domains)
	Tdomains = map(d->convert_numtype(T, d), domains)
	VcatDomain{N,T,DIM,typeof(Tdomains)}(Tdomains)
end

tointernalpoint(d::VcatDomain{N,T,DIM}, x) where {N,T,DIM} =
	convert_fromcartesian(x, Val{DIM}())
toexternalpoint(d::VcatDomain{N,T,DIM}, y) where {N,T,DIM} =
	convert_tocartesian(y, Val{DIM}())

convert(::Type{EuclideanDomain{N,T}}, d::VcatDomain{N,T}) where {N,T} = d
convert(::Type{EuclideanDomain{N,T}}, d::VcatDomain{N,S}) where {N,T,S} =
	VcatDomain{N,T}(elements(d))

_boundary(d::VcatDomain, domains) =
	UnionDomain(VcatDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))

"""
A `VectorProductDomain` is a product domain of arbitrary dimension with element
type `Vector{T}`, with all member domains having element type `T`.
"""
struct VectorProductDomain{T,D} <: ProductDomain{Vector{T}}
	domains	::	Vector{D}
end

VectorProductDomain(domains::Domain...) =
	VectorProductDomain{numtype(domains...)}(collect(domains))

function VectorProductDomain(domains::Vector)
	T = mapreduce(numtype, promote_type, domains)
	VectorProductDomain{T}(domains)
end

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

convert(::Type{VectorDomain{T}}, d::VectorProductDomain{T}) where {T} = d
convert(::Type{VectorDomain{T}}, d::VectorProductDomain{S}) where {S,T} =
	VectorProductDomain{T}(elements(d))

_boundary(d::VectorProductDomain, domains) =
	UnionDomain([VectorProductDomain([domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...]) for i in 1:length(domains)])


"""
A `TupleProductDomain` is a product domain that concatenates the elements of
its member domains in a tuple.
"""
struct TupleProductDomain{T,DD} <: ProductDomain{T}
	domains	::	DD
end

TupleProductDomain(domains::Domain...) = TupleProductDomain(domains)
TupleProductDomain(domains...) = TupleProductDomain(map(Domain, domains)...)
function TupleProductDomain(domains)
	T = Tuple{map(eltype, domains)...}
	TupleProductDomain{T}(domains)
end

TupleProductDomain{T}(domains...) where {T} = TupleProductDomain{T}(domains)
function TupleProductDomain{T}(domains) where {T<:Tuple}
	Tdomains = map((t,d) -> convert(Domain{t},d), tuple(T.parameters...), domains)
	TupleProductDomain{T,typeof(Tdomains)}(Tdomains)
end

convert(::Type{Domain{T}}, d::TupleProductDomain{T}) where {T} = d
convert(::Type{Domain{T}}, d::TupleProductDomain{S}) where {S,T} =
	TupleProductDomain{T}(elements(d))

_boundary(d::TupleProductDomain, domains) =
	UnionDomain(TupleProductDomain(domains[1:i-1]..., boundary(domains[i]), domains[i+1:end]...) for i in 1:length(domains))
