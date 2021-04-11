
"""
A product map is diagonal and acts on each of the components of x separately:
`y = f(x)` becomes `y_i = f_i(x_i)`.
"""
abstract type ProductMap{T} <: CompositeLazyMap{T} end

elements(m::ProductMap) = m.maps

VcatMapElement = Union{Map{<:SVector},Map{<:Number}}

ProductMap(maps::Tuple) = ProductMap(maps...)
ProductMap(maps::SVector) = ProductMap(maps...)
ProductMap(maps...) = TupleProductMap(maps...)
ProductMap(maps::VcatMapElement...) = VcatMap(maps...)
ProductMap(maps::AbstractVector) = VectorProductMap(maps)

ProductMap{T}(maps...) where {T} = _TypedProductMap(T, maps...)
_TypedProductMap(::Type{T}, maps...) where {T<:Tuple} = TupleProductMap(maps...)
_TypedProductMap(::Type{SVector{N,T}}, maps...) where {N,T} = VcatMap{N,T}(maps...)
_TypedProductMap(::Type{T}, maps...) where {T<:AbstractVector} = VectorProductMap{T}(maps...)


isconstant(m::ProductMap) = mapreduce(isconstant, &, elements(m))
islinear(m::ProductMap) = mapreduce(islinear, &, elements(m))
isaffine(m::ProductMap) = mapreduce(isaffine, &, elements(m))

matrix(m::ProductMap) = toexternalmatrix(m, map(matrix, elements(m)))
vector(m::ProductMap) = toexternalpoint(m, map(vector, elements(m)))
constant(m::ProductMap) = toexternalpoint(m, map(constant, elements(m)))

jacobian(m::ProductMap, x) =
	toexternalmatrix(m, map(jacobian, elements(m), tointernalpoint(m, x)))
jacobian(m::ProductMap) = LazyJacobian(m)

similarmap(m::ProductMap, ::Type{T}) where {T} = ProductMap{T}(elements(m))

tointernalpoint(m::ProductMap, x) = x
toexternalpoint(m::ProductMap, y) = y

applymap(m::ProductMap, x) =
	toexternalpoint(m, map(applymap, elements(m), tointernalpoint(m, x)))

productmap(map1, map2) = productmap1(map1, map2)
productmap1(map1, map2) = productmap2(map1, map2)
productmap2(map1, map2) = ProductMap(map1, map2)
productmap(map1::ProductMap, map2::ProductMap) =
	ProductMap(elements(map1)..., elements(map2)...)
productmap1(map1::ProductMap, map2) = ProductMap(elements(map1)..., map2)
productmap2(map1, map2::ProductMap) = ProductMap(map1, elements(map2)...)

for op in (:inv, :leftinverse, :rightinverse)
    @eval $op(m::ProductMap) = ProductMap(map($op, elements(m)))
	@eval $op(m::ProductMap, x) = toexternalpoint(m, map($op, elements(m), tointernalpoint(m, x)))
end

size(m::ProductMap) = reduce((x,y) -> (x[1]+y[1],x[2]+y[2]), map(size,elements(m)))

==(m1::ProductMap, m2::ProductMap) = all(map(isequal, elements(m1), elements(m2)))

"""
A `VcatMap` is a product map with domain and codomain vectors
concatenated (`vcat`) into a single vector.
"""
struct VcatMap{N,T,DIM,MAPS} <: ProductMap{SVector{N,T}}
    maps    ::  MAPS
end

VcatMap(maps::Union{Tuple,Vector}) = VcatMap(maps...)
function VcatMap(maps...)
	T = numtype(maps...)
	M,N = reduce((x,y) -> (x[1]+y[1],x[2]+y[2]), map(size,maps))
	@assert M==N
	VcatMap{N,T}(maps...)
end

mapdim(map) = size(map,2)

VcatMap{N,T}(maps::Union{Tuple,Vector}) where {N,T} = VcatMap{N,T}(maps...)
function VcatMap{N,T}(maps...) where {N,T}
	DIM = map(mapdim,maps)
	VcatMap{N,T,DIM}(maps...)
end

VcatMap{N,T,DIM}(maps...) where {N,T,DIM} = VcatMap{N,T,DIM,typeof(maps)}(maps)

size(m::VcatMap{N}) where {N} = (N,N)

tointernalpoint(m::VcatMap{N,T,DIM}, x) where {N,T,DIM} =
	convert_fromcartesian(x, Val{DIM}())
toexternalpoint(m::VcatMap{N,T,DIM}, y) where {N,T,DIM} =
	convert_tocartesian(y, Val{DIM}())

matsize(A::AbstractArray) = size(A)
matsize(A::Number) = (1,1)

# The Jacobian is block-diagonal
function toexternalmatrix(m::VcatMap{N,T}, matrices) where {N,T}
	A = zeros(T, N, N)
	l = 0
	for el in matrices
		m,n = matsize(el)
		@assert m==n
		A[l+1:l+m,l+1:l+n] .= el
		l += n
	end
	SMatrix{N,N}(A)
end


"""
A `VectorProductMap` is a product map where all components are univariate maps,
with inputs and outputs collected into a `Vector`.
"""
struct VectorProductMap{T<:AbstractVector,M} <: ProductMap{T}
    maps    ::  Vector{M}
end

VectorProductMap(maps::Map...) = VectorProductMap(maps)
VectorProductMap(maps) = VectorProductMap(collect(maps))
function VectorProductMap(maps::Vector)
	T = mapreduce(numtype, promote_type, maps)
	VectorProductMap{Vector{T}}(maps)
end

VectorProductMap{T}(maps::Map...) where {T} = VectorProductMap{T}(maps)
VectorProductMap{T}(maps) where {T} = VectorProductMap{T}(collect(maps))
function VectorProductMap{T}(maps::Vector) where {T}
	Tmaps = convert.(Map{eltype(T)}, maps)
	VectorProductMap{T,eltype(Tmaps)}(Tmaps)
end

# the Jacobian is a diagonal matrix
toexternalmatrix(m::VectorProductMap, matrices) = Diagonal(matrices)


"""
A `TupleProductMap` is a product map with all components collected in a tuple.
There is no vector-valued function associated with this map.
"""
struct TupleProductMap{T,MM} <: ProductMap{T}
    maps    ::  MM
end

TupleProductMap(maps::Vector) = TupleProductMap(maps...)
TupleProductMap(maps::Map...) = TupleProductMap(maps)
TupleProductMap(maps...) = TupleProductMap(map(Map, maps)...)
TupleProductMap(map::Map) = TupleProductMap((map,))
function TupleProductMap(maps::Tuple)
	T = Tuple{map(eltype, maps)...}
	TupleProductMap{T}(maps)
end

TupleProductMap{T}(maps::Vector) where {T} = TupleProductMap{T}(maps...)
TupleProductMap{T}(maps::Map...) where {T} = TupleProductMap{T}(maps)
TupleProductMap{T}(maps...) where {T} = TupleProductMap{T}(maps)
function TupleProductMap{T}(maps) where {T <: Tuple}
	Tmaps = map((t,d) -> convert(Map{t},d), tuple(T.parameters...), maps)
	TupleProductMap{T,typeof(Tmaps)}(Tmaps)
end
