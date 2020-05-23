
"""
A product map is diagonal and acts on each of the components of x separately:
`y = f(x)` becomes `y_i = f_i(x_i)`.
"""
abstract type ProductMap{T} <: CompositeLazyMap{T} end

VcatMapElement = Union{Map{<:SVector},Map{<:Number}}

ProductMap(maps::VcatMapElement...) = VcatProductMap(maps...)
ProductMap{T}(maps::VcatMapElement...) where {N,S,T <: SVector{N,S}} = VcatProductMap{S,N}(maps)

elements(m::ProductMap) = m.maps

convert(::Type{Map{T}}, m::ProductMap{T}) where {T} = m
convert(::Type{Map{T}}, m::ProductMap{S}) where {S,T} = ProductMap{T}(m.maps...)

tointernalpoint(m::ProductMap, x) = x
toexternalpoint(m::ProductMap, y) = y

applymap(m::ProductMap, x) =
	toexternalpoint(m, map(applymap, elements(m), tointernalpoint(m, x)))

tensorproduct(map1::Map, map2::Map) = ProductMap(map1, map2)
tensorproduct(map1::ProductMap, map2::Map) = ProductMap(elements(map1)..., map2)
tensorproduct(map1::Map, map2::ProductMap) = ProductMap(map1, elements(map2)...)
tensorproduct(map1::ProductMap, map2::ProductMap) = ProductMap(elements(map1)..., elements(map2)...)

for op in (:inv, :leftinverse, :rightinverse)
    @eval $op(m::ProductMap) = ProductMap(map($op, elements(m))...)
	@eval $op(m::ProductMap, x) = toexternalpoint(m, map($op, elements(m), tointernalpoint(m, x)))
end

size(m::ProductMap) = reduce((x,y) -> (x[1]+y[1],x[2]+y[2]), map(size,elements(m)))

==(m1::ProductMap, m2::ProductMap) = all(map(isequal, elements(m1), elements(m2)))


"""
A `VcatProductMap` is a product map with domain and codomain vectors
concatenated (`vcat`) into a single vector.
"""
struct VcatProductMap{T,N,DIM,MAPS} <: ProductMap{SVector{N,T}}
    maps    ::  MAPS
end

VcatProductMap(maps...) = VcatProductMap(maps)
function VcatProductMap(maps)
	T = mapreduce(numtype, promote_type, maps)
	VcatProductMap{T}(maps)
end

function VcatProductMap{T}(maps) where {T}
	M,N = reduce((x,y) -> (x[1]+y[1],x[2]+y[2]), map(size,maps))
	@assert M==N
	VcatProductMap{T,N}(maps)
end

mapdim(map) = size(map,2)

function VcatProductMap{T,N}(maps) where {T,N}
	DIM = map(mapdim,maps)
	VcatProductMap{T,N,DIM}(maps)
end

VcatProductMap{T,N,DIM}(maps) where {T,N,DIM} = VcatProductMap{T,N,DIM,typeof(maps)}(maps)

size(m::VcatProductMap{T,N}) where {T,N} = (N,N)

tointernalpoint(m::VcatProductMap{T,N,DIM}, x) where {T,N,DIM} =
	convert_fromcartesian(x, Val{DIM}())
toexternalpoint(m::VcatProductMap{T,N,DIM}, y) where {T,N,DIM} =
	convert_tocartesian(y, Val{DIM}())

convert(::Type{Map{SVector{N,T}}}, m::VcatProductMap{T,M,N}) where {M,N,T} = m
convert(::Type{Map{SVector{N,T}}}, m::VcatProductMap) where {N,T} = VcatProductMap{T,N}(m.maps)

isconstant(m::VcatProductMap) = mapreduce(isconstant, &, elements(m))
islinear(m::VcatProductMap) = mapreduce(islinear, &, elements(m))
isaffine(m::VcatProductMap) = mapreduce(isaffine, &, elements(m))

matsize(A::AbstractArray) = size(A)
matsize(A::Number) = (1,1)

function toexternalmatrix(m::VcatProductMap{T,N}, matrices) where {T,N}
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

matrix(m::VcatProductMap) = toexternalmatrix(m, map(matrix, elements(m)))

vector(m::VcatProductMap) = toexternalpoint(m, map(vector, elements(m)))

constant(m::VcatProductMap) = toexternalpoint(m, map(constant, elements(m)))

jacobian(m::VcatProductMap, x) =
	toexternalmatrix(m, map(jacobian, elements(m), tointernalpoint(m, x)))

jacobian(m::VcatProductMap) = LazyJacobian(m)
