
abstract type AbstractProductMap{T} <: CompositeLazyMap{T} end

"""
A product map is diagonal and acts on each of the components of x separately:
`y = f(x)` becomes `y_i = f_i(x_i)`.
"""
struct ProductMap{T,MAPS} <: AbstractProductMap{T}
    maps    ::  MAPS
end

function ProductMap(maps...)
    T = Tuple{map(domaintype, maps)...}
    ProductMap{T}(maps...)
end

ProductMap{T}(maps...) where {T} = ProductMap{T,typeof(maps)}(maps)

# TODO: make proper conversion
convert(::Type{Map{T}}, m::ProductMap{T}) where {T} = m
convert(::Type{Map{T}}, m::ProductMap{S}) where {S,T} = ProductMap{T}(m.maps...)

applymap(m::ProductMap, x) = map(applymap, elements(m), x)

tensorproduct(map1::Map, map2::Map) = ProductMap(map1, map2)
tensorproduct(map1::ProductMap, map2::Map) = ProductMap(elements(map1)..., map2)
tensorproduct(map1::Map, map2::ProductMap) = ProductMap(map1, elements(map2)...)
tensorproduct(map1::ProductMap, map2::ProductMap) = ProductMap(elements(map1)..., elements(map2)...)

for op in (:inv, :leftinverse, :rightinverse, :jacobian)
    @eval $op(m::ProductMap) = ProductMap(map($op, elements(m))...)
end

# isconstant(m::ProductMap) = mapreduce(isconstant, &, elements(m))
# islinear(m::ProductMap) = mapreduce(islinear, &, elements(m))
# isaffine(m::ProductMap) = mapreduce(isaffine, &, elements(m))

size(m::ProductMap) = reduce((x,y) -> (x[1]+y[1],x[2]+y[2]), map(size,elements(m)))

==(m1::ProductMap, m2::ProductMap) = all(map(isequal, elements(m1), elements(m2)))
