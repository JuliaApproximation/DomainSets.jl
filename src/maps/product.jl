
abstract type AbstractProductMap{T} <: Map{T} end

"""
A product map is diagonal and acts on each of the components of x separately:
`y = f(x)` becomes `y_i = f_i(x_i)`.
"""
struct ProductMap{T,MAPS} <: AbstractProductMap{T}
    # maps has an indexable and iterable type, for example a tuple of maps
    maps    ::  MAPS
end

function ProductMap(maps...)
    MAPS = typeof(maps)
    T = Tuple{map(domaintype, maps)...}
    ProductMap{T,MAPS}(maps)
end

elements(m::ProductMap) = m.maps

applymap(m::ProductMap, x) = map(applymap, elements(m), x)

isreal(m::ProductMap) = all(map(isreal, elements(m)))

tensorproduct(map1::Map, map2::Map) = ProductMap(map1, map2)
tensorproduct(map1::ProductMap, map2::Map) = ProductMap(elements(map1)..., map2)
tensorproduct(map1::Map, map2::ProductMap) = ProductMap(map1, elements(map2)...)
tensorproduct(map1::ProductMap, map2::ProductMap) = ProductMap(elements(map1)..., elements(map2)...)

for op in (:inv, :leftinv, :rightinv, :jacobian)
    @eval $op(m::ProductMap) = ProductMap(map($op, elements(m))...)
end
