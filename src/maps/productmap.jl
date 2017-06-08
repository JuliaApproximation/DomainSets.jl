# productmap.jl

"""
A product map is diagonal and acts on each of the components of x separately:
y = f(x) becomes y_i = f_i(x_i)
"""
struct ProductMap{MAPS,T,S} <: AbstractMap{T,S}
    # maps has an indexable and iterable type, for example a tuple of maps
    maps    ::  MAPS
end

function ProductMap(maps...)
    MAPS = typeof(maps)
    S = map(domaintype, maps)
    T = map(rangetype, maps)
    ProductMap{MAPS,T,S}(maps)
end

elements(pm::ProductMap) = pm.maps

applymap(pm::ProductMap, x) = map(apply_map, elements(pm), x)

tensorproduct(map1::AbstractMap, map2::AbstractMap) = ProductMap(map1, map2)
tensorproduct(map1::ProductMap, map2::AbstractMap) = ProductMap(elements(map1)..., map2)
tensorproduct(map1::AbstractMap, map2::ProductMap) = ProductMap(map1, elements(map2)...)
tensorproduct(map1::ProductMap, map2::ProductMap) = ProductMap(elements(map1)..., elements(map2)...)

inv(pm::ProductMap) = ProductMap(map(inv, elements(pm)))

for op in (:is_linear, :isreal)
    @eval $op(dmap::ProductMap) = reduce(&, map($op, elements(dmap)))
end


# TODO: implement jacobian
# jacobian(map::ProductMap, x) =
